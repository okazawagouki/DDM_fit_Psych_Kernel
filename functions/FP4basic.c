#include <math.h>
#include <string.h>
#include <immintrin.h>
#include "mex.h"
#include "matrix.h"

#define CACHE_LINE_SIZE 64

#ifdef __GNUC__
    #define LIKELY(x) __builtin_expect(!!(x), 1)
    #define UNLIKELY(x) __builtin_expect(!!(x), 0)
#endif

void fpe_slover_cc (double *v_lut, double *v_init, int vlen, double dmesh, int stepnum, double dt, double *mu, double *sigma, double alpha,
                        double lb_margin, double ub_margin, double *lb_change, double *ub_change,
                        bool static_mu, bool static_noise, int num_output,
                        double *v_final, double *prob_notcross, double *prob_cross, double *prob_gtzero, double *full_prob);

static void make_tridiag_factor (double B, double C, double alpha, double dt, double dmesh, int vlen, double *v_lut,
                                double *dl, double *d, double *du, double *buf_diag, double *buf_scale);

static void factor_coeff_matrix (double *dl, double *d, double *du, 
                            double *buf_diag, double *buf_scale, int vlen);

static inline void solve_linear_eq (double *x, double *b, double *buf_diag, double *buf_scale, int vlen);

static void get_prob_info (double *ptr_v, int *bound, int t, int vlen, int stepnum, int num_output,
        double *prob_notcross, double *prob_cross, double *prob_gtzero, double *full_prob);

/*******************
Matlab call:
    [ufinal, Pt, Ptlow, Pthigh, Pg0, Pxt] = FP4 ( xmesh, uinit, k, sigma, [lb_change ub_change], [lb_margin; ub_margin], dt, alpha);

Description:
    Propagates the probability density <uinit> on the spatial grid <xmesh> across time, using the diffusion
    parameters k (drift rate) and sigma (diffusion standard deviation). The duration of the propagation is
    defined by dt and size of [lb_change ub_change] (see below). Don't include Pxt in the output arguments
    to keep memory management efficient and fast.

Input:
    0 [vmesh]
        A column vector specifying the spatial grid for probability density propagation in space
        the values in xmesh are monotonically increasing and range from low values (lower bound - lb_margin + dx) to
        high values (upper bound + ub_margin - dx). xmesh should be regular. define it as:
        lower bound - lb_margin + dx : dx : upper bound + ub_margin - dx

    1 [vinit]
        A column vector specifying the initial probabilities on the grid it is usually a delta function, such as:
        uinit = zeros(size(xmesh));   uinit(xmesh==0) = 1;

    2 [mu]
        Either a scalar specifying the constant drift rate or a vector specifying time varying drift rate.
        In case of vector, the length should be equal to the time step of [lb_change, ub_change].

    3 [sigma]
        Either a scalar specifying the constant SD of the diffusion process or a vector specifying time varying SD.
        In case of vector, the length should be equal to the time step of [lb_change, ub_change].

    4 [lb_change, ub_change]
        A matrix consisting of two columns, each column shows how the corresponding margin changes across time.

        The number of rows in this matrix defines the number of time steps for the propagation of the probability
        density. stepnum*dt is equal to the total propagation time
        
        Make sure that the change in the bound height is valid: bound height should change so that the upper
        bound stay above zero and the lower bound stay below zero all the time, also neither of them can exceed the xmesh.

    5 [lb_margin, ub_margin]
        A column vector specifying the margin of the lower and upper boundaries.
        
        To find the probability of crossing the upper and lower boundaries separately we should reduce xinit and increase xend by
        many folds of sigma and calculate how much of the probability falls beyond the actual target boundaries.
        
        [lb_margin] and [ub_margin] define the distance of the actual boundaries from xinit and xend.
        The starting and ending point of xmesh. When lb_margin and ub_margin are zero xinit and xend correspond to the boundaries.

    6 [dt]
        The size of each time step, specifying the resolution of the temporal grid

    7 [alpha]
        Leak parameter. the inverse of leak time constant (tau). If set to 0, it becomes the diffusion process.
        default: 0 (diffusion process)


output:
    0 [ufinal]
        A column vector, final probability density at the end of the propagation

    1 [Pt]
        A column vector, it is the survivor function: the probability of NOT crossing the bounds up to each moment during the propagation

    2 [Ptlow, Pthigh]
        A matrix with two columns, the first column corresponds to the probability of crossing the lower bound
        at each moment. the second column corresponds to the probability of crossing the upper bound at each moment

    3 [Pg0]
        A column vector, the probability of staying above zero at each moment. it is useful when the diffusion process terminates before bound crossing
        I have written the program so that it includes also half of the mass in the spatial grid location that corresponds to zero

    4 [Pxt]
        A matrix with columns corresponding to the probability density at each moment during the propagation. 
        Its calculation and storage requires a huge amount of memory. don't ask for it if you are not sure that your computer can handle it or not.

*******************/

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *v_lut;
    double *v_init;
    double *mu;
    double *sigma;
    double dmesh, dt;
    double *bound_change, *bound_margin;
    double *lb_change, *ub_change;
    double *change_pos, *change_neg; 
    double ub_margin, lb_margin;
    double alpha = 0.0;

    int stepnum, vlen;
    bool static_mu, static_noise;
    bool no_output = false;

    switch (nrhs) {
        case 8:
            alpha = mxGetScalar(prhs[7]);

        case 7:
            vlen = (int)mxGetM(prhs[0]);
            stepnum = (int)mxGetM(prhs[4]);

            // Verify parameter range
            if (mxGetM(prhs[0]) != mxGetM(prhs[1]))
                mexErrMsgIdAndTxt("FP4new:ParamError", "Mesh value lut has different size of v_init.");
            if (mxGetN(prhs[4])!=2)
                mexErrMsgIdAndTxt("FP4new:ParamError", "Boundary margin and change should have 2 columns.");
            if (mxGetNumberOfElements(prhs[5])!=2)
                mexErrMsgIdAndTxt("FP4new:ParamError", "Boundary margin should be a vector with two elements.");
            if ((mxGetM(prhs[2])!=stepnum) && (mxGetM(prhs[2])!=1))
                mexErrMsgIdAndTxt("FP4new:ParamError", "Drift rate should be either scalar or the same length as bound.");
            if ((mxGetM(prhs[3])!=stepnum) && (mxGetM(prhs[3])!=1))
                mexErrMsgIdAndTxt("FP4new:ParamError", "Noise should be either scalar or the same length as bound.");

            // Get arguments from prhs
            v_lut = mxGetDoubles(prhs[0]);
            v_init = mxGetDoubles(prhs[1]);
            mu = mxGetDoubles(prhs[2]);
            sigma = mxGetDoubles(prhs[3]);
            bound_change = mxGetDoubles(prhs[4]);
            bound_margin = mxGetDoubles(prhs[5]);
            lb_margin = bound_margin[0];
            ub_margin = bound_margin[1];
            dt = mxGetScalar(prhs[6]);
            break;

        default:
            mexErrMsgIdAndTxt("FP4new:ParamError", "Illegal number of inputs");
            break;
    }

    dmesh = v_lut[1]-v_lut[0];

    // Check if the bound is valid
    lb_change = bound_change;
    ub_change = bound_change+stepnum;
    for (int i=0; i<stepnum; i++) {
        if ( lb_change[i]+lb_margin < dmesh || lb_change[i]+lb_margin+v_lut[0] > -dmesh ||
            -ub_change[i]+ub_margin < dmesh || -ub_change[i]+ub_margin > v_lut[vlen-1]-dmesh ) {
                mexErrMsgIdAndTxt("FP4new:ParamError", "Bound change crossed zero or exceed xmesh range");
        }
    }

    // Infer noise, evidence and leak mode
    static_mu = mxGetM(prhs[2])==1 ? true : false;
    static_noise = mxGetM(prhs[3])==1 ? true : false;
    
    // Handle outputs
    mxArray *Ufinal;
    mxArray *Pt;
    mxArray *Ptlh;
    mxArray *Pg0;
    mxArray *Pxt;

    double *v_final = NULL;
    double *prob_notcross = NULL;
    double *prob_gtzero = NULL;
    double *prob_cross = NULL; //[low high]
    double *full_prob = NULL;

    // Prealloc output
    switch (nlhs) {
        case 5:
            Pxt = mxCreateDoubleMatrix(vlen, stepnum, mxREAL);
            full_prob = mxGetDoubles(Pxt);

        case 4:
            Pg0 = mxCreateDoubleMatrix(stepnum+1, 1, mxREAL);
            prob_gtzero = mxGetDoubles(Pg0);

        case 3:
            Ptlh = mxCreateDoubleMatrix(stepnum+1, 2, mxREAL);  // follows C convention
            prob_cross = mxGetDoubles(Ptlh);

        case 2:
            Pt = mxCreateDoubleMatrix(stepnum+1, 1, mxREAL);
            prob_notcross = mxGetDoubles(Pt);

        case 1:
            Ufinal = mxCreateDoubleMatrix(vlen, 1, mxREAL);
            v_final = mxGetDoubles(Ufinal);
            break;

        case 0:
            mexWarnMsgIdAndTxt("FP4new:NoOutput", "No output, quit directly");
            no_output = true;
            break;

        default:
            mexErrMsgIdAndTxt("FP4new:OutputError", "Illegal number of outputs");
            break;
    }

    if (no_output) { return; }

    /***** Main computation *****/
    fpe_slover_cc(v_lut, v_init, vlen, dmesh, stepnum, dt, mu, sigma, alpha,
                    lb_margin, ub_margin, lb_change, ub_change,
                    static_mu, static_noise, nlhs,
                    v_final, prob_notcross, prob_cross, prob_gtzero, full_prob);

    // Save output
    switch (nlhs) {
        case 5:
            plhs[4] = Pxt;

        case 4:
            plhs[3] = Pg0;

        case 3:
            plhs[2] = Ptlh;

        case 2:
            plhs[1] = Pt;

        case 1:
            plhs[0] = Ufinal;
            break;
    }
}

void fpe_slover_cc (double *v_lut, double *v_init, int vlen, double dmesh, int stepnum, double dt, double *mu, double *sigma, double alpha,
                        double lb_margin, double ub_margin, double *lb_change, double *ub_change,
                        bool static_mu, bool static_noise, int num_output,
                        double *v_final, double *prob_notcross, double *prob_cross, double *prob_gtzero, double *full_prob)
{
    // Alloc spaces for solving the tridiagonal system
    // Memory alignment pattern for vinit (19 elements double array, cache line 64 byte)
    // |<-     cache line A     ->|<-     cache line B     ->|<-     cache line C     ->|<-     cache line D     ->|
    // |     (>=1 element)  1  2  | 3  4  5  6  7  8  9  10  |  11 12 13 14 15 16 17 18 | 19  (only one element)   |
    double *v_curr, *v_next;
    double *dl, *d, *du;
    double *buf_diag, *buf_scale;
    int *bound;

    int cache_num_vlen = (vlen-1)%(CACHE_LINE_SIZE/sizeof(double))==0 ? (vlen-1)/(CACHE_LINE_SIZE/sizeof(double))+1 
                            : (vlen-1)/(CACHE_LINE_SIZE/sizeof(double))+2;
    int cache_num_bound = 2*(stepnum+1)%(CACHE_LINE_SIZE/sizeof(int))==0 ? 2*(stepnum+1)/(CACHE_LINE_SIZE/sizeof(int)) 
                            : 2*(stepnum+1)/(CACHE_LINE_SIZE/sizeof(int))+1;
    int cache_num_d = vlen%(CACHE_LINE_SIZE/sizeof(double))==0 ? vlen/(CACHE_LINE_SIZE/sizeof(double))
                            : vlen/(CACHE_LINE_SIZE/sizeof(double))+1;
    int cache_num_buf = 2*(vlen-1)%(CACHE_LINE_SIZE/sizeof(double))==0 ? 2*(vlen-1)/(CACHE_LINE_SIZE/sizeof(double))+1 
                            : 2*(vlen-1)/(CACHE_LINE_SIZE/sizeof(double))+2;
    int vlen_size = cache_num_vlen*CACHE_LINE_SIZE;
    int bound_size = cache_num_bound*CACHE_LINE_SIZE;
    int buf_size = cache_num_buf*CACHE_LINE_SIZE;
    int d_size = cache_num_d*CACHE_LINE_SIZE;

#ifdef _MSC_VER         
    v_curr = (double *)_aligned_malloc(vlen_size, CACHE_LINE_SIZE);
    v_next = (double *)_aligned_malloc(vlen_size, CACHE_LINE_SIZE);
    buf_scale = (double *)_aligned_malloc(vlen_size, CACHE_LINE_SIZE);
    dl = (double *)_aligned_malloc(d_size, CACHE_LINE_SIZE);
    d = (double *)_aligned_malloc(d_size, CACHE_LINE_SIZE);
    du = (double *)_aligned_malloc(d_size, CACHE_LINE_SIZE);
    buf_diag = (double *)_aligned_malloc(buf_size, CACHE_LINE_SIZE);
    bound = (int *)_aligned_malloc(bound_size, CACHE_LINE_SIZE);
#else                   
    v_curr = (double *)aligned_alloc(CACHE_LINE_SIZE, vlen_size);
    v_next = (double *)aligned_alloc(CACHE_LINE_SIZE, vlen_size);
    buf_scale = (double *)aligned_alloc(CACHE_LINE_SIZE, vlen_size);
    dl = (double *)aligned_alloc(CACHE_LINE_SIZE, d_size);
    d = (double *)aligned_alloc(CACHE_LINE_SIZE, d_size);
    du = (double *)aligned_alloc(CACHE_LINE_SIZE, d_size);
    buf_diag = (double *)aligned_alloc(CACHE_LINE_SIZE, buf_size);
    bound = (int *)aligned_alloc(CACHE_LINE_SIZE, bound_size);
#endif

    // For aligned data, we use another pointer to do regular indexing
    double *ptr_vcurr, *ptr_vnext, *ptr_buf_scale;
    double *ptr_buf_diag;
    int offset_v, offset_buf;
    offset_v = (vlen-1)%(CACHE_LINE_SIZE/sizeof(double))==0 ? 0
                    :(CACHE_LINE_SIZE/sizeof(double))-(vlen-1)%(CACHE_LINE_SIZE/sizeof(double));
    offset_buf = 2*(vlen-1)%(CACHE_LINE_SIZE/sizeof(double))==0 ? 0
                    :(CACHE_LINE_SIZE/sizeof(double))-2*(vlen-1)%(CACHE_LINE_SIZE/sizeof(double));
    ptr_vcurr = v_curr + offset_v;
    ptr_vnext = v_next + offset_v;
    ptr_buf_scale = buf_scale + offset_v;
    ptr_buf_diag = buf_diag + offset_buf;

    // Make boundary information, offset by 1
    lb_margin = floor(lb_margin/dmesh+0.5);
    ub_margin = floor(ub_margin/dmesh+0.5);
    bound[0] = (int)(lb_margin+floor(lb_change[0]/dmesh+0.5));
    bound[1] = (int)(ub_margin-floor(ub_change[0]/dmesh+0.5));
    for (int t=1; t<=stepnum; t++) {
        bound[2*t] = (int)(lb_margin+floor(lb_change[t-1]/dmesh+0.5));
        bound[2*t+1] = (int)(ub_margin-floor(ub_change[t-1]/dmesh+0.5));
    }

    // Copy v_init into correct position
    memcpy(ptr_vcurr, v_init, vlen*sizeof(double));

    // Calculate bound crossing information at time point zero
    get_prob_info(ptr_vcurr, bound, 0, vlen, stepnum, num_output, prob_notcross, prob_cross, prob_gtzero, full_prob);

    // Main computation, param saved for computation
    double last_param[2] = {0.0, 0.0};
    double B, C;
    double *temp_v;
    for (int t=1; t<=stepnum; t++) {
        B = static_mu ? -mu[0] : -mu[t-1];
        C = static_noise ? 0.5*sigma[0]*sigma[0] : 0.5*sigma[t-1]*sigma[t-1];

#ifdef __GNUC__
        if (UNLIKELY( !(B==last_param[0] && C==last_param[1] && t>1) )) {
            // Make a new tridiag coeff and factor it in advanc (less likely to happen)
            make_tridiag_factor(B, C, alpha, dt, dmesh, vlen, v_lut, dl, d, du, ptr_buf_diag, ptr_buf_scale);
            last_param[0] = B;
            last_param[1] = C;
        }
#else
        if (!(B==last_param[0] && C==last_param[1] && t>1)) {
            // Make a new tridiag coeff and factor it in advance
            make_tridiag_factor(B, C, alpha, dt, dmesh, vlen, v_lut, dl, d, du, ptr_buf_diag, ptr_buf_scale);
            last_param[0] = B;
            last_param[1] = C;
        }
#endif

        // Solve linear systems
        solve_linear_eq(ptr_vnext, ptr_vcurr, ptr_buf_diag, ptr_buf_scale, vlen);

        // Get probability and clear crossed-bound value
        get_prob_info(ptr_vnext, bound, t, vlen, stepnum, num_output, prob_notcross, prob_cross, prob_gtzero, full_prob);

        // Swap v vectors
        temp_v = ptr_vnext;
        ptr_vcurr = temp_v;
        ptr_vnext = ptr_vcurr; 
    }

    // Always save v_final
    memcpy(v_final, ptr_vcurr, vlen*sizeof(double));

    // Free memory
#ifdef _MSC_VER         
    _aligned_free(v_curr);
    _aligned_free(v_next);
    _aligned_free(buf_scale);
    _aligned_free(dl);
    _aligned_free(d);
    _aligned_free(du);
    _aligned_free(buf_diag);
    _aligned_free(bound);
#else                   
    free(v_curr);
    free(v_next);
    free(buf_scale);
    free(dl);
    free(d);
    free(du);
    free(buf_diag);
    free(bound);
#endif

}


static void factor_coeff_matrix (double *dl, double *d, double *du, 
                            double *buf_diag, double *buf_scale, int vlen)
{
    register double scale;

    for (int i=1; i<vlen; i++) {
        scale = dl[i-1]/d[i-1];
        d[i] = d[i]-scale*du[i-1];
        buf_scale[i] = scale;
    }

#ifndef _MSC_VER
    #pragma omp simd aligned(d:64)
#endif
    for (int i=0; i<vlen; ++i) {
        d[i] = 1.0/d[i];
    }

#ifndef _MSC_VER
    #pragma omp simd aligned(d:64)
#endif
    for (int i=0; i<vlen; ++i) {
        // [d,du]
        buf_diag[i*2] = d[i];
        buf_diag[i*2+1] = du[i];
    }
}


static inline void solve_linear_eq (double *x, double *b, double *buf_diag, double *buf_scale, int vlen)
{   
    register double last_elem;
    register double x1, x2, x3, x4;
    register int i;

    // forward elimination
    last_elem = b[0];
    for (i=1; i<vlen; i++) {
        last_elem = b[i]-buf_scale[i]*last_elem;
        b[i] = last_elem;
    }

    // backward subsitution
    x1 = b[vlen-1]*buf_diag[2*(vlen-1)];
    x[vlen-1] = x1;
    i = vlen-2;
    _mm_prefetch(b+i-7, _MM_HINT_T0);
    _mm_prefetch(buf_diag+2*(i-7), _MM_HINT_T0);
    
    while (i-3 >= 0) {
        // prepare cache line for write in this iter
        #ifdef _MSC_VER
            _mm_prefetch(x+i-3, _MM_HINT_T0);
        #else
            _mm_prefetch(x+i-3, _MM_HINT_ET0);
        #endif
        
        // manual loop unrolling
        x4 = -(buf_diag[2*i+1]*buf_diag[2*i])*x1+(b[i]*buf_diag[2*i]);
        x[i] = x4;
        x3 = -(buf_diag[2*(i-1)+1]*buf_diag[2*(i-1)])*x4+(b[i-1]*buf_diag[2*(i-1)]);
        x[i-1] = x3;
        x2 = -(buf_diag[2*(i-2)+1]*buf_diag[2*(i-2)])*x3+(b[i-2]*buf_diag[2*(i-2)]);
        x[i-2] = x2;
        x1 = -(buf_diag[2*(i-3)+1]*buf_diag[2*(i-3)])*x2+(b[i-3]*buf_diag[2*(i-3)]);
        x[i-3] = x1;

        // prepare cache line for the next iter
        _mm_prefetch(b+i-7, _MM_HINT_T0);
        _mm_prefetch(buf_diag+2*(i-7), _MM_HINT_T0);     
        
        i -= 4;
    }

    for (; i>=0; i--) {
        x1 = (b[i]*buf_diag[2*i])-(buf_diag[2*i+1]*buf_diag[2*i])*x1;
        x[i] = x1;
    }
}


static void get_prob_info (double *ptr_v, int *bound, int t, int vlen, int stepnum, int num_output,
        double *prob_notcross, double *prob_cross, double *prob_gtzero, double *full_prob)
{
    register double temp_prob = 0.0;
    register int offset = vlen-bound[2*t+1];    // for upper bound
    register bool prob_cleared = false;

    switch (num_output) {
        case 5:
            // Do nothing for now
            
        case 4:
            // Probability greater than zero
            temp_prob = 0.0;
            #ifndef _MSC_VER
                #pragma omp simd reduction(+:temp_prob)
            #endif
            for (int i=vlen/2+1; i<vlen; i++) {
                temp_prob += fabs(ptr_v[i]);
            }
            prob_gtzero[t] = temp_prob + 0.5*ptr_v[vlen/2];

        case 3:
            // Crossing the lower (negative) bound
            temp_prob = 0.0;
            #ifndef _MSC_VER
                #pragma omp simd reduction(+:temp_prob)
            #endif
            for (int i=0; i<bound[2*t]; i++) {
                temp_prob += fabs(ptr_v[i]);
            }
            prob_cross[t] = temp_prob;

            // Crossing upper (positive) bound
            temp_prob = 0.0;
            #ifndef _MSC_VER
                #pragma omp simd reduction(+:temp_prob)
            #endif
            for (int i=0; i<bound[2*t+1]; i++) {
                temp_prob += fabs(ptr_v[i+offset]);
            }
            prob_cross[stepnum+1+t] = temp_prob;

            // Clear bound-crossed value
            #ifndef _MSC_VER
                #pragma omp simd
            #endif
            for (int i=0; i<bound[2*t]; i++) {
                ptr_v[i] = 0.0;
            }

            #ifndef _MSC_VER
                #pragma omp simd
            #endif
            for (int i=0; i<bound[2*t+1]; i++) {
                ptr_v[i+offset] = 0.0;
            }
            prob_cleared = true;

        case 2:
            // Not crossing bound
            temp_prob = 0.0;
            #ifndef _MSC_VER
                #pragma omp simd reduction(+:temp_prob)
            #endif
            for (int i=bound[2*t]; i<vlen-bound[2*t+1]; i++) {
                temp_prob += fabs(ptr_v[i]);
            }
            prob_notcross[t] = temp_prob;

        case 1:
            // Clear bound-crossed value ONLY
            if ( !prob_cleared ) {
                #ifndef _MSC_VER
                    #pragma omp simd
                #endif
                for (int i=0; i<bound[2*t]; i++) {
                    ptr_v[i] = 0.0;
                }

                #ifndef _MSC_VER
                    #pragma omp simd
                #endif
                for (int i=0; i<bound[2*t+1]; i++) {
                    ptr_v[i+offset] = 0.0;
                }
                prob_cleared = true;
            }
            break;
    }

    // Copy the whole v_vector
    if (num_output == 5 && t > 0) {
        #ifndef _MSC_VER
            #pragma omp simd
        #endif
        for (int i=0; i<vlen; i++) {
            full_prob[(t-1)*vlen+i] = fabs(ptr_v[i]);
        }
    }

    prob_cleared = false;

}

static void make_tridiag_factor (double B, double C, double alpha, double dt, double dmesh, int vlen, double *v_lut,
                                double *dl, double *d, double *du, double *buf_diag, double *buf_scale)
{   
    double w, W;
    double W_plus, W_minus;
    double temp_dl, temp_d, temp_du;

    const double coeff_1 = dt*C/(1.0*dmesh*dmesh);
    const double coeff_2 = dt/(1.0*dmesh);
    const double coeff_3 = C/dmesh;

    if (alpha==0.0) {
        // no leakage
        if (C != 0.0) {
            w = 0.5*B*dmesh/C;
            W = w==0.0 ? 1.0 : w/sinh(w);
            W_plus = W * exp(w);
            W_minus = W * exp(-w);
            temp_dl = coeff_1 * W_minus;
            temp_d = 1.0 + coeff_2*(coeff_3*W_plus + coeff_3*W_minus);
            temp_du = coeff_1 * W_plus;
        } else {
            temp_dl = 0.0;
            temp_d = 1.0;
            temp_du = 0.0;
        }

        #ifndef _MSC_VER
            #pragma omp simd aligned(dl,d,du:64)
        #endif
        for (int i=0; i<vlen; i++) {
            dl[i] = temp_dl;
            d[i] = temp_d;
            du[i] = temp_du;
        }
    } else {
        // leaky condition
        #ifndef _MSC_VER
            #pragma omp simd aligned(dl,d,du:64)
        #endif
        for (int i=0; i<vlen; i++) {
            double Bm = alpha*(v_lut[i]-dmesh*0.5)+B;
            w = 0.5*Bm*dmesh/C;
            W = w==0.0 ? 1.0 : w/sinh(w);
            W_plus = W * exp(w);
            W_minus = W * exp(-w);
            dl[i] = coeff_1 * W_minus;
            d[i] = 1.0 + coeff_2*coeff_3*(W_plus+W_minus);
            du[i] = coeff_1 * W_plus;
        }
    }

    // factor coeffs in advance
    factor_coeff_matrix (dl, d, du, buf_diag, buf_scale, vlen);
}