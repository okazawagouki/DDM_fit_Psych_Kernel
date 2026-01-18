# Fitting drift diffusion model to behavioral data (shared version)

---

## Overview

The code fits the drift diffusion model (DDM) to choice and reaction time data using a maximum likelihood estimation. On each iteration of parameter search, it solves the Fokker-Planck equation with finite difference method to compute a joint probability distribution of choice and reaction time. The model parameters include sensitivity, bound height, mean and SD of non-decision time, and urgency. Urgency is modeled using a hyperbolic function with two parameters: asymptotic value and half time (Churchland 2008). The code is developed in [Kiani lab@NYU](https://www.cns.nyu.edu/kianilab/Home.html).

## Requirement
The code is tested using Matlab 2018a (w/ statistics toolbox).  
FP4.mex should be compiled for the environment you run the code.


## Package

[tutorial\_DDM.m](./tutorial_DDM.m) - main code to run DDM fit.

[tutorial\_DDM\_PsychKernel.m](./tutorial_DDM_PsychKernel.m) - main code to run DDM fit when there is fluctuation in stimulus strength. With the fluctuation desion, psychophysical kernel can be calculated.

example_data.mat - example behavioral data of one subject from Purcell (2016) and Okazawa (2018).

example_PK_data.mat - example behavioral data (with fluctuation in stimulus strength) of three subjects from Luo (2025).

## Some cautions


1. The fittings in the tutorial could take more than half an hour, although they use a speeded option (it runs fits with relaxed convergence criteria). If you turn this option off, the fits would take even longer.
2. The pacakge contains example data and example results of model fits. If you run the tutorial, you will see the results of the model fits.
3. This tutorial runs only one fit using one initial starting point, but ideally you want to use ten or more starting points and choose the one that converges with the highest log likelihood.
4. Urgency parameters are particularly sensitive to the starting point. This makes sense because a broad range of parameters produce similar results. We use a hyperbolic function (Ua * (t/(t + Uh)), Ua: asymptotic value, Uh: half time); using this function, a flat decision bound can be created by either setting Ua = 0 or making Ua large but also making Uh extremely large or small. To obtain a reasonable convergence, you should choose multiple initial starting points within reasonable parameter ranges.
5. The confidence intervals of the parameters should be estimated using bootstrap. For each bootstrap, you should run fits from multiple starting points to find the best parameters. Thus, if you perform 100 bootstrap estimates, you end up running 100 x 10 fits, which could be very time consuming. This would be only possible with high computing cluster.


## Distribution

At this point, the package is not publicly distributed.

## Reference

When referring to a DDM fit in a manuscript, we usually cite the following papers.

**Fokker-Planck equation**:

Karlin, S. & Taylor, H. E. A (1981) *Second Course in Stochastic Processes* (Elsevier, Amsterdam).

**DDM fit**:

Kiani, R., and Shadlen, M.N. (2009). Representation of confidence associated with a decision by neurons in the parietal cortex. *Science* 324, 759–764.

Okazawa, G., Sha, L., Purcell, B.A., and Kiani, R. (2018). Psychophysical reverse correlation reflects both sensory and decision-making processes. *Nat. Commun.* 9, 3479.

**Urgency**:

Churchland, A. K., Kiani, R. & Shadlen, M. N. (2008) Decision-making with multiple alternatives. *Nat. Neurosci.* 11, 693–702.

Purcell, B. A. & Kiani, R. Neural mechanisms of post-error adjustments of decision policy in parietal cortex (2016). *Neuron* 89, 658–671.


---
Okazawa lab @ ION

[http://english.cebsit.cas.cn/lab/okazawagoki/research/](http://english.cebsit.cas.cn/lab/okazawagoki/research/)


[okazawa@nyu.edu](mailto:okazawa@nyu.edu)

Package developed by Luo Tianlin and Gouki Okazawa




