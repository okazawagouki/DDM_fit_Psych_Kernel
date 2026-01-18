function [model_param, modelLL, trial_data, exitflg] = fit_DDM_Kernel(trial_data, guess, fixed, options)
    
    default.fitwhat = 'RESPRT'; % currently support RESPRT only
    default.parfor = true;
    default.T0_sigma = nan;
    default.urgency = 'hyperbolic';
    default.max_t = 5000;
    default.max_dx = 0.05;
    default.dt = 0.1;
    default.dt_seq = [];
    default.sigma = 1;
    default.TolX = 1e-3;
    default.TolFun = 1e-2;
    default.fluc_t = 100; % duration of each noise frame(msec)
    default.pos_cho = 2; % correct choice for positive coherence
    default.nfeature = 3;
    default.scale_factor = [0.5, 50, 500, 100, 500, 50, 1];
    default.explore = 0;
    
    default.maxiter = 500*sum(fixed==0);
    default.temporary_file = []; % if empty, do not use temp file
    default.urg_half_range = [0 Inf]; % possible range of urgency half time (ms)
    
    
    if ~exist('options', 'var')
        opt = default;
    else
        opt = safeStructAssign(default, options);
    end
    
    if ~isnan(opt.T0_sigma) && fixed(opt.nfeature+3) == 0
        warning('T0_sigma option is not NaN. Assuming that T0_sigma is fixed');
        fixed(opt.nfeature+3) = 1;
    end
        
    if length(opt.scale_factor) == 7 && opt.nfeature > 1
        % scale factor only defined for one feature
        opt.scale_factor = [opt.scale_factor(1) * ones(1, opt.nfeature), opt.scale_factor(2:end)];
    end
    
    model_param = struct('init', guess, 'fixed', fixed, 'final', [], 'se', []);
    guess = guess ./ opt.scale_factor(1:length(guess)); % scaling initial guess
    
    call_num = 1;
    if all(fixed==1)
        fprintf('\n\nAll parameters are fixed, no optimization will be performed.\n\n');
        gen_fitdata = true;
        modelLL = -fitmodel_MLEerr([]);
        return;
    end
    
    gen_fitdata = false; % during optimization, no need to generate fit data (making it false speed up)
    options = optimset('Display', 'none', 'MaxFunEvals', 300*sum(fixed==0), 'MaxIter', opt.maxiter, 'TolX', opt.TolX, 'TolFun', opt.TolFun);
    [p, fval, exitflg] = fminsearch(@fitmodel_MLEerr, guess(fixed==0), options);
    model_param.final = getParam(p, guess, fixed);
    modelLL = -fval;
    
    gen_fitdata = true; % now generate fit data using fit parameters
    fitmodel_MLEerr(p);
    
    
        %error function
    function err = fitmodel_MLEerr(param)
        param = getParam(param,guess,fixed);
        err = fitHist(opt.temporary_file, param, [], 1);
        if ~isempty(err) && ~gen_fitdata
            call_num = call_num + 1;
            return;
        end
        tic;
        
        k = param(1:opt.nfeature);
        idx = opt.nfeature;
        B = param(idx + 1);
        sigma = opt.sigma;
        T0 = param(idx + 2);
        T0_sigma = param(idx + 3);
        if strcmp(opt.urgency,'hyperbolic')
            urg_half = param(idx + 4);
            urg_asymp = param(idx + 5);
        elseif strcmp(opt.urgency,'sigmoid')
            urg_half = param(idx + 4);
            urg_asymp = param(idx + 5);
            urg_slope = param(idx + 6);
        else
            error('unknown function for urgency!');
        end
        if B <=0 || T0 <= 0 || T0_sigma <= 0 || urg_asymp < 0 || urg_asymp > B || ...
                (urg_asymp > 0 && (urg_half < opt.urg_half_range(1) || urg_half > opt.urg_half_range(2)))
                % parameters should not exceed reasonable ranges
            err = Inf;
            fitHist(opt.temporary_file, param, err, 1, toc);
            return;
        end
        
        conv_kernel_T0 = normpdf(0:T0+3*T0_sigma, T0, T0_sigma)';
        conv_kernel_T0 = conv_kernel_T0/sum(conv_kernel_T0);
        
        max_dx = opt.max_dx;
        if ~isempty(opt.dt_seq)
            ind = find(opt.dt_seq(:,2) > call_num,1);
            dt = opt.dt_seq(ind,1);
        else
            dt = opt.dt;
        end
        
        if strcmp(opt.urgency,'hyperbolic')
            urgency = (0:opt.max_t)'*urg_asymp./((0:opt.max_t)'+urg_half);
        elseif strcmp(opt.urgency,'sigmoid')
            tr = (0:opt.max_t)'-urg_half;
            urgency = urg_asymp*((1-exp(-urg_slope.*tr))./(1+exp(-urg_slope.*tr))+1)/2;
        end
        bound_change = [min(urgency,B-6*max_dx) max(-urgency,-B+6*max_dx)];
            %make the spatial grid for the probability density
        bound = repmat([-B B],[size(bound_change,1),1]) + bound_change;
        bound_range = [min(bound(:)) max(bound(:))];
        dx = min(max_dx, diff(bound_range)/100);
        b_margin = [-B-bound_range(1)+4*sigma; bound_range(2)-B+4*sigma];
        xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';
            % adjust dx so that the mesh has zero and is therefore symmetric 
        while ~any(xmesh==0)
            [~, I] = min(abs(xmesh));
            delta_mesh = xmesh(I);
            b_margin = b_margin + delta_mesh;
            xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)'; % recompute the mesh
        end
        if mod(length(xmesh),2)~=1
            error('the length of xmesh must be odd');
        end
        
            %the initial distribution is a delta function
        delta = zeros(size(xmesh));
        [~,I] = min(abs(xmesh));
        delta(I(1)) = 1;
        
        
        trial_num = length(trial_data.Coh);
        coh = trial_data.Coh;
        fluc = trial_data.Fluc;
        targ_cor = trial_data.Stim;
        response = trial_data.Resp;
        rt = trial_data.RT;
        
            %calculate the expected probability correct and the expected decision time 
        fit_pc = nan(trial_num, 1);
        fit_pchoice1 = nan(trial_num, 1); % p(choice == 1)
        fit_rt = nan(trial_num, 1);
        LL = nan(trial_num,1);
        
        b_change = interp1(0:opt.max_t, bound_change, dt:dt:opt.max_t);
        
        muid = ceil((0:dt:opt.max_t)/opt.fluc_t);
        muid = muid(2:end);
        
        if opt.parfor
            pfl = Inf;
        else
            pfl = 0;
        end
        
        parfor (c = 1:trial_num, pfl)
%         for c =1:trial_num
            mu = k(1) * fluc{c}(:,1);
            for f=2:opt.nfeature
                mu = mu + k(f) * fluc{c}(:,f);
            end
            mu = mu(muid(muid<=length(mu))); %#ok<*PFBNS>
            mu = mu(:);
            
            if gen_fitdata && ~opt.explore
                % for generating predicted choice and RT, the diffusion up
                % until max reaction time should be obtained.
                mu_ext = ones(size(b_change,1)-length(mu),1) * coh(c) * sum(k);
                    % generate coherence beyond fluctuation. Since the
                    % average of fluctuation is coh, we put average here to
                    % generate fit data.
                [~,~,Ptb] = FP4basic(xmesh, delta, [mu; mu_ext], sigma, b_change, b_margin, dt, 0);
            else
                % if you just want to get the likelihood of particular RT,
                % diffusion process beyond this RT is unnecessary. This
                % significantly speeds up the computaiton.
                [~,~,Ptb] = FP4basic(xmesh, delta, mu, sigma, b_change(1:length(mu),:), b_margin, dt);
            end
            
                %store the arrays for each coherence level, use 1 ms time resolution 
            if dt==1
                Ptb_coh = Ptb(1:end-1,:);
            else
                Ptb_coh = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
            end
                %add the non-decision time
            Ptb_coh(:,1) = convT0(Ptb_coh(:,1), conv_kernel_T0);
            Ptb_coh(:,2) = convT0(Ptb_coh(:,2), conv_kernel_T0);
            
                %calculate the expected probability correct and the expected decision time 
            if gen_fitdata && ~opt.explore
                t1 = 1:opt.max_t;
                if strcmp(opt.fitwhat,'RESPRT')
                    fit_rt(c) = t1*(Ptb_coh(:,2)+Ptb_coh(:,1))/sum(Ptb_coh(:,2)+Ptb_coh(:,1));
                    if opt.pos_cho == 2
                        fit_pchoice1(c) = sum(Ptb_coh(:,1))/(sum(Ptb_coh(:,1))+sum(Ptb_coh(:,2)));
                    else
                        fit_pchoice1(c) = sum(Ptb_coh(:,2))/(sum(Ptb_coh(:,1))+sum(Ptb_coh(:,2)));
                    end
                    %%expected pc dictated by the true experiment design
                    if targ_cor(c) == opt.pos_cho
                        fit_pc(c) = sum(Ptb_coh(:,2))/(sum(Ptb_coh(:,1))+sum(Ptb_coh(:,2)));
                    else
                        fit_pc(c) = sum(Ptb_coh(:,1))/(sum(Ptb_coh(:,1))+sum(Ptb_coh(:,2)));
                    end
                end
            end            
            P = Ptb_coh(:);
            minP = min(P(P~=0));
            
            %calculate the error 
            switch opt.fitwhat
                case 'RESPRT'
                    if response(c) == opt.pos_cho
                        LL(c) = log(max(Ptb_coh(floor(rt(c)),2),minP));
                    else
                        LL(c) = log(max(Ptb_coh(floor(rt(c)),1),minP));
                    end
            end
        end
        
        % save fit result to trial_data
        if gen_fitdata && ~opt.explore
            trial_data.fit_pc = fit_pc;
            trial_data.fit_pchoice1 = fit_pchoice1;
            trial_data.fit_rt = fit_rt;
            trial_data.LL = LL;
        end
        
        if opt.explore
            trial_data.LL = LL;
        end
        
        err = -sum(LL);
        if ~gen_fitdata
            fitHist(opt.temporary_file, param, err, 1, toc);
        end
        
        %print progress report!
        fprintf('\n\n\n****************************************\n');
        fprintf('run %d\n', call_num);
        ks = sprintf('%g ', k);
        if strcmp(opt.urgency,'hyperbolic')
            fprintf('\tk= %s\n\tB= %g\n\tT0= %g\n\tT0_sigma= %g\n\turg_half= %g\n\turg_asymp= %g\n', ...
                    ks, B, T0, T0_sigma, urg_half, urg_asymp);
        elseif strcmp(opt.urgency,'sigmoid')
            fprintf('\tk= %s\n\tB= %g\n\tT0= %g\n\tT0_sigma= %g\n\turg_half= %g\n\turg_asymp= %g\n\turg_slope= %g\n', ...
                    ks, B, T0, T0_sigma, urg_half, urg_asymp, urg_slope);
        end
        
        fprintf('err: %f\n', err);
        fprintf('number of processed trials: %d\n', trial_num);
        fprintf('Time of one cycle: %1.3fs\n', toc);
        call_num = call_num + 1;
    end
                %this function retrieves the full parameter set given the adjustable and
                %fixed parameters
    function param2 = getParam (param1, guess, fixed)
        if isempty(param1)
            param2 = guess;                      % get parameters from guess
            param2 = param2 .* opt.scale_factor(1:length(guess)); % scale back parameters
        else
            param2(fixed==0) = param1;           % get adjustable parameters from param1
            param2(fixed==1) = guess(fixed==1);  % get fixed parameters from guess
            param2 = param2 .* opt.scale_factor(1:length(guess)); % scale back parameters
            if ~isnan(opt.T0_sigma)
                param2(opt.nfeature + 3) = opt.T0_sigma * param2(opt.nfeature + 2);
            end
        end
    end


end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function P = convT0(P, kernel)
    len = length(P);
    P = conv(P, kernel);
    P = P(1:len);
end
