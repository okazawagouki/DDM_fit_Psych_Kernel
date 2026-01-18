function [model_param, modelLL, trial_data] = fit_DDM(trial_data, guess, fixed, options)
    
    default.fitwhat = 'PCRT';
        % fit to: PC (choice), RT (reaction time), PCRT (both choice and RT)
    default.T0_sigma = 0.3;
        % if set, SD of non-decision time is set to Mean * T0_sigma.
        % if nan, SD is also fit to data.
    default.urgency = 'hyperbolic';
        % shape of urgency signal (collapsing bound). hyperbolic of sigmoid
        % use hyperbolic unless there is any specific reason.
    default.feedback = false; % show error trace figure
    default.max_t = []; % maximum RT (msec). If not set, use max RT.
    default.max_dx = 0.05; % sampling precision of decision variable in finite difference method
    default.dt = 0.1; % sampling precision of time (ms) in finite difference method
    default.TolX = 1e-3; % termination criterion for parameters change
    default.TolFun = 1e-2; % termination criterion for logLL change
    default.urg_half_range = [300 1e4]; % possible range of urgency half time (ms)
    default.scale_factor = [0.5, 50, 500, 100, 500, 50, 1];
    
    if ~exist('options', 'var') || isempty(options)
        opt = default;
    else
        opt = safeStructAssign(default, options);
    end
    
    model_param = struct('init', guess, 'fixed', fixed, 'final', [], 'se', []);
    
    if isempty(opt.max_t)
        opt.max_t = ceil(max(trial_data.RT));
    end
    
    trial_num = length(trial_data.Stim);
    Coh = abs(trial_data.Coh);
    coh_set = unique(Coh);
    LX = cell(length(coh_set), 2);
    for cc = 1 : length(coh_set)
        Lcoh = Coh==coh_set(cc);
        LX{cc,1} = Lcoh & trial_data.Resp==trial_data.Stim;
        LX{cc,2} = Lcoh & trial_data.Resp~=trial_data.Stim;
    end
    
    call_num = 1;
    if all(fixed==1)
        fh = -1;
        fprintf('\n\nall parameters are fixed, no optimization will be performed\n\n');
        modelLL = -fitmodel_MLEerr([]);
        return;
    end
    
    if opt.feedback
        fh = figure;
        hold on;
        xlabel('Call number' );
        ylabel('-LL');
        set(gca, 'xlim', [0 50]);
        err_hist = [];
        drawnow;
    end
    
    guess = guess ./ opt.scale_factor(1:length(guess)); % scaling initial guess
    options = optimset('Display', 'final', 'MaxFunEvals', 500*sum(fixed==0), 'MaxIter', 500*sum(fixed==0), 'TolX', opt.TolX, 'TolFun', opt.TolFun);
    [p, fval] = fminsearch(@fitmodel_MLEerr, guess(fixed==0), options);
    model_param.final = getParam(p, guess, fixed);
    modelLL = -fval;
    
    opt.feedback = false;
    
    guess = model_param.final ./ opt.scale_factor(1:length(guess));
    fixed = ones(size(guess));
    fitmodel_MLEerr([]);
    
    
        %error function
    function err = fitmodel_MLEerr(param)
        param = getParam(param,guess,fixed);
        
        k = param(1);
        B = abs(param(2));  %don't accept negative bound heights
        if B<1
            err = Inf;
            return;
        end
        sigma = 1;
        T0 = abs(param(3));
        T0_sigma = abs(param(4));
        if strcmp(opt.urgency,'hyperbolic')
            urg_half = abs(param(5));
            urg_asymp = abs(param(6));
        elseif strcmp(opt.urgency,'sigmoid')
            urg_half = abs(param(5));
            urg_asymp = abs(param(6));
            urg_slope = param(7);
        else
            error('unknown function for urgency!');
        end
        if urg_asymp < 0 || urg_asymp > B || ...
                (urg_asymp > 0 && (urg_half < opt.urg_half_range(1) || urg_half > opt.urg_half_range(2)))
                % Urgency parameters should not exceed reasonable ranges
            err = Inf;
            return;
        end
        
        conv_kernel_T0 = normpdf(0:T0+3*T0_sigma, T0, T0_sigma)';
        conv_kernel_T0 = conv_kernel_T0/sum(conv_kernel_T0);
        
        if strcmp(opt.urgency,'hyperbolic')
            urgency = (0:opt.max_t)'*urg_asymp./((0:opt.max_t)'+urg_half);
        elseif strcmp(opt.urgency,'sigmoid')
            tr = (0:opt.max_t)'-urg_half;
            urgency = urg_asymp*((1-exp(-urg_slope.*tr))./(1+exp(-urg_slope.*tr))+1)/2;
        end
        bound_change = [min(urgency,B-6*opt.max_dx) max(-urgency,-B+6*opt.max_dx)];
            %make the spatial grid for the probability density
        bound = repmat([-B B],[size(bound_change,1),1]) + bound_change;
        bound_range = [min(bound(:)) max(bound(:))];
        dx = min(opt.max_dx, diff(bound_range)/100);
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
        
        Ptb_coh = zeros(opt.max_t, 2, length(coh_set));                % time * bound * coh
        for c = 1 : length(coh_set)
            mu = k*coh_set(c);
            uinit = delta;
            b_change = interp1(0:opt.max_t, bound_change, opt.dt:opt.dt:opt.max_t);
            [~,~,Ptb] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin, opt.dt);
                %store the arrays for each coherence level, use 1 ms time resolution 
            Ptb_coh(:,:,c) = [local_sum(Ptb(2:end,1),round(1/opt.dt)), local_sum(Ptb(2:end,2),round(1/opt.dt))];
                %add the non-decision time
            Ptb_coh(:,1,c) = convT0(Ptb_coh(:,1,c), conv_kernel_T0);
            Ptb_coh(:,2,c) = convT0(Ptb_coh(:,2,c), conv_kernel_T0);
        end
        
            %calculate the expected probability correct and the expected decision time 
        t = 1 : opt.max_t;
        
        trial_data.fit_pc = nan(trial_num, 1);
        trial_data.fit_rt = nan(trial_num, 1);
        for c = 1 : length(coh_set)
            L = LX{c,1}|LX{c,2};
            trial_data.fit_pc(L) = sum(Ptb_coh(:,2,c)) / (sum(Ptb_coh(:,1,c))+sum(Ptb_coh(:,2,c)));
            trial_data.fit_rt(L) = t*(Ptb_coh(:,2,c)+Ptb_coh(:,1,c))/sum(Ptb_coh(:,2,c)+Ptb_coh(:,1,c));
        end
        trial_data.fit_pchoice1 = nan(trial_num, 1);
        trial_data.fit_pchoice1(trial_data.Stim==2) = 1 - trial_data.fit_pc(trial_data.Stim==2);
        trial_data.fit_pchoice1(trial_data.Stim==1) = trial_data.fit_pc(trial_data.Stim==1);
        
        
            %calculate the error 
        P = Ptb_coh(:);
        if all(P(:)==0)
            err = Inf;
            return;
        end
        minP = min(P(P~=0));
        LL = 0;
        n = 0;
        for c = 1 : length(coh_set)
            switch opt.fitwhat
                case 'PCRT'
                    LL = LL + ...
                         sum(log(max(Ptb_coh(round(trial_data.RT(LX{c,1})),2,c),minP))) + ...
                         sum(log(max(Ptb_coh(round(trial_data.RT(LX{c,2})),1,c),minP)));
                    n = n + sum(LX{c,1}) + sum(LX{c,2});
                case 'PC'
                    if coh_set(c)==0
                        continue;
                    end
                    LL = LL + ...
                         sum(log(max(trial_data.fit_pc(LX{c,1}),minP))) + ...
                         sum(log(max(1-trial_data.fit_pc(LX{c,2}),minP)));
                    n = n + sum(LX{c,1}) + sum(LX{c,2});
                case 'RT'
                    Prt_coh = Ptb_coh(:,2,c) + Ptb_coh(:,1,c);
                    LL = LL + ...
                         sum(log(max(Prt_coh(round(trial_data.RT(LX{c,1}|LX{c,2}))),minP)));
                    n = n + sum(LX{c,1}) + sum(LX{c,2});
            end
        end
        
        err = -LL;
        
        %print progress report!
        fprintf('\n\n\n****************************************\n');
        fprintf('run %d\n', call_num);
        if strcmp(opt.urgency,'hyperbolic')
            fprintf('\tk= %g\n\tB= %g\n\tT0= %g\n\tT0_sigma= %g\n\turg_half= %g\n\turg_asymp= %g\n', ...
                    k, B, T0, T0_sigma, urg_half, urg_asymp);
        elseif strcmp(opt.urgency,'sigmoid')
            fprintf('\tk= %g\n\tB= %g\n\tT0= %g\n\tT0_sigma= %g\n\turg_half= %g\n\turg_asymp= %g\n\turg_slope= %g\n', ...
                    k, B, T0, T0_sigma, urg_half, urg_asymp, urg_slope);
        end
        fprintf('err: %f\n', err);
        fprintf('number of processed trials: %d\n', n);
        
        if opt.feedback
            err_hist = [err_hist; err];
            if ishandle(fh)
                figure(fh);
            else
                fh = figure;
            end
            plot(err_hist, 'k', 'marker', '.', 'markers', 10);
            set(gca, 'xlim', [0 max(50, length(err_hist)+10)]);
            drawnow;
        end
        call_num = call_num + 1;
    end

                %this function retrieves the full parameter set given the adjustable and
                %fixed parameters
    function param2 = getParam ( param1 , guess , fixed )
        param2(fixed==0) = param1;              %get adjustable parameters from param1
        param2(fixed==1) = guess(fixed==1);     %get fixed parameters from guess
        param2 = param2 .* opt.scale_factor(1:length(param2));    %scale back parameters
        if ~isnan(opt.T0_sigma)
            param2(4) = opt.T0_sigma * param2(3);
        end
    end
end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function P = convT0(P, kernel)
    len = length(P);
    P = conv(P, kernel);
    P = P(1:len);
end




function newd = local_sum ( oldd , sumlen )
    newd_len = floor(length(oldd)/sumlen);
    oldd_len = newd_len * sumlen;
    newd = sum(reshape(oldd(1:oldd_len),[sumlen newd_len]))';
end


function targ_struct = safeStructAssign(targ_struct, ref_struct, exclude_fields, cutdown)
    if nargin<3 || isempty(exclude_fields)
        exclude_fields = {};
    elseif ~iscell(exclude_fields)
        exclude_fields = {exclude_fields};
    end

    if nargin<4
        cutdown = true;
    end

        %ensure targ_struct and ref_struct are structures and that targ_struct has the same
        %size as ref_struct
    if ~isstruct(targ_struct) || ~isstruct(ref_struct)
        return;
    else
        if length(targ_struct)<length(ref_struct)
            for i = length(targ_struct)+1 : length(ref_struct)
                targ_struct(i) = targ_struct(end);
            end
        elseif length(targ_struct)>length(ref_struct) && cutdown
            targ_struct = targ_struct(1:length(ref_struct));        
        end
    end

        %override existing values of targ_struct with those in ref_struct
    names = fieldnames(ref_struct);
    for i = 1 : length(names)
        if isfield(targ_struct,names{i}) && ~ismember(names{i},exclude_fields)
            for j = 1 : length(ref_struct)
                targ_struct(j).(names{i}) = ref_struct(j).(names{i});
            end
        end
    end
end





