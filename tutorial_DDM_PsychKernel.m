clear;

addpath functions/;


speeded = true;
    % if you want to run fit a bit faster in this tutorial, turn this on. 

load('data/example_PK_data.mat');
% example data of motion direction discrimination task
%
% trial_data.Coh: Coherence level (signed)
% trial_data.Stim: Stimulus category (1 or 2)
% trial_data.Fluc: stim fluctuation (the same unit as coherence)
% trial_data.Resp: Subject's choice (1 or 2)
% trial_data.RT: Reaction time (msec)


%% ######### fit data #########

%% define random initial point
seed = 20; % change this for each initial point
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));

k = rand(1,3) * .2; % sensitivity
B = 20 + rand * 20; % bound height
T0 = 450 + rand * 150; % mean non-decision time (ms)
T0_sigma = T0  * .3;   % SD of non-decision time (ms)
urg_half = 1;   % unused.
urg_asymp = 0;  % unused.
guess = [k, B, T0, T0_sigma, urg_half, urg_asymp];
fixed = [zeros(1, length(k)), 0 0 0 1 1]; % fixed parameter (here, fix T0_sigma and urgency parameters)

%% set options

clear opt;

opt.fitwhat = 'RESPRT';
    % fit to: PC (choice), RT (reaction time), PCRT (both choice and RT)
opt.nfeature = 3;

opt.fluc_t = 1000/75 * 8; % duration of each noise cycle (msec)
    % in this example, 75Hz monitor with 8 frame of each stimulus.
    
opt.temporary_file = 'data/example_PK_fit_temporary.mat';

if speeded
    opt.TolX = 5e-1; % termination criterion for parameters change
    opt.TolFun = 1; % termination criterion for logLL change
else
    opt.TolX = 1e-2; % termination criterion for parameters change
    opt.TolFun = 1e-1; % termination criterion for logLL change
end

%% run fit

fprintf('Running DDM fit.\n');
[model_param, modelLL, trial_data] = fit_DDM_Kernel(trial_data, guess, fixed, opt);
save('data/example_PK_fit.mat', 'model_param', 'modelLL', 'trial_data');


%% ######### show result #########

show_Psych_Chrono('data/example_PK_fit.mat');


%% ######### generate model psych kernel #########


%% set options

clear opt;

opt.nfeature = 3;
opt.iters = 10000;
opt.t_max = 5000;
opt.fluc_t = 1000/75 * 8; % duration of each noise cycle (msec)
opt.coh = 0;
opt.sigma = 0.2; % stimulus fluctuation is 20% SD (0.2)
opt.dec_noise = 1;


%% set model param

S = load('data/example_PK_fit.mat');
model_param = S.model_param.final;

opt.k = model_param(1:opt.nfeature);
idx = opt.nfeature;
opt.B = abs(model_param(idx + 1))' * [-1 1];
opt.non_dec_time = abs(round(model_param(idx + 2)));
opt.non_dec_time_sd = abs(round(model_param(idx + 3)));

opt.bt_mediant = ceil(median(S.trial_data.RT(S.trial_data.Coh==0))/opt.fluc_t);

fprintf('Median stimulus frame len: %d\n', opt.bt_mediant);

%% run

Nsim = 10; % as one simulation makes a large matrix, split into multiple simulations

sim_all = cell(Nsim,1);
for n=1:Nsim
    fprintf('Running simulation %d/%d...\n', n, Nsim);
    if n == 1
        sim = DDM_Kernel_sim(opt);
    else
        sim1 = DDM_Kernel_sim(opt);
        for iF=1:opt.nfeature
            sim.kernel{iF}.motion{1} = sim.kernel{iF}.motion{1} + sim1.kernel{iF}.motion{1};
            sim.kernel{iF}.motion{2} = sim.kernel{iF}.motion{2} + sim1.kernel{iF}.motion{2};
            sim.kernel{iF}.choice{1} = sim.kernel{iF}.choice{1} + sim1.kernel{iF}.choice{1};
            sim.kernel{iF}.choice{2} = sim.kernel{iF}.choice{2} + sim1.kernel{iF}.choice{2};
        end
    end
end
for iF=1:opt.nfeature
    sim.kernel{iF}.motion{1} = sim.kernel{iF}.motion{1}/Nsim;
    sim.kernel{iF}.motion{2} = sim.kernel{iF}.motion{2}/Nsim;
    sim.kernel{iF}.choice{1} = sim.kernel{iF}.choice{1}/Nsim;
    sim.kernel{iF}.choice{2} = sim.kernel{iF}.choice{2}/Nsim;
end

save('data/example_PK_model_kernel.mat', 'sim');


%% ######### show result #########


clear opt;

opt.max_coh = .15;
opt.smoothing = 3;
opt.fluc_t = 8/75;
opt.positive_cat = 2;
opt.nfeature = 3;
fh = show_Kernel('data/example_PK_model_kernel.mat', trial_data, opt);




