clear;

addpath functions/;

speeded = true;
    % if you want to run fit a bit faster in this tutorial, turn this on. 

load('data/example_data.mat');
% example data of motion direction discrimination task
%
% trial_data.Coh: Coherence level (signed)
% trial_data.Stim: Stimulus category (1 or 2)
% trial_data.Resp: Subject's choice (1 or 2)
% trial_data.RT: Reaction time (msec)


%% ######### fit data #########

%% define random initial point
seed = 1; % change this for each initial point
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));

k = .2 + rand * .2; % sensitivity
B = 20 + rand * 20; % bound height
T0 = 250 + rand * 150; % mean non-decision time (ms)
T0_sigma = T0 * .3;   % unused. SD of non-decision time (see below)
urg_half = 1;   % unused.
urg_asymp = 0;  % unused.
guess = [k, B, T0, T0_sigma, urg_half, urg_asymp];
fixed = [0 0 0 1 1 1]; % fixed parameter (here, fix T0_sigma and urgency parameters)

%% set options
opt.fitwhat = 'PCRT';
    % fit to: PC (choice), RT (reaction time), PCRT (both choice and RT)
opt.T0_sigma = 0.3;
    % if set, SD of non-decision time is set to Mean * T0_sigma.
    % if nan, SD is also fit to data.
opt.feedback = true; % show error trace figure

if speeded
    opt.TolX = 5e-1; % termination criterion for parameters change
    opt.TolFun = 1; % termination criterion for logLL change
else
    opt.TolX = 1e-2; % termination criterion for parameters change
    opt.TolFun = 1e-1; % termination criterion for logLL change
end

%% run fit

fprintf('Running DDM fit.\n');
[model_param, modelLL, trial_data] = fit_DDM(trial_data, guess, fixed, opt);
save('data/example_fit.mat', 'model_param', 'modelLL', 'trial_data');


%% ######### show result #########

show_Psych_Chrono('data/example_fit.mat');









