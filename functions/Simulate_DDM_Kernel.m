function sim = Simulate_DDM_Kernel(bestfitfile, optdef)

%% simulation setting
def.Ntrials = 1e5; % total number of trials
def.nfeature = 3;
opt = safeStructAssign(def, optdef);

%% model parameters

% default parameter for Okazawa (2021) J Neurosci
p.iters = 10000;
p.t_max = 5000;
p.fluc_t = 8000/75;
p.w = 1;
p.coh = 0;
p.nfeature = 3;
p.sigma = 0.2; % stimulus fluctuation is 20% SD (0.2)
p.dec_noise = 1;
p.error_no_reach = false;
p = safeStructAssign(p, optdef);

Nsim = opt.Ntrials/p.iters;
if Nsim < 1
    p.iters = opt.Ntrials;
    Nsim = 1;
else
    Nsim = ceil(Nsim);
    p.iters = opt.Ntrials/Nsim;
end

%% load fit parameters
S = load(bestfitfile);
model_param = S.fitresult.model_param.final;

if size(S.fitresult.trial_data.fluc{1},2) ~= opt.nfeature
    error('The number of features in fluc does not match nfeature.');
end

p.k = model_param(1:opt.nfeature);
idx = opt.nfeature;
p.B = abs(model_param(idx + 1))' * [-1 1];
p.non_dec_time = abs(round(model_param(idx + 2)));
p.non_dec_time_sd = abs(round(model_param(idx + 3)));

p.bt_mediant = ceil(median(S.fitresult.trial_data.rt(S.fitresult.trial_data.coh==0) * 1e3)/p.fluc_t);
fprintf('Median stimulus frame len: %d\n', p.bt_mediant);

%% run simulation
sim_all = cell(Nsim,1);
for n=1:Nsim
    fprintf('Running simulation %d/%d trials...\n', p.iters * n, p.iters * Nsim);
    sim_all{n} = DDM_Kernel_sim(p);
end

% average kernel across N repetitions
sim = sim_all{1};
for iF=1:opt.nfeature % features
    for n=2:Nsim
        sim.kernel{iF}.motion{1} = sim.kernel{iF}.motion{1} + sim_all{n}.kernel{iF}.motion{1};
        sim.kernel{iF}.motion{2} = sim.kernel{iF}.motion{2} + sim_all{n}.kernel{iF}.motion{2};
        sim.kernel{iF}.choice{1} = sim.kernel{iF}.choice{1} + sim_all{n}.kernel{iF}.choice{1};
        sim.kernel{iF}.choice{2} = sim.kernel{iF}.choice{2} + sim_all{n}.kernel{iF}.choice{2};
    end
    sim.kernel{iF}.motion{1} = sim.kernel{iF}.motion{1}/Nsim;
    sim.kernel{iF}.motion{2} = sim.kernel{iF}.motion{2}/Nsim;
    sim.kernel{iF}.choice{1} = sim.kernel{iF}.choice{1}/Nsim;
    sim.kernel{iF}.choice{2} = sim.kernel{iF}.choice{2}/Nsim;
end



