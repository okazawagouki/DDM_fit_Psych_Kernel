function sim = DDM_Kernel_sim(p)
tic;
% p is a struct to define parameter. See below for the default parameters.

%% default parameters
def.iters = 10000;   % number of simulated trials per stimulus strength
def.t_max = 5000;    % number of simulated time steps in each trial
def.dt = 1;          % time unit, in ms
def.fluc_t = 1;      % duration of one noise frame

def.nfeature = 3;
def.coh = 0;        % coherence (if not scalar, coh is randomly chosen for eath trial)
def.k = [1 1 1];    % sensitivity parameter (drift rate = coh * k)
def.k0 = 0;         % base drift of the accumulators, equivalent of **bias**
def.B = [-30 30];   % lower and upper bounds (1 x 2) or (t_max x 2)
def.sigma = 1;      % standard deviation of fluctuation of stimuli
def.dec_noise = 0;  % noise added to decision variable (decision noise)
def.w = 1;          % weight (1 x 1) or (t_max x 1)
def.non_dec_time = 0;    % average non-decision time
def.non_dec_time_sd = 0; % SD of non-decision time
def.seed = sum(clock * 100);
def.error_no_reach = true; % end the program with error when less than 95% of trials reach the bound.
def.bt_mediant = nan;

%% setup
if nargin < 1
    p = def;
else
    p = safeStructAssign(def, p);
end

    %fixed bound height over time. change it to simulate collapsing bounds
if size(p.B,1) == 1
    p.bound_height = repmat(p.B,[p.t_max,1]);
else
    p.bound_height = p.B;
    p = rmfield(p, 'B');
end

    %fixed weights over time. change it to simulate dynamic weighting of sensory evidence 
p.weight = ones(1, p.t_max) * p.w;



sim = struct();

if ~isnan(p.seed)
    %fix the random seed to reduce variability in the model exploration phase. 
    %For final simulations, you should use random seed
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',p.seed));
end

%% run simulation
    %loop through stimulus set
if isscalar(p.coh)
    coh = p.coh;
else
    idx = ceil(rand(1, p.iters) * length(p.coh));
    p.coh = p.coh(:)';
    coh = ones(ceil(p.t_max/p.fluc_t),1) * p.coh(idx);
end
momentary_coh = coh * p.dt;
s = p.sigma*sqrt(p.dt);

  % stimulus fluctuation
  % Sb .. stimulus fluctuation subtracted mean coherence level. Temporal resolution = fluc_t ms.
  % S .. real stimulus values used in DDM. Temporal resolution = 1 ms.
  % for the kernel analysis, Sb will be used.
Sb = cell(p.nfeature,1);
S = cell(p.nfeature,1);
idx = ceil((1:p.t_max)/p.fluc_t);
for n=1:p.nfeature
    Sb{n} = normrnd(momentary_coh, s, [ceil(p.t_max/p.fluc_t), p.iters]); % Sb = base
    S{n} = Sb{n}(idx,:)'; % S = S of 1 ms resolution
    S{n} = S{n} * p.k(n) + p.k0; % add sensitivity
    Sb{n} = (Sb{n} - momentary_coh)';
end

    %apply temporal weights
V = zeros(size(S{1}));
for n=1:p.nfeature
    V = V + S{n} .* repmat(p.weight, [p.iters 1]);
end
    %add the initial value
V(:,1) = V(:,1);

    %add decision noise
V = V + randn(size(V)) * p.dec_noise;
    %calculate the cumulative evidence
V = cumsum(V, 2);

    %find the choice and bound crossing time for each trial 
choice = nan(p.iters, 1);
bound_crossing = zeros(p.iters, 1);
RT = nan(p.iters, 1);
    %first find trials in which bound crossing took place
[I, J] = find(V>=repmat(p.bound_height(:,2)',[p.iters 1]) | V<=repmat(p.bound_height(:,1)',[p.iters 1]));
    %now find the first bound crossing time on these trials, ultimately we'd like to be
    %have an index to the trials with bound crossing. Also, J must be be the first time
    %of bound crossing on those trials 

    % Sort by I first, then by J to ensure minimum J comes first for each I
[~, ind] = sortrows([I, J], [1, 2]);
I = I(ind);
J = J(ind);
[I, ind] = unique(I, 'first');
J = J(ind);
bound_crossing(I) = 1;

RT(I) = J + p.non_dec_time + round(randn(size(J)) * p.non_dec_time_sd);
RT(RT<=0) = 0; % clip at 0
idx = RT >= p.t_max;
RT(idx) = NaN; % remove trials whose decision time is larger than stimulus duration or smaller than 0.
bound_crossing(idx) = 0;

L = V(sub2ind(size(V),I,J)) >= p.bound_height(J,2);
choice(I(L)) = 1;
L = V(sub2ind(size(V),I,J)) <= p.bound_height(J,1);
choice(I(L)) = 2;

    %determine what happens to decision variable and evidence traces after bound crossing 
for trial = find(bound_crossing)'
    t = ceil(RT(trial)/p.fluc_t) + 1;
    if t < 1, t = 1; end
    
    for n=1:p.nfeature
        if t <= size(Sb{n},2)
            Sb{n}(trial, t:end) = NaN;
        end
    end
end

%% save the results
sim.choice = nansum(choice==2)/sum(~isnan(choice));
sim.medianRT = nanmedian(RT);
    %average decision variable (path) aligned to stimulus onset
bt_median = floor(nanmedian(RT));
if isnan(bt_median) % no trial crossed the bound
    bt_median = p.t_max;
end

if isnan(p.bt_mediant)
    bt_mediant = ceil(bt_median/p.fluc_t);
else
    bt_mediant = p.bt_mediant;
end

    %psychophysical kernel aligned to stimulus onset
for n=1:p.nfeature
    sim.kernel{n}.motion{1} = nanmean(Sb{n}(choice==1,1:bt_mediant),1);
    sim.kernel{n}.motion{2} = nanmean(Sb{n}(choice==2,1:bt_mediant),1);
end

    %average decision variable and psychophysical kernel aligned to choice 
S_choice = repmat({nan(p.iters,bt_mediant)}, p.nfeature, 1);
for trial = 1 : p.iters
    if bound_crossing(trial)==1
        st = max(0,round(RT(trial)/p.fluc_t) - bt_mediant)+1;
        en = min(floor(p.t_max/p.fluc_t), round(RT(trial)/p.fluc_t));
        valid_portion = st:en;
        for n=1:p.nfeature
            S_choice{n}(trial,end-length(valid_portion)+1:end) = Sb{n}(trial,valid_portion);
        end
    end
end
for n=1:p.nfeature
    sim.kernel{n}.choice{1} = nanmean(S_choice{n}(choice==1,:),1);
    sim.kernel{n}.choice{2} = nanmean(S_choice{n}(choice==2,:),1);
end
sim.bound_crossed_percent = sum(bound_crossing)/length(bound_crossing) * 100;

sim.fluc_t = p.fluc_t;
sim.param = p;


t_iter = toc;
fprintf('percent cross bound: %1.1f%%, median RT = %1.3f ms, simulation time = %1.1f s\n', ...
    sim.bound_crossed_percent, nanmedian(RT), t_iter);
if sim.bound_crossed_percent < 95
    if p.error_no_reach
        error('Bound crossed percent is less than 95%. You should use longer t_max');
    else
        warning('Bound crossed percent is less than 95%. You should use longer t_max');
    end
end



