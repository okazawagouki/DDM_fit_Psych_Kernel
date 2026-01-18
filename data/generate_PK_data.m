clear;


len = 1000 * 8/75;



S = load('/Users/gouki/GoogleDrive/for__Gouki/ITintegration_psych/Analysis/IntermediateData/task2/lc_task2.mat');
D = S.trial_data;

in = D.fixed_seed == 0;

trial_data.Coh = D.Coh(in);
trial_data.RT = D.RT(in);
trial_data.Stim = D.Stim(in);
trial_data.Resp = D.Resp(in);
trial_data.Fluc = D.Coh_fluc(in);

for n=1:length(trial_data.Fluc)
    L = ceil(trial_data.RT(n) / len);
    trial_data.Fluc{n} = trial_data.Fluc{n}(:, 1:L)';
end

save example_PK_data.mat trial_data;


