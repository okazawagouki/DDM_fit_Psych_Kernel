function show_Psych_Chrono(fname)


S = load(fname);


Coh = S.trial_data.Coh * 1e2;
Cor = S.trial_data.Resp == 2;
RT = S.trial_data.RT;
mCor = 1 - S.trial_data.fit_pchoice1;
mRT = S.trial_data.fit_rt;

uCoh = unique(Coh);

c0 = 1.5;
uCoh_mod = uCoh;
uCoh_mod(uCoh_mod==0) = c0; % replace 0% for log-scale plot

h =  findobj('type','figure');
n = length(h);
figure('color', 'w', 'position', [100 + n * 50, 300 - n * 50, 400, 250]);

subplot(1,2, 1);
hold on;
M = calcGroupMean(mCor, Coh, uCoh, 'binary');
plot(uCoh_mod, M, 'color', [.5 .5 .5], 'linew', 2);

[M, SE] = calcGroupMean(Cor, Coh, uCoh, 'binary');
errorbar_mod(uCoh_mod, M, SE);
axis square;
xlabel('Coherence (%)');
ylabel('P(response=2)');

subplot(1,2, 2);
hold on;
M = calcGroupMean(mRT, Coh, uCoh);
plot(uCoh_mod, M, 'color', [.5 .5 .5], 'linew', 2);

[M, SE] = calcGroupMean(RT, Coh, uCoh);
errorbar_mod(uCoh_mod, M, SE);

axis square;
xlabel('Coherence (%)');
ylabel('RT (ms)');


end


function [m, mse] = calcGroupMean(v, g, groups, data_type, mn_type)
    if nargin<4 || isempty(data_type)
        data_type = 'continuous';
    end
    if nargin<5 || isempty(mn_type)
        mn_type = 'mean';
    end

    switch mn_type
        case 'mean'
            m = arrayfun(@(s) mean(v(g==s)), groups);
        case 'median'
            m = arrayfun(@(s) median(v(g==s)), groups);
    end
    switch data_type
        case 'binary'
            mse = arrayfun(@(s) sqrt(m(groups==s).*(1-m(groups==s))./sum(g==s)), groups);
        case 'continuous'
            mse = arrayfun(@(s) std(v(g==s))/sqrt(sum(g==s)), groups);
    end
end

function errorbar_mod(x,y,se, col)
    if ~exist('col', 'var')
        col = 'k';
    end
    plot(x,y, '.', 'color', col, 'markersize', 12);
    for n=1:length(x)
        plot([1 1] * x(n), [-1 1] * se(n) + y(n), 'color', col);
    end
end


