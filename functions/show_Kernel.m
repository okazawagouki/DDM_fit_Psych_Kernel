function fh = show_Kernel(kernelsim, trial_data, opt)
def.max_coh = 0; % maximum coherence of trials used for kernel analysis
def.smoothing = 1; % smoothing of kernel, e.g. fspecial('average', [1, 3])
def.fluc_t = .1; % duration of each fluctuation (sec)
def.positive_cat = 2; % category (choice) of positive coherence
def.nfeature = 3;
def.ylim = [-0.01 0.1];

opt = safeStructAssign(def, opt);

coh = trial_data.Coh;
fluc = trial_data.Fluc;
resp = trial_data.Resp;
rt = trial_data.RT / 1000; % convert from msec to sec

if ~isempty(opt.max_coh)
    ind = abs(coh) <= opt.max_coh;
    coh = coh(ind);
    fluc = fluc(ind);
    resp = resp(ind);
    rt = rt(ind);
    fprintf('Using coherence:');
    fprintf('%g, ', unique(abs(coh)));
    fprintf('\n');
end

%% data
framelen = ceil(median(rt)/opt.fluc_t);

nfrm = cellfun(@(x)size(x,1), fluc);
if any(nfrm * opt.fluc_t < rt)
    error('some stimulus fluctuations are shorter than RT');
end


stkernel = cell(opt.nfeature,1);
stkernel_se = cell(opt.nfeature,1);
rekernel = cell(opt.nfeature,1);
rekernel_se = cell(opt.nfeature,1);

for f=1:opt.nfeature
    stfluc = nan(length(fluc), framelen);
    refluc = nan(length(fluc), framelen);
    for n=1:length(fluc)
        rtf = ceil(rt(n)/opt.fluc_t); % frame num corresponding to RT
        ln = min(rtf, framelen); % length of frame used
        if opt.smoothing>1
            fluc{n}(:,f) = nanconv(fluc{n}(:,f), fspecial('average', [opt.smoothing, 1]), 'same');
        end
        stfluc(n, 1:ln) = fluc{n}(1:ln,f) - coh(n);
        refluc(n, end-ln+1:end) = fluc{n}(rtf-ln+1:rtf,f) - coh(n);
    end

    stkernel{f} = nanmean(stfluc(resp == opt.positive_cat,:)) - nanmean(stfluc(resp == 3 - opt.positive_cat,:));
    stkernel_se{f} = calc_diff_SE(stfluc(resp==1,:), stfluc(resp==2,:));
    rekernel{f} = nanmean(refluc(resp == opt.positive_cat,:)) - nanmean(refluc(resp == 3 - opt.positive_cat,:));
    rekernel_se{f} = calc_diff_SE(refluc(resp==1,:), refluc(resp==2,:));
end


t = (0:framelen-1) * opt.fluc_t;

%% model

S = load(kernelsim);
skernel = S.sim.kernel;

sframelen = min(framelen, length(skernel{1}.motion{1}));

stkernel_sim = cell(opt.nfeature,1);
rekernel_sim = cell(opt.nfeature,1);
for f=1:opt.nfeature
    stkernel_sim{f} = skernel{f}.motion{1}(1:sframelen) - skernel{f}.motion{2}(1:sframelen);
    rekernel_sim{f} = skernel{f}.choice{1}(end-sframelen+1:end) - skernel{f}.choice{2}(end-sframelen+1:end);
    if opt.smoothing>1
        stkernel_sim{f}(:) = nanconv(stkernel_sim{f}(:), fspecial('average', [opt.smoothing, 1]), 'same');
        rekernel_sim{f}(:) = nanconv(rekernel_sim{f}(:), fspecial('average', [opt.smoothing, 1]), 'same');
    end
end
simt = (0:sframelen-1) * opt.fluc_t;

%% figure
fh = figure('color', 'w', 'position', [100 100 300 150 * opt.nfeature]);
ax = nan(opt.nfeature, 2);
for f=1:opt.nfeature
    ax(f, 1) = subplot(opt.nfeature, 2, f * 2 - 1);
    hold on;
    plot_trace(t, stkernel{f}, stkernel_se{f}, [1 0 0]);
    plot(simt, stkernel_sim{f}, '-', 'color', [.3 .3 .3], 'linew', 3);
    axis square;
    plot(xlim, [0 0], 'color', [.5 .5 .5]);
    ylim(opt.ylim)
    if f == opt.nfeature
        xlabel('Time from stim on');
    end
    
    ax(f, 2) = subplot(opt.nfeature, 2, f * 2);
    hold on;
    plot_trace(fliplr(-t), rekernel{f}, rekernel_se{f}, [1 0 0]);
    plot(fliplr(-simt), rekernel_sim{f}, '-', 'color', [.3 .3 .3], 'linew', 3);
    axis square;
    plot(xlim, [0 0], 'color', [.5 .5 .5]);
    ylim(opt.ylim)
    if f == opt.nfeature
        xlabel('Time from response');
    end
    set(gca, 'YAxisLocation', 'right');
end


end



function C = nanconv(A,B,shape)

    if (nargin < 3)
      shape = 'full';
    end

    AO = ones(size(A));
    I = isnan(A);
    A(I) = 0;
    AO(I) = 0;
    C = convn(A,B,shape) ./ convn(AO,B,shape);
    C(I) = nan;

end


function se = calc_diff_SE(data1, data2)
    se1 = nanstd(data1,0,1) ./ sqrt(sum(~isnan(data1),1));
    se2 = nanstd(data2,0,1) ./ sqrt(sum(~isnan(data2),1));
    se = sqrt(se1.^2 + se2.^2);
end


function lh = plot_trace(x, y, se, color)
    hold on;


    if isempty(se)
        se = zeros(size(y));
    end

    palecol = color * .5 + [1 1 1] * .5;
    x2 = [x(:); flipud(x(:))];
    y2 = [y(:)-se(:); flipud(y(:)+se(:))];
    h = fill(x2,y2,palecol);
    set(h, 'EdgeColor', palecol);

        %draw main trace
    lh = plot(x,y,'Color', color, 'LineWidth', 1);
end

