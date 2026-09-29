function fig_seismograms(models, nev)
% FIG_SEISMOGRAMS  Observed FORGE seismograms vs model predictions.
% For each model the moment tensor of each event is fitted on P-10 ms..S+60 ms
% (direct waves only); the S coda after that window is a pure prediction of
% the medium.  Outputs:
%   figs/fig_seis_overlay.png   observed vs predicted traces, several events x models
%   figs/fig_seis_envelopes.png observed vs predicted log envelopes per event
%   figs/fig_coda_ratio.png     per-event predicted/observed coda ratio, all events
%   results/seis_compare.mat, results/table_coda_ratio.csv
if nargin < 1, models = {'H','L','A','E15','E50','G15','G50','P3','P6','P9'}; end
if nargin < 2, nev = 150; end
maxNumCompThreads(4);
here = fileparts(mfilename('fullpath'));
models = models(cellfun(@(m) exist(fullfile(here,'..','sgt',[m '.mat']),'file')==2, models));
G = forge_setup(); fs = 1/G.dt; t = (0:G.nt-1)'*G.dt - G.t0;
[b,a] = butter(3,[40 80]/(fs/2));
R = [sqrt(sum((G.evxyz - G.sens(1,:)).^2,2)) sqrt(sum((G.evxyz - G.sens(2,:)).^2,2))];
% ---------------- observed (as in S08) ----------------
ev = load(fullfile(here,'..','data','forge2022_events.mat')); ev = ev.ev;
Df = load(fullfile(here,'..','data','forge_features.mat'));
tr = (-0.05:1/ev.fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
td = (0:size(ev.W,1)-1)'/ev.fs - ev.pre; hch = {[2 3],[5 6]};
both = intersect(Df.D(1,1).k, Df.D(2,1).k);
snr = min(Df.NZ(1).snr(ismember(Df.D(1,1).k,both)), Df.NZ(2).snr(ismember(Df.D(2,1).k,both)));
[~,o] = sort(snr,'descend'); both = both(o(1:min(nev,numel(o))));
[~, ie] = ismember(both, G.ev);
OBS = zeros(G.nt, 4, numel(ie));
for i = 1:numel(ie), for g = 1:2
    x = conv2(double(ev.W(:,hch{g},both(i))), ric, 'same');
    OBS(:,2*g-1:2*g,i) = interp1(td - G.off(g), x, t, 'linear', 0);
end, end
OBS = reshape(filtfilt(b, a, reshape(OBS, G.nt, [])), G.nt, 4, []);
win = false(G.nt, 4, numel(ie));
for i = 1:numel(ie), for g = 1:2
    win(:,2*g-1:2*g,i) = repmat(t >= R(ie(i),g)/G.vp-0.010 & t <= R(ie(i),g)/G.vs+0.060, 1, 2);
end, end
basis = mt_basis();
rotf = fullfile(here,'..','data','sensor_rotation.mat');
PRED = zeros(G.nt, 4, numel(ie), numel(models), 'single'); VR = zeros(numel(ie), numel(models));
for m = 1:numel(models)
    L = load(fullfile(here,'..','sgt',[models{m} '.mat']),'E');
    GF = zeros(G.nt, 4, 6, numel(ie));
    for i = 1:numel(ie), for g = 1:2, for n = 1:2, for j = 1:6
        GF(:,2*g-2+n,j,i) = cumsum(sgt_synth(L.E{g,n}(:,:,ie(i)), basis{j}))*G.dt;
    end, end, end, end
    GF = reshape(filtfilt(b, a, reshape(GF, G.nt, [])), G.nt, 4, 6, []);
    if ~exist(rotf,'file')          % tool orientation from homogeneous model, once
        assert(strcmp(models{m},'H'), 'first model must be H to fix sensor rotation');
        ROT = cell(2,1);
        for g = 1:2
            best = -inf;
            for hand = [1 -1], for phi = 0:5:355
                Rm = [cosd(phi) sind(phi); -hand*sind(phi) hand*cosd(phi)];
                v = arrayfun(@(i) mt_invert_one(OBS(:,:,i), GF(:,:,:,i), win(:,:,i), g, Rm, fs), 1:40);
                if median(v) > best, best = median(v); ROT{g} = Rm; end
            end, end
        end
        save(rotf, 'ROT');
    end
    load(rotf, 'ROT');
    for i = 1:numel(ie)
        [VR(i,m), ~, pr, ~, sh] = mt_invert_one(OBS(:,:,i), GF(:,:,:,i), win(:,:,i), 0, ROT, fs);
        for g = 1:2, pr(:,2*g-1:2*g) = circshift(pr(:,2*g-1:2*g), sh(g)); end   % back to observed time
        PRED(:,:,i,m) = pr;
    end
    fprintf('%s: median VR %.3f\n', models{m}, median(VR(:,m)));
end
% ---------------- coda ratios (energy S+30..100 ms / S window) ----------------
co = NaN(numel(ie), 2); cp = NaN(numel(ie), 2, numel(models));
env = @(x) movmean(sum(abs(hilbert(x)).^2, 2), 20);
for i = 1:numel(ie), for g = 1:2
    tS = R(ie(i),g)/G.vs; c = 2*g-1:2*g;
    iS = t >= tS-0.008 & t < tS+0.060; iC = t >= tS+0.030 & t < tS+0.100;
    e = env(OBS(:,c,i)); co(i,g) = mean(e(iC))/max(e(iS));
    for m = 1:numel(models)
        e = env(double(PRED(:,c,i,m))); cp(i,g,m) = mean(e(iC))/max(e(iS));
    end
end, end
lr = log10(cp ./ co);                     % log10 predicted/observed coda ratio
T = table(models', squeeze(median(lr(:,1,:),1)), squeeze(median(lr(:,2,:),1)), ...
    squeeze(mean(abs(lr(:,1,:)) < 0.5, 1)), squeeze(mean(abs(lr(:,2,:)) < 0.5, 1)), median(VR)', ...
    'VariableNames', {'model','med_log10_pred_over_obs_A','med_log10_pred_over_obs_B','frac_within_x3_A','frac_within_x3_B','VR_med'});
disp(T);
writetable(T, fullfile(here,'..','results','table_coda_ratio.csv'));
save(fullfile(here,'..','results','seis_compare.mat'), 'models','ie','VR','co','cp','both','-v7.3');

% ---------------- Fig: trace overlays ----------------
show = {'H','A','E15','G15','P6'}; show = show(ismember(show, models));
lbl = struct('H','homogeneous','L','layered','A','VTI','E15','exponential a=15 m','E50','exponential a=50 m', ...
    'G15','Gaussian a=15 m','G50','Gaussian a=50 m','P3','pink \sigma=0.09','P6','pink \sigma=0.13','P9','pink \sigma=0.18');
sel = pick_events(R(ie,:), 4);
save(fullfile(here,'..','results','seis_examples.mat'), 'sel', 'show', 'models', 'ie', 'VR', 't', 'R');
f1 = figure('Position',[30 30 1700 950], 'Visible','off'); np = 0;
for r = 1:numel(sel)
    i = sel(r); g = 1; c = 1; tS = R(ie(i),g)/G.vs; tP = R(ie(i),g)/G.vp;
    w = t > tP-0.02 & t < tS+0.14;
    o = OBS(w,c,i); sc = max(abs(o));
    for q = 1:numel(show)
        m = find(strcmp(models, show{q}));
        subplot(numel(sel), numel(show), (r-1)*numel(show)+q); hold on
        patch([tS+0.03 tS+0.10 tS+0.10 tS+0.03]-tP, [-1.2 -1.2 1.2 1.2], [0.93 0.93 0.93], 'EdgeColor','none');
        plot(t(w)-tP, o/sc, 'k', 'LineWidth', 0.9);
        plot(t(w)-tP, double(PRED(w,c,i,m))/sc, 'r', 'LineWidth', 0.9);
        xline(tS-tP,':'); ylim([-1.2 1.2]); xlim([-0.02 tS-tP+0.14]); box on
        np = np + 1; panel_label(gca, np);
        if q == 1, ylabel(sprintf('ev %d, R=%.0f m', G.ev(ie(i)), R(ie(i),g))); end
        if r == numel(sel), xlabel('time after P (s)'); end
        text(0.98, 0.9, sprintf('VR %.2f', VR(i,m)), 'Units','normalized','HorizontalAlignment','right','FontSize',8);
    end
end
exportgraphics(f1, fullfile(here,'..','figs','fig_seis_overlay.png'), 'Resolution', 140); close(f1);

% ---------------- Fig: envelopes per event ----------------
f2 = figure('Position',[30 30 1600 800], 'Visible','off');
cols = lines(numel(show));
for r = 1:numel(sel)
    for g = 1:2
        i = sel(r); c = 2*g-1:2*g; tS = R(ie(i),g)/G.vs;
        subplot(2, numel(sel), (g-1)*numel(sel)+r);
        e = env(OBS(:,c,i)); pk = max(e(t>tS-0.01 & t<tS+0.06));
        semilogy(t-tS, e/pk, 'k', 'LineWidth', 2.2); hold on
        for q = 1:numel(show)
            m = find(strcmp(models, show{q}));
            e = env(double(PRED(:,c,i,m))); pk = max(e(t>tS-0.01 & t<tS+0.06));
            semilogy(t-tS, e/pk, 'Color', cols(q,:), 'LineWidth', 1.1);
        end
        xlim([-0.06 0.13]); ylim([1e-4 2]); grid on
        panel_label(gca, (g-1)*numel(sel)+r, 'tr');
        if r == 1, ylabel('E / E_{S peak}'); end
        xlabel('t - t_S (s)');
    end
end
legend(['observed', cellfun(@(s) lbl.(s), show, 'uni', 0)], 'Location','southwest','FontSize',7);
exportgraphics(f2, fullfile(here,'..','figs','fig_seis_envelopes.png'), 'Resolution', 140); close(f2);

% per-event coda-ratio figure: fig_coda_ratio.m (reads seis_compare.mat)
end

function sel = pick_events(R, n)
% events spread in distance to sensor A
[~, o] = sort(R(:,1)); sel = o(round(linspace(1, numel(o), n+2))); sel = sel(2:end-1);
end
