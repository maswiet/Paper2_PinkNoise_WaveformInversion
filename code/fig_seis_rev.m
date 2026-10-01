function fig_seis_rev(models)
% FIG_SEIS_REV  Observed vs predicted sensor-A seismograms (40-80 Hz) with the revised,
% non-overlapping windows: moment tensor fitted on [P-10 ms, S+15 ms] (light blue), guard band,
% coda window [S+40, S+100 ms] (grey). Three events (near, middle, far) x numel(models) columns.
if nargin < 1, models = {'H','L6','P1','P3'}; end
maxNumCompThreads(4);
here = fileparts(mfilename('fullpath'));
G = forge_setup(); fs = 1/G.dt; t = (0:G.nt-1)'*G.dt - G.t0;
[b,a] = butter(3,[40 80]/(fs/2));
ev = load(fullfile(here,'..','data','forge2022_events.mat')); ev = ev.ev;
S19 = load(fullfile(here,'..','results','s19_H.mat'));
tr = (-0.05:1/ev.fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
td = (0:size(ev.W,1)-1)'/ev.fs - ev.pre; hch = {[2 3],[5 6]};
[~, ie] = ismember(S19.ev, G.ev); RA = S19.R(:,1);
[~, o] = sort(RA); pick = o(round([0.2 0.5 0.8]*numel(o)));
L = load(fullfile(here,'..','data','sensor_rotation.mat')); ROT = L.ROT; basis = mt_basis();
f1 = figure('Position',[30 30 1650 820], 'Visible','off'); np = 0;
for r = 1:3
    i = pick(r); k = ie(i);
    obs = zeros(G.nt, 4);
    for g = 1:2
        x = conv2(double(ev.W(:,hch{g},S19.ev(i))), ric, 'same');
        obs(:,2*g-1:2*g) = interp1(td - G.off(g), x, t, 'linear', 0);
    end
    obs = filtfilt(b, a, obs);
    tP = S19.R(i,:)/G.vp; tS = S19.R(i,:)/G.vs;
    win = false(G.nt,4); for g = 1:2, win(:,2*g-1:2*g) = repmat(t >= tP(g)-0.010 & t <= tS(g)+0.015, 1, 2); end
    for q = 1:numel(models)
        Lm = load(fullfile(here,'..','sgt',[models{q} '.mat']),'E');
        gf = zeros(G.nt,4,6);
        for g = 1:2, for n = 1:2, for j = 1:6
            gf(:,2*g-2+n,j) = cumsum(sgt_synth(Lm.E{g,n}(:,:,k), basis{j}))*G.dt;
        end, end, end
        gf = reshape(filtfilt(b, a, reshape(gf, G.nt, [])), G.nt, 4, 6);
        Rf = mt_fit_variants(obs, gf, win, ROT, fs, {'full'}, []);
        o1 = Rf.full.obsS(:,1); p1 = Rf.full.pred(:,1); sc = max(abs(o1(t > tP(1)-0.02 & t < tS(1)+0.14)));
        np = np + 1; subplot(3, numel(models), np); hold on
        patch(([tP(1)-0.010 tS(1)+0.015 tS(1)+0.015 tP(1)-0.010]) - tP(1), [-1.3 -1.3 1.3 1.3], [0.85 0.92 1], 'EdgeColor','none');
        patch(([tS(1)+0.040 tS(1)+0.100 tS(1)+0.100 tS(1)+0.040]) - tP(1), [-1.3 -1.3 1.3 1.3], [0.9 0.9 0.9], 'EdgeColor','none');
        plot(t - tP(1), o1/sc, 'k', 'LineWidth', 0.8); plot(t - tP(1), p1/sc, 'r', 'LineWidth', 0.8);
        xline(tS(1)-tP(1), ':'); xlim([-0.02 tS(1)-tP(1)+0.13]); ylim([-1.3 1.3]); box on
        if q == 1, ylabel(sprintf('R = %.0f m', S19.R(i,1))); end
        if r == 3, xlabel('time after P (s)'); end
        text(0.98, 0.93, sprintf('VR %.2f', Rf.full.vr), 'Units','normalized', 'HorizontalAlignment','right', 'FontSize', 8);
        panel_label(gca, np);
    end
end
exportgraphics(f1, fullfile(here,'..','figs','fig_seis_rev.png'), 'Resolution', 140); close(f1);
end
