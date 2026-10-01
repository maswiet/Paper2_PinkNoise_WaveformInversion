% S28  Robustness of the per-event coda test to location / S-timing error (reviewer points 11, 26):
% the observed 40-80 Hz S-envelope peak is compared with the predicted S time; the coda-prediction
% statistics are recomputed for events whose observed S peak lies within +-10 ms (and +-5 ms)
% of the prediction on BOTH sensors (a well-timed subset).
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
G = forge_setup(); fs = 1/G.dt; t = (0:G.nt-1)'*G.dt - G.t0; [b,a] = butter(3,[40 80]/(fs/2));
ev = load(fullfile(here,'..','data','forge2022_events.mat')); ev = ev.ev;
S0 = load(fullfile(out,'s19_H.mat'));
tr = (-0.05:1/ev.fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
td = (0:size(ev.W,1)-1)'/ev.fs - ev.pre; hch = {[2 3],[5 6]};
ne = numel(S0.ev); dS = NaN(ne,2);
for i = 1:ne
    for g = 1:2
        x = filtfilt(b, a, interp1(td - G.off(g), conv2(double(ev.W(:,hch{g},S0.ev(i))), ric, 'same'), t, 'linear', 0));
        e = sum(abs(hilbert(x)).^2, 2); tS = S0.R(i,g)/G.vs;
        m = t >= tS - 0.025 & t <= tS + 0.060; [~, j] = max(e .* m); dS(i,g) = t(j) - tS;
    end
end
fprintf('observed S-peak minus predicted S: median %.1f / %.1f ms, |.|<10 ms: %.0f%% / %.0f%%\n', ...
    1e3*median(dS), 100*mean(abs(dS) < 0.010));
w10 = all(abs(dS) < 0.010, 2); w5 = all(abs(dS) < 0.005, 2);
d = dir(fullfile(out, 's19_*.mat')); rows = {};
for k = 1:numel(d)
    S = load(fullfile(out, d(k).name));
    for sub = {'all', 'within10ms', 'within5ms'}
        switch sub{1}, case 'all', m = true(ne,1); case 'within10ms', m = w10; otherwise, m = w5; end
        rows(end+1,:) = {S.model, sub{1}, nnz(m), median(S.lr(m,1,1,1),'omitnan'), median(S.lr(m,2,1,1),'omitnan')}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'model','subset','n','med_log10_pred_obs_A','med_log10_pred_obs_B'});
writetable(T, fullfile(out, 'table_s_timing_subsets.csv'));
disp(T(ismember(T.model, {'H','L','A','P1','P3','P6','E50','G50'}), :));
save(fullfile(out, 's_timing.mat'), 'dS', 'w10', 'w5');
