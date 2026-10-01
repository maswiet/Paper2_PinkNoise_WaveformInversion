% S24  (a) Observed 40-80 Hz coda/S ratio by stimulation stage (reviewer point 19; exploratory),
%          controlling for distance by a linear fit in log10 R.
%      (b) Statistics of the random fields actually simulated (reviewer point 13): sigma of lnV
%          after truncation, clipped fraction, and the radially averaged 3-D power spectrum.
maxNumCompThreads(2);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
load ../data/forge2022_events.mat
S = load(fullfile(out,'s19_H.mat')); [~, loc] = ismember(S.ev, (1:numel(ev.M))');
st = ev.stage(S.ev); rows = {};
for g = 1:2
    y = log10(S.obsratio(:,g,1)); x = log10(S.R(:,g)); ok = isfinite(y);
    X = [ones(nnz(ok),1) x(ok) st(ok)==3]; c = X \ y(ok);
    res = y(ok) - X*c; se = sqrt(sum(res.^2)/(nnz(ok)-3) * diag(inv(X'*X)));
    for s = unique(st(ok))'
        m = ok & st == s;
        rows(end+1,:) = {char('A'+g-1), s, nnz(m), median(y(m)), median(S.R(m,g))}; %#ok<SAGROW>
    end
    fprintf('sensor %c: stage-3 offset in log10(coda/S) after distance correction: %+.2f +- %.2f\n', char('A'+g-1), c(3), se(3));
end
T = cell2table(rows, 'VariableNames', {'sensor','stage','n','median_log10_coda_over_S','median_R_m'});
disp(T); writetable(T, fullfile(out,'table_stage_coda.csv'));
% ---- (b) field statistics ----
G = forge_setup(); rows = {};
spec = {'pink s.13 p1.5','P6'; 'pink s.18 p1.5','P9'; 'exp a15 s.13','E15'; 'exp a50 s.13','E50'; ...
        'Gauss a15 s.13','G15'; 'Gauss a50 s.13','G50'};
for i = 1:size(spec,1)
    M = model_catalog(spec{i,2}, G); lv = log(double(M.vp(:))/G.vp);
    clipfrac = mean(abs(lv) >= max(abs(lv))*0.999);
    rows(end+1,:) = {spec{i,1}, std(lv), clipfrac*100, min(M.vp(:)), max(M.vp(:)), min(M.vs(:)), max(M.vs(:))}; %#ok<SAGROW>
end
T2 = cell2table(rows, 'VariableNames', {'medium','sigma_lnV_after_truncation','pct_cells_at_truncation','vp_min','vp_max','vs_min','vs_max'});
disp(T2); writetable(T2, fullfile(out,'table_field_stats.csv'));
