% S29  Effective scattering strength sigma_eff per family and sensor from the per-event coda test
% (reviewer point 3): sigma at which the median log10(predicted/observed coda ratio) crosses zero,
% interpolated linearly in log sigma between simulated levels (realisations averaged).
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
fam = struct( ...
 'pink15', {{0.02,{'P4'}; 0.045,{'P1','R1','R2','T1'}; 0.09,{'P3','R3'}; 0.13,{'P6','R6'}; 0.18,{'P9','R5'}}}, ...
 'pink18', {{0.045,{'P2'}; 0.09,{'P8','X_pink_0.09_1.8_33'}; 0.13,{'X_pink_0.13_1.8_11','X_pink_0.13_1.8_33'}; 0.18,{'P11','X_pink_0.18_1.8_33'}}}, ...
 'exp15', {{0.09,{'X_exp_0.09_15_11','X_exp_0.09_15_33'}; 0.13,{'E15','E15b'}; 0.18,{'X_exp_0.18_15_11','X_exp_0.18_15_33'}}}, ...
 'exp50', {{0.09,{'X_exp_0.09_50_11','X_exp_0.09_50_33'}; 0.13,{'E50','X_exp_0.13_50_33'}; 0.18,{'X_exp_0.18_50_11','X_exp_0.18_50_33'}}}, ...
 'gau15', {{0.09,{'X_gau_0.09_15_11','X_gau_0.09_15_33'}; 0.13,{'G15','X_gau_0.13_15_33'}; 0.18,{'X_gau_0.18_15_11','X_gau_0.18_15_33'}}}, ...
 'gau50', {{0.09,{'X_gau_0.09_50_11','X_gau_0.09_50_33'}; 0.13,{'G50','G50b'}; 0.18,{'X_gau_0.18_50_11','X_gau_0.18_50_33'}}});
need = {};
for f = fieldnames(fam)'
    L = fam.(f{1}); for i = 1:size(L,1), need = [need L{i,2}]; end %#ok<AGROW>
end
need = unique(need); miss = need(~cellfun(@(n) exist(fullfile(out, ['s19_' n '.mat']), 'file') == 2, need));
miss = miss(cellfun(@(n) exist(fullfile(here,'..','sgt',[n '.mat']), 'file') == 2, miss));
if ~isempty(miss), s19_coda_prediction(miss, 12); end
rows = {}; curves = {};
for f = fieldnames(fam)'
    L = fam.(f{1}); sg = cell2mat(L(:,1)); med = NaN(numel(sg), 2); nrel = zeros(numel(sg),1);
    for i = 1:numel(sg)
        v = [];
        for n = L{i,2}
            fn = fullfile(out, ['s19_' n{1} '.mat']);
            if exist(fn, 'file'), S = load(fn); v(end+1,:) = [median(S.lr(:,1,1,1),'omitnan') median(S.lr(:,2,1,1),'omitnan')]; end %#ok<AGROW>
        end
        if ~isempty(v), med(i,:) = mean(v, 1); nrel(i) = size(v,1); end
        curves(end+1,:) = {f{1}, sg(i), nrel(i), med(i,1), med(i,2)}; %#ok<SAGROW>
    end
    se = NaN(1,2);
    for g = 1:2
        y = med(:,g); ok = isfinite(y); x = log(sg(ok)); y = y(ok);
        j = find(y(1:end-1) <= 0 & y(2:end) > 0, 1);
        if ~isempty(j), se(g) = exp(x(j) + (0 - y(j))*(x(j+1) - x(j))/(y(j+1) - y(j)));
        elseif all(y > 0), se(g) = -min(sg);                   % below the smallest simulated sigma
        elseif all(y < 0), se(g) = -max(sg) - 1;               % above the largest (flag)
        end
    end
    rows(end+1,:) = {f{1}, se(1), se(2)}; %#ok<SAGROW>
end
T = cell2table(rows, 'VariableNames', {'family','sigma_eff_A','sigma_eff_B'});
C = cell2table(curves, 'VariableNames', {'family','sigma','n_realisations','med_log10_pred_obs_A','med_log10_pred_obs_B'});
disp(C); disp(T);
writetable(T, fullfile(out, 'table_sigma_eff.csv')); writetable(C, fullfile(out, 'table_sigma_curves.csv'));
