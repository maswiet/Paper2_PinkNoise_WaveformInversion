% S27  Summaries of the per-event coda-prediction test (S19) for all models.
% results/table_coda_prediction.csv      primary windows, all events and top-150 subset
% results/table_coda_window_sensitivity.csv  fit end x coda window grid (median log10 pred/obs)
% results/table_mt_variants.csv          full / deviatoric / DC, condition number, decomposition
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
d = dir(fullfile(out, 's19_*.mat')); x3 = log10(3);
rows = {}; rw = {}; rm = {};
for i = 1:numel(d)
    S = load(fullfile(out, d(i).name)); lr = S.lr;
    for sub = {'all','top150'}
        if strcmp(sub{1},'all'), m = true(size(lr,1),1); else, m = S.top150; end
        a = lr(m,1,1,1); b = lr(m,2,1,1);
        rows(end+1,:) = {S.model, sub{1}, nnz(m), median(a,'omitnan'), median(b,'omitnan'), ...
            100*mean(abs(a) < x3,'omitnan'), 100*mean(abs(b) < x3,'omitnan'), median(S.vr(m,1),'omitnan')}; %#ok<SAGROW>
    end
    for f = 1:numel(S.fitEnd)
        for q = 1:size(S.coda,1)
            rw(end+1,:) = {S.model, 1e3*S.fitEnd(f), sprintf('%d-%d', round(1e3*S.coda(q,:))), ...
                S.coda(q,1) < S.fitEnd(f), median(lr(:,1,f,q),'omitnan'), median(lr(:,2,f,q),'omitnan')}; %#ok<SAGROW>
        end
    end
    vn = {'full','dev','dc'};
    for v = 1:3
        rm(end+1,:) = {S.model, vn{v}, median(S.lrv(:,1,v),'omitnan'), median(S.lrv(:,2,v),'omitnan'), ...
            median(S.vrv(:,v),'omitnan'), median(S.cond,'omitnan'), prctile(S.cond,90), ...
            median(abs(S.decomp(:,1)),'omitnan'), median(abs(S.decomp(:,2)),'omitnan'), median(S.decomp(:,3),'omitnan')}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'model','events','n','med_log10_pred_obs_A','med_log10_pred_obs_B','pct_within3_A','pct_within3_B','VR'});
writetable(T, fullfile(out, 'table_coda_prediction.csv')); disp(T(strcmp(T.events,'all'),:));
W = cell2table(rw, 'VariableNames', {'model','fit_end_ms','coda_ms','overlap','med_A','med_B'});
writetable(W, fullfile(out, 'table_coda_window_sensitivity.csv'));
V = cell2table(rm, 'VariableNames', {'model','mt','med_A','med_B','VR','cond_median','cond_p90','pct_ISO_abs','pct_CLVD_abs','pct_DC'});
writetable(V, fullfile(out, 'table_mt_variants.csv')); disp(V(ismember(V.model, {'H','P3','P6'}),:));
disp(W(ismember(W.model, {'H','P3'}),:));
