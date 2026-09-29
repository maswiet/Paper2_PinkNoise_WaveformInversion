% S12  Equal-sigma shape-resolution test (truths R6 pink, E15b exponential,
% G50b Gaussian; candidates at sigma 0.13) -> results/table_shape_test.csv
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
truths = {'R6','E15b','G50b'}; cand = {'P6','E15','E50','G15','G50','H','A'};
rows = {};
for i = 1:numel(truths)
    C = load(fullfile(out, sprintf('compare_%s_forge.mat', truths{i}))); R = C.RES;
    S = readtable(fullfile(out, sprintf('table_shape_%s.csv', truths{i})));
    for m = 1:numel(cand)
        j = find(strcmp({R.name}, cand{m})); k = find(strcmp(S.model, cand{m}));
        rows(end+1,:) = {truths{i}, cand{m}, R(j).PHI, S.D2_rms(k), S.D1_misfit(k)}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'truth','model','PHI','D2_rms','D1_misfit'});
writetable(T, fullfile(out,'table_shape_test.csv')); disp(T);
% coherence curves and D1 medians for data and models
ms = {'data','H','A','P3','P6','P9','E15','E50','G15','G50'}; rows = {};
for m = 1:numel(ms)
    Z = load(fullfile(out, sprintf('shape_%s_forge.mat', ms{m})));
    rows(end+1,:) = {ms{m}, median(Z.dlev,'omitnan'), Z.coh(1,1), Z.coh(1,end), Z.coh(2,1), Z.coh(2,end)}; %#ok<SAGROW>
end
T = cell2table(rows, 'VariableNames', {'model','D1_median','cohA_0_5m','cohA_80_160m','cohB_0_5m','cohB_80_160m'});
writetable(T, fullfile(out,'table_coherence.csv')); disp(T);
