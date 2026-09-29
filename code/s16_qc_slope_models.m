% S16  Qc^-1 frequency slope (25 -> 50 Hz, sensor A) for FORGE and all simulated
% realisations, grouped by family; figure + CSV for the manuscript.
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
D = load(fullfile(out,'coda_multioctave.mat'));
q = D.res(1).QI(:,1:2); q(~D.res(1).OK(:,1:2)) = NaN;
dq = median(q, 'omitnan'); dslope = log2(dq(2)/dq(1));
% bootstrap over events for the data slope
rng(5); bs = zeros(1000,1); n = size(q,1);
for b = 1:1000, i = randi(n, n, 1); m = median(q(i,:), 'omitnan'); bs(b) = log2(m(2)/m(1)); end
fam = {'smooth', {'H','L','A'}; 'pink s=.09-.18', {'P3','P6','R6','R6c','R6d','P9'}; ...
       'exponential a=15', {'E15','E15b','E15c','E15d'}; 'exponential a=50', {'E50'}; ...
       'Gaussian a=15', {'G15'}; 'Gaussian a=50', {'G50','G50b','G50c','G50d'}};
rows = {};
for f = 1:size(fam,1)
    for m = fam{f,2}
        Z = load(fullfile(out, sprintf('codamo_%s_5m.mat', m{1})));
        mq = median(Z.QI{1}, 'omitnan');
        rows(end+1,:) = {fam{f,1}, m{1}, mq(1), mq(2), log2(mq(2)/mq(1))}; %#ok<SAGROW>
    end
end
rows(end+1,:) = {'FORGE data', 'data', dq(1), dq(2), dslope};
T = cell2table(rows, 'VariableNames', {'family','model','Qc_inv_25Hz','Qc_inv_50Hz','slope_25_50'});
T.slope_lo95 = NaN(height(T),1); T.slope_hi95 = T.slope_lo95;
T.slope_lo95(end) = prctile(bs, 2.5); T.slope_hi95(end) = prctile(bs, 97.5);
disp(T); writetable(T, fullfile(out,'table_qc_slope.csv'));
% figure
f1 = figure('Position',[40 40 820 420], 'Visible','off'); hold on
col = [0.35 0.35 0.35; 0.85 0.33 0.10; 0.2 0.45 0.75; 0.45 0.65 0.9; 0.3 0.6 0.3; 0.55 0.8 0.5];
for fi = 1:size(fam,1)
    s = T.slope_25_50(strcmp(T.family, fam{fi,1}));
    plot(s, fi*ones(size(s)), 'o', 'MarkerSize', 8, 'MarkerFaceColor', col(fi,:), 'MarkerEdgeColor', 'k');
end
patch([T.slope_lo95(end) T.slope_hi95(end) T.slope_hi95(end) T.slope_lo95(end)], [0.4 0.4 size(fam,1)+0.6 size(fam,1)+0.6], ...
    [0 0 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
xline(dslope, 'k-', 'LineWidth', 2);
set(gca, 'YTick', 1:size(fam,1), 'YTickLabel', fam(:,1), 'YDir', 'reverse'); ylim([0.4 size(fam,1)+0.6]); box on; grid on
xlabel('Q_c^{-1} frequency slope, log_2[Q_c^{-1}(50 Hz)/Q_c^{-1}(25 Hz)]');
exportgraphics(f1, fullfile(here,'..','figs','fig_qc_slope.png'), 'Resolution', 150); close(f1);
