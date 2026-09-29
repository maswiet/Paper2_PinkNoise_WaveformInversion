function T = fig_shape_diagnostics(obs, models, tag)
% FIG_SHAPE_DIAGNOSTICS  Compare shape-sensitive diagnostics (S11) of an
% observation ('data' or a synthetic truth) with candidate models.
% Misfit: D1 = |median difference|; D2 = rms difference of median coherence curves (both sensors).
if nargin < 3, tag = obs; end
here = fileparts(mfilename('fullpath'));
ld = @(n) load(fullfile(here,'..','results',sprintf('shape_%s_forge.mat', n)));
O = ld(obs); n = numel(models);
d1 = NaN(n,1); d2 = NaN(n,1);
for m = 1:n
    Z = ld(models{m});
    d1(m) = abs(median(Z.dlev,'omitnan') - median(O.dlev,'omitnan'));
    v = Z.coh - O.coh; d2(m) = sqrt(mean(v(:).^2,'omitnan'));
end
T = table(models(:), d1, d2, d1/0.5 + d2/0.1, 'VariableNames', {'model','D1_misfit','D2_rms','combined'});
T = sortrows(T, 'combined');
fprintf('\nObserved: %s\n', obs); disp(T);
writetable(T, fullfile(here,'..','results',sprintf('table_shape_%s.csv', tag)));
% ---- figure ----
col = @(s) [0.35 0.35 0.35]*any(strcmp(s,{'H','L','A'})) + [0.2 0.45 0.75]*any(s(1)=='EG') + ...
           [0.85 0.33 0.10]*any(s(1)=='PRT');
f = figure('Position',[40 40 1500 460], 'Visible','off');
subplot(1,3,1); hold on
for m = 1:n
    Z = ld(models{m}); q = prctile(Z.dlev, [25 50 75]);
    plot([q(1) q(3)], [m m], '-', 'Color', col(models{m}), 'LineWidth', 7); plot(q(2), m, 'kd', 'MarkerFaceColor','w');
end
q = prctile(O.dlev, [25 50 75]); patch([q(1) q(3) q(3) q(1)], [0.4 0.4 n+0.6 n+0.6], [0 0 0], 'FaceAlpha', 0.08, 'EdgeColor','none');
xline(q(2), 'k-', 'LineWidth', 2);
set(gca,'YTick',1:n,'YTickLabel',models,'YDir','reverse'); ylim([0.4 n+0.6]); box on
xlabel('coda level 40-80 Hz minus 20-40 Hz (log_{10})'); panel_label(gca, 1, 'tr');
c = (O.edges(1:end-1) + O.edges(2:end))/2; c(1) = 3;
for g = 1:2
    subplot(1,3,1+g); hold on
    for m = 1:n
        Z = ld(models{m}); plot(c, Z.coh(g,:), '-o', 'Color', col(models{m}), 'MarkerSize', 3);
    end
    plot(c, O.coh(g,:), 'k-s', 'LineWidth', 2.5, 'MarkerFaceColor','k');
    set(gca,'XScale','log'); grid on; ylim([0.3 1]); xlabel('inter-event separation (m)'); ylabel('median coda coherence');
    panel_label(gca, 1+g, 'tr');
end
exportgraphics(f, fullfile(here,'..','figs',sprintf('fig_shape_%s.png', tag)), 'Resolution', 140); close(f);
end
