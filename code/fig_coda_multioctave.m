% Figure: FORGE S-coda measures in octave bands (from S14 results, no recomputation)
load ../results/coda_multioctave.mat
f1 = figure('Position',[30 30 1500 460], 'Visible','off');
lab = {'coda level log_{10}(E_{coda}/E_{S})', 'Q_c^{-1} (single scattering)', 'rms S-pulse width (ms)'};
fld = {'CL','QI','SW'}; sc = [1 1 1e3];
for p = 1:3
    subplot(1,3,p); hold on
    for g = 1:2
        Y = res(g).(fld{p})*sc(p); Y(~res(g).OK) = NaN;
        Q = prctile(Y, [25 50 75]);
        errorbar(fc*(1 + 0.04*(g-1)), Q(2,:), Q(2,:)-Q(1,:), Q(3,:)-Q(2,:), '-o', 'LineWidth', 1.5, 'MarkerFaceColor','w');
    end
    set(gca,'XScale','log'); if p == 2, set(gca,'YScale','log'); end
    yl = ylim;
    patch([200 2000 2000 200], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');                            % resonance-affected band
    xline(36, ':', 'HandleVisibility', 'off'); xline([330 480 920], 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    ylim(yl); xlim([15 2000]); grid on; box on
    xlabel('octave-band centre frequency (Hz)'); ylabel(lab{p});
    if p == 1, legend('sensor A','sensor B','Location','southwest'); end
    panel_label(gca, p, 'tr');
end
exportgraphics(f1, '../figs/fig_coda_multioctave.png', 'Resolution', 140); close(f1);
