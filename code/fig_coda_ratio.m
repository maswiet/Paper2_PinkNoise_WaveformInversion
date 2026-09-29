% Per-event predicted/observed S-coda ratio (from fig_seismograms results)
S = load('../results/seis_compare.mat');
lr = log10(S.cp ./ S.co); models = S.models;
lab = struct('H','homogeneous','L','layered','A','VTI','E15','exponential a=15','E50','exponential a=50', ...
    'G15','Gaussian a=15','G50','Gaussian a=50','P3','pink \sigma=.09','P6','pink \sigma=.13','P9','pink \sigma=.18');
f = figure('Position',[30 30 1150 500], 'Visible','off');
for g = 1:2
    subplot(1,2,g); hold on
    for m = 1:numel(models)
        v = lr(:,g,m); v = v(isfinite(v)); q = prctile(v, [10 25 50 75 90]);
        if startsWith(models{m},'P'), col = [0.85 0.33 0.10];
        elseif startsWith(models{m},{'E','G'}), col = [0.2 0.45 0.75];
        else, col = [0.35 0.35 0.35]; end
        plot([q(1) q(5)], [m m], '-', 'Color', col);
        plot([q(2) q(4)], [m m], '-', 'Color', col, 'LineWidth', 7);
        plot(q(3), m, 'kd', 'MarkerFaceColor', 'w');
    end
    xline(0,'k--'); xline([-log10(3) log10(3)],'k:');
    set(gca,'YTick',1:numel(models),'YTickLabel',cellfun(@(s) lab.(s), models,'uni',0),'YDir','reverse');
    xlabel('log_{10}(predicted / observed coda ratio), per event'); box on; xlim([-2.5 1.5]); ylim([0.4 numel(models)+0.6]);
    panel_label(gca, g, 'tr');
end
exportgraphics(f, '../figs/fig_coda_ratio.png', 'Resolution', 150); close(f);
disp('done');
