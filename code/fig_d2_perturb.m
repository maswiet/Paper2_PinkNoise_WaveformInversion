% FIG_D2_PERTURB  Additional-file figure of the D2 location stress test (window S+40..100 ms), drawn
% from results/table_d2_perturb.csv (written by S31): perturbation 5 m (solid), 20 m (dashed), 50 m (dotted).
here = fileparts(mfilename('fullpath'));
TP = readtable(fullfile(here,'..','results','table_d2_perturb.csv'));
names = {'data','H','L','GRAD3','A','L6','FZ','P3','P6','E15','G50'};
cls = containers.Map(names, {'data','smooth','smooth','smooth','smooth','det','det','pink','pink','exp','gau'});
col = containers.Map({'data','smooth','det','pink','exp','gau'}, {[0 0 0],[.6 .6 .6],[.55 .25 .65],[.85 .33 .1],[.2 .45 .75],[.3 .6 .3]});
edges = [0 5 10 20 40 80 160]; cc = (edges(1:end-1) + edges(2:end))/2; cc(1) = 3;
ls = {'-', '--', ':'}; SP = [5 20 50];
f = figure('Visible','off');
for g = 1:2
    ax = subplot(1,2,g); hold on
    for m = [2:numel(names) 1]
        for k = 1:3
            r = TP(strcmp(TP.model, names{m}) & strcmp(TP.window, 'S+40-100 ms') & strcmp(TP.sensor, char('A'+g-1)) & TP.perturb_m == SP(k), :);
            if isempty(r), continue; end
            plot(cc, r.median_coh, ls{k}, 'Color', col(cls(names{m})), 'LineWidth', 0.7 + 0.8*strcmp(names{m},'data'));
        end
    end
    set(ax, 'XScale','log', 'XTick', [3 10 30 100]); grid on; box on; ylim([0.4 0.85]); xlim([2 200]);
    xlabel('inter-event separation (m)'); ylabel(sprintf('median coda coherence, sensor %c', 'A'+g-1));
    panel_label(ax, g, 'tr');
    if g == 1
        h = gobjects(3,1); for k = 1:3, h(k) = plot(NaN, NaN, ls{k}, 'Color', 'k', 'LineWidth', 1); end
        lg = legend(h, {'\sigma_{loc} = 5 m', '\sigma_{loc} = 20 m', '\sigma_{loc} = 50 m'}, 'Location', 'southwest'); lg.Box = 'off';
    end
end
pub_export(f, fullfile(here,'..','figs','fig_d2_perturb.png'), 17.4, 7.5, 7); close(f);
