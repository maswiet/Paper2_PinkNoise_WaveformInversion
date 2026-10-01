% S25  Uncertainty of the inter-event coda coherence D2 (reviewer points 10, 17, 26):
% (a) event-level (cluster) bootstrap: events resampled with replacement, each pair weighted by
%     the product of its events' multiplicities (pairs sharing an event are not independent);
% (b) hypocentre perturbation: catalogue positions perturbed by isotropic N(0, s), s = 5, 10 m,
%     before binning by separation (100 draws).
% Output: results/table_d2_uncertainty.csv, figs/fig_d2_uncertainty.png
maxNumCompThreads(2);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
G = forge_setup();
names = {'data','H','A','L','L6','GRAD3','FZ','P3','P6','E15','G50'};
names = names(cellfun(@(n) exist(fullfile(out, sprintf('shape_%s_forge.mat', n)), 'file') == 2, names));
edges = [0 5 10 20 40 80 160]; nb = numel(edges)-1; B = 500; rng(21);
pos = containers.Map('KeyType','double','ValueType','any');
for i = 1:numel(G.ev), pos(G.ev(i)) = G.evxyz(i,:); end
rows = {}; CUR = struct();
for m = 1:numel(names)
    Z = load(fullfile(out, sprintf('shape_%s_forge.mat', names{m})));
    for g = 1:2
        d = Z.pairs{g}(:,1); c = Z.pairs{g}(:,2); pid = Z.pairidx{g};
        [ue, ~, ic] = unique(pid(:)); p1 = ic(1:end/2); p2 = ic(end/2+1:end); nu = numel(ue);
        med = binmed(d, c, ones(size(c)), edges);
        bs = NaN(B, nb);
        for b = 1:B
            cnt = accumarray(randi(nu, nu, 1), 1, [nu 1]);
            bs(b,:) = binmed(d, c, cnt(p1).*cnt(p2), edges);
        end
        xyz = cell2mat(arrayfun(@(e) pos(e), ue, 'uni', 0));
        pert = NaN(2, nb);
        for si = 1:2
            s = 5*si; pp = NaN(100, nb);
            for r = 1:100
                X = xyz + s*randn(size(xyz)); dd = sqrt(sum((X(p1,:) - X(p2,:)).^2, 2));
                pp(r,:) = binmed(dd, c, ones(size(c)), edges);
            end
            pert(si,:) = mean(pp, 1, 'omitnan');
        end
        np = histcounts(d, edges);
        ne_bin = arrayfun(@(k) numel(unique(pid(d > edges(k) & d <= edges(k+1), :))), 1:nb);
        for k = 1:nb
            rows(end+1,:) = {names{m}, char('A'+g-1), edges(k), edges(k+1), med(k), prctile(bs(:,k),2.5), ...
                prctile(bs(:,k),97.5), pert(1,k), pert(2,k), np(k), ne_bin(k)}; %#ok<SAGROW>
        end
        CUR(m,g).med = med; CUR(m,g).lo = prctile(bs,2.5); CUR(m,g).hi = prctile(bs,97.5); CUR(m,g).pert = pert;
    end
end
T = cell2table(rows, 'VariableNames', {'model','sensor','d_lo_m','d_hi_m','median_coh','lo95','hi95', ...
    'median_coh_loc5m','median_coh_loc10m','n_pairs','n_events'});
writetable(T, fullfile(out,'table_d2_uncertainty.csv'));
disp(T(strcmp(T.model,'data') | strcmp(T.model,'H'), :));
% slope of coherence vs log10 separation (5-160 m) with bootstrap CI, per model and sensor
f1 = figure('Position',[40 40 1150 440], 'Visible','off');
cc = (edges(1:end-1) + edges(2:end))/2; cc(1) = 3;
col = containers.Map({'data','H','A','L','L6','GRAD3','FZ','P3','P6','E15','G50'}, {[0 0 0],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.55 .25 .65],[.5 .5 .5],[.55 .25 .65],[.85 .33 .1],[.85 .33 .1],[.2 .45 .75],[.3 .6 .3]});
for g = 1:2
    subplot(1,2,g); hold on
    for m = 1:numel(names)
        C = CUR(m,g); lw = 1.2 + 1.3*strcmp(names{m},'data');
        errorbar(cc*(1 + 0.025*(m-6)), C.med, C.med - C.lo, C.hi - C.med, '-o', 'Color', col(names{m}), 'LineWidth', lw, 'MarkerSize', 3);
        if strcmp(names{m},'data')
            plot(cc, C.pert(2,:), 'k--', 'LineWidth', 1);
        end
    end
    set(gca,'XScale','log'); grid on; box on; ylim([0.3 1]); xlim([2 200]);
    xlabel('inter-event separation (m)'); ylabel('median coda coherence'); panel_label(gca, g, 'tr');
end
exportgraphics(f1, fullfile(here,'..','figs','fig_d2_uncertainty.png'), 'Resolution', 140); close(f1);

function med = binmed(d, c, w, edges)
med = NaN(1, numel(edges)-1);
for k = 1:numel(edges)-1
    m = d > edges(k) & d <= edges(k+1) & w > 0;
    if nnz(m) < 20, continue; end
    [cs, o] = sort(c(m)); ww = w(m); ww = cumsum(ww(o)); med(k) = cs(find(ww >= ww(end)/2, 1));
end
end
