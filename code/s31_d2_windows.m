% S31  Inter-event coda coherence D2 in three windows (second review): S+15..100 ms (original),
% S+40..100 ms (same as the late-energy test) and S+50..110 ms.
% (a) cluster bootstrap over events (B = 500) of the binned median coherence;
% (b) location stress test: the binning coordinates of DATA AND MODELS are perturbed by isotropic
%     N(0, s) with s = 5, 10, 20, 30, 50 m (50 draws), so that models are degraded like the data;
% (c) separation-independent statistic: median coherence of all pairs within 160 m, with cluster
%     bootstrap - insensitive to how location error moves pairs between bins.
% Output: results/table_d2_windows.csv, table_d2_perturb.csv, table_d2_allpairs.csv,
%         figs/fig_d2_windows.png (main text), figs/fig_d2_perturb.png (Additional file 1)
maxNumCompThreads(4);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
G = forge_setup();
names = {'data','H','L','GRAD3','A','L6','FZ','P3','P6','E15','G50'};
cls = containers.Map(names, {'data','smooth','smooth','smooth','smooth','det','det','pink','pink','exp','gau'});
WIN = {'', 'S+15-100 ms'; '_w40', 'S+40-100 ms'; '_w50', 'S+50-110 ms'};
edges = [0 5 10 20 40 80 160]; nb = numel(edges)-1; B = 500; SP = [5 10 20 30 50]; ND = 50;
pos = containers.Map('KeyType','double','ValueType','any');
for i = 1:numel(G.ev), pos(G.ev(i)) = G.evxyz(i,:); end
rows = {}; rowsP = {}; rowsA = {}; CUR = struct();
for w = 1:size(WIN,1)
    rng(21);
    for m = 1:numel(names)
        f = fullfile(out, sprintf('shape_%s_forge%s.mat', names{m}, WIN{w,1}));
        if ~exist(f, 'file'), fprintf('missing %s\n', f); continue; end
        Z = load(f);
        for g = 1:2
            d = Z.pairs{g}(:,1); c = Z.pairs{g}(:,2); pid = Z.pairidx{g};
            [ue, ~, ic] = unique(pid(:)); p1 = ic(1:end/2); p2 = ic(end/2+1:end); nu = numel(ue);
            med = binmed(d, c, ones(size(c)), edges);
            bs = NaN(B, nb); ba = NaN(B, 1); inr = d <= edges(end);
            for b = 1:B
                cnt = accumarray(randi(nu, nu, 1), 1, [nu 1]); wgt = cnt(p1).*cnt(p2);
                bs(b,:) = binmed(d, c, wgt, edges);
                ba(b) = wmed(c(inr), wgt(inr));
            end
            np = histcounts(d, edges);
            ne_bin = arrayfun(@(k) numel(unique(pid(d > edges(k) & d <= edges(k+1), :))), 1:nb);
            for k = 1:nb
                rows(end+1,:) = {names{m}, cls(names{m}), WIN{w,2}, char('A'+g-1), edges(k), edges(k+1), med(k), ...
                    prctile(bs(:,k),2.5), prctile(bs(:,k),97.5), np(k), ne_bin(k)}; %#ok<SAGROW>
            end
            rowsA(end+1,:) = {names{m}, cls(names{m}), WIN{w,2}, char('A'+g-1), median(c(inr)), prctile(ba,2.5), prctile(ba,97.5), nnz(inr)}; %#ok<SAGROW>
            xyz = cell2mat(arrayfun(@(e) pos(e), ue, 'uni', 0));
            for si = 1:numel(SP)
                pp = NaN(ND, nb);
                for r = 1:ND
                    X = xyz + SP(si)*randn(size(xyz)); dd = sqrt(sum((X(p1,:) - X(p2,:)).^2, 2));
                    pp(r,:) = binmed(dd, c, ones(size(c)), edges);
                end
                pm = mean(pp, 1, 'omitnan');
                for k = 1:nb
                    rowsP(end+1,:) = {names{m}, cls(names{m}), WIN{w,2}, char('A'+g-1), SP(si), edges(k), edges(k+1), pm(k)}; %#ok<SAGROW>
                end
            end
            CUR(w,m,g).med = med; CUR(w,m,g).lo = prctile(bs,2.5); CUR(w,m,g).hi = prctile(bs,97.5); %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'model','class','window','sensor','d_lo_m','d_hi_m','median_coh','lo95','hi95','n_pairs','n_events'});
writetable(T, fullfile(out, 'table_d2_windows.csv'));
TP = cell2table(rowsP, 'VariableNames', {'model','class','window','sensor','perturb_m','d_lo_m','d_hi_m','median_coh'});
writetable(TP, fullfile(out, 'table_d2_perturb.csv'));
TA = cell2table(rowsA, 'VariableNames', {'model','class','window','sensor','median_coh_all','lo95','hi95','n_pairs'});
writetable(TA, fullfile(out, 'table_d2_allpairs.csv')); disp(TA);
% ---------------- noise floor: FORGE records in a noise-only window (S+550..610 ms) ----------------
NF = NaN(2, nb); NA = NaN(2, 3); fnz = fullfile(out, 'shape_data_forge_w550.mat');
if exist(fnz, 'file')
    Z = load(fnz);
    for g = 1:2
        d = Z.pairs{g}(:,1); c = Z.pairs{g}(:,2); pid = Z.pairidx{g}; NF(g,:) = binmed(d, c, ones(size(c)), edges);
        [ue, ~, ic] = unique(pid(:)); p1 = ic(1:end/2); p2 = ic(end/2+1:end); nu = numel(ue); inr = d <= edges(end);
        ba = NaN(B,1); rng(22);
        for b = 1:B, cnt = accumarray(randi(nu, nu, 1), 1, [nu 1]); w = cnt(p1).*cnt(p2); ba(b) = wmed(c(inr), w(inr)); end
        NA(g,:) = [median(c(inr)) prctile(ba,2.5) prctile(ba,97.5)];
        TA = [TA; {'data (noise window)', 'noise', 'S+550-610 ms', char('A'+g-1), NA(g,1), NA(g,2), NA(g,3), nnz(inr)}]; %#ok<AGROW>
    end
    writetable(TA, fullfile(out, 'table_d2_allpairs.csv'));
end
% ---------------- figure: 2 sensors x 3 windows ----------------
col = containers.Map({'data','smooth','det','pink','exp','gau'}, {[0 0 0],[.55 .55 .55],[.55 .25 .65],[.85 .33 .1],[.2 .45 .75],[.3 .6 .3]});
cc = (edges(1:end-1) + edges(2:end))/2; cc(1) = 3;
f1 = figure('Visible','off'); np_ = 0;
for g = 1:2
    for w = 1:size(WIN,1)
        np_ = np_ + 1; ax = subplot(2, 3, np_); hold on
        for m = numel(names):-1:1
            if size(CUR,1) < w || size(CUR,2) < m || isempty(CUR(w,m,g).med), continue; end
            C = CUR(w,m,g); isd = strcmp(names{m},'data'); lw = 0.7 + 0.9*isd;
            errorbar(cc*(1 + 0.022*(m-6)), C.med, C.med - C.lo, C.hi - C.med, '-o', 'Color', col(cls(names{m})), ...
                'LineWidth', lw, 'MarkerSize', 2 + 1.5*isd, 'MarkerFaceColor', col(cls(names{m})), 'CapSize', 2);
        end
        if any(isfinite(NF(g,:))), plot(cc, NF(g,:), 'k:', 'LineWidth', 1.1); end
        set(ax, 'XScale','log'); grid on; box on; ylim([0.35 1]); xlim([2 200]); set(ax, 'XTick', [3 10 30 100]);
        if g == 2, xlabel('inter-event separation (m)'); else, set(ax, 'XTickLabel', []); end
        if w == 1, ylabel(sprintf('median coda coherence, sensor %c', 'A'+g-1)); else, set(ax, 'YTickLabel', []); end
        text(0.04, 0.95, WIN{w,2}, 'Units','normalized', 'FontSize', 7, 'VerticalAlignment', 'top');
        panel_label(ax, np_, 'tr');
        set(ax, 'Position', [0.09 + 0.305*(w-1), 0.56 - 0.40*(g-1) + 0.02, 0.285, 0.38]);
    end
end
ks = {'data','smooth','det','pink','exp','gau'};
lab = {'FORGE','smooth (H, L, GRAD3, A)','6-m layering, fracture zones','pink','exponential','Gaussian'};
hh = gobjects(7,1);
for q = 1:6, hh(q) = plot(NaN, NaN, '-o', 'Color', col(ks{q}), 'MarkerFaceColor', col(ks{q}), 'MarkerSize', 3, 'LineWidth', 0.7 + 0.9*(q==1)); end
hh(7) = plot(NaN, NaN, 'k:', 'LineWidth', 1.1);
lg = legend(hh, [lab {'FORGE noise window'}], 'Orientation', 'horizontal'); lg.NumColumns = 4; lg.Box = 'off';
lg.ItemTokenSize = [14 8]; lg.Position = [0.08 0.005 0.90 0.08];
pub_export(f1, fullfile(here,'..','figs','fig_d2_windows.png'), 17.4, 13, 7); close(f1);
% ---------------- figure: perturbation stress test (drawn from the CSV) ----------------
fig_d2_perturb;

function med = binmed(d, c, w, edges)
med = NaN(1, numel(edges)-1);
for k = 1:numel(edges)-1
    m = d > edges(k) & d <= edges(k+1) & w > 0;
    if nnz(m) < 20, continue; end
    med(k) = wmed(c(m), w(m));
end
end

function v = wmed(c, w)
m = w > 0; c = c(m); w = w(m);
if isempty(c), v = NaN; return; end
[cs, o] = sort(c); ww = cumsum(w(o)); v = cs(find(ww >= ww(end)/2, 1));
end
