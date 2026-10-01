function fig_rev_panels()
% FIG_REV_PANELS  Figures for the revised manuscript:
%  fig_coda_ratio_rev.png  per-event log10(predicted/observed coda ratio), all smooth baselines and
%                          the stochastic families (non-overlapping windows, 416 events)
%  fig_sigma_eff.png       median log10(pred/obs) versus sigma per family and sensor (S29)
%  fig_phi_family.png      Phi with event-bootstrap 95 % intervals and holdout scores (S20)
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results'); fd = fullfile(here,'..','figs');
x3 = log10(3);
% ---------------- per-event coda ratio ----------------
grp = {'H','smooth'; 'L','smooth'; 'GRAD3','smooth'; 'A','smooth'; 'L6','det'; 'FZ','det'; ...
       'P1','pink'; 'P3','pink'; 'P6','pink'; 'VSD','pinkvs'; 'VSO','pinkvs'; ...
       'X_exp_0.09_15_11','exp'; 'E15','exp'; 'X_exp_0.09_50_11','exp'; 'E50','exp'; ...
       'X_gau_0.09_15_11','gau'; 'G15','gau'; 'X_gau_0.09_50_11','gau'; 'G50','gau'};
lab = containers.Map( ...
    {'H','L','L6','GRAD3','FZ','A','P1','P3','P6','VSD','VSO','X_exp_0.09_15_11','E15','X_exp_0.09_50_11','E50', ...
     'X_gau_0.09_15_11','G15','X_gau_0.09_50_11','G50'}, ...
    {'H','L (50 m)','L6 (6 m)','GRAD3','FZ','A (VTI)','pink .045','pink .09','pink .13','pink .13 Vs-ind.','pink Vs-dom.', ...
     'exp a15 .09','exp a15 .13','exp a50 .09','exp a50 .13','Gau a15 .09','Gau a15 .13','Gau a50 .09','Gau a50 .13'});
col = containers.Map({'smooth','det','pink','pinkvs','exp','gau'}, {[0.45 0.45 0.45],[0.55 0.25 0.65],[0.85 0.45 0.1],[0.95 0.65 0.3],[0.2 0.45 0.8],[0.2 0.65 0.35]});
have = cellfun(@(n) exist(fullfile(out, ['s19_' n '.mat']), 'file') == 2, grp(:,1)); grp = grp(have,:);
f1 = figure('Position',[40 40 1500 520], 'Visible','off');
for g = 1:2
    subplot(1,2,g); hold on
    for m = 1:size(grp,1)
        S = load(fullfile(out, ['s19_' grp{m,1} '.mat']), 'lr'); v = S.lr(:,g,1,1); v = v(isfinite(v));
        q = prctile(v, [10 25 50 75 90]); c = col(grp{m,2});
        plot([m m], q([1 5]), '-', 'Color', c);
        patch(m + [-0.3 0.3 0.3 -0.3], q([2 2 4 4]), c, 'FaceAlpha', 0.35, 'EdgeColor', c);
        plot(m, q(3), 'd', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    end
    yline(0, 'k-'); yline([-x3 x3], 'k:');
    set(gca, 'XTick', 1:size(grp,1), 'XTickLabel', cellfun(@(n) lab(n), grp(:,1), 'UniformOutput', false), 'XTickLabelRotation', 50);
    xlim([0.4 size(grp,1)+0.6]); ylim([-3.2 2.2]); grid on; box on
    ylabel('log_{10}(predicted / observed coda ratio)');
    panel_label(gca, g);
end
exportgraphics(f1, fullfile(fd, 'fig_coda_ratio_rev.png'), 'Resolution', 140); close(f1);
% ---------------- sigma_eff curves ----------------
fn = fullfile(out, 'table_sigma_curves.csv');
if exist(fn, 'file')
    C = readtable(fn); fams = unique(C.family, 'stable');
    fc = containers.Map({'pink15','pink18','exp15','exp50','gau15','gau50'}, ...
        {[0.85 0.45 0.1],[0.95 0.7 0.35],[0.2 0.45 0.8],[0.55 0.7 0.95],[0.2 0.65 0.35],[0.6 0.85 0.6]});
    fl = containers.Map({'pink15','pink18','exp15','exp50','gau15','gau50'}, ...
        {'pink p = 1.5','pink p = 1.8','exponential a = 15 m','exponential a = 50 m','Gaussian a = 15 m','Gaussian a = 50 m'});
    sm = {'H','L','GRAD3','A'}; sv = NaN(numel(sm),2);
    for k = 1:numel(sm)
        f = fullfile(out, ['s19_' sm{k} '.mat']);
        if exist(f,'file'), S = load(f,'lr'); sv(k,:) = [median(S.lr(:,1,1,1),'omitnan') median(S.lr(:,2,1,1),'omitnan')]; end
    end
    f2 = figure('Position',[40 40 1250 480], 'Visible','off');
    for g = 1:2
        subplot(1,2,g); hold on
        yl = [min(sv(:,g)) max(sv(:,g))];
        patch([0.015 0.2 0.2 0.015], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor','none', 'DisplayName', 'smooth media (range)');
        for k = 1:numel(fams)
            r = C(strcmp(C.family, fams{k}), :); y = r{:, 3+g}; ok = isfinite(y);
            plot(r.sigma(ok), y(ok), '-o', 'Color', fc(fams{k}), 'MarkerFaceColor', fc(fams{k}), 'LineWidth', 1.3, 'DisplayName', fl(fams{k}));
        end
        dl = {'L6','6-m layering (L6)','-.'; 'FZ','fracture zones (FZ)','--'};
        for q = 1:2
            f = fullfile(out, ['s19_' dl{q,1} '.mat']);
            if exist(f,'file'), S = load(f,'lr'); yline(median(S.lr(:,g,1,1),'omitnan'), dl{q,3}, 'Color', [0.55 0.25 0.65], 'LineWidth', 1.3, 'DisplayName', dl{q,2}); end
        end
        yline(0, 'k-', 'HandleVisibility','off'); yline([-x3 x3], 'k:', 'HandleVisibility','off');
        set(gca, 'XScale', 'log'); xlim([0.015 0.2]); set(gca, 'XTick', [0.02 0.045 0.09 0.13 0.18]);
        grid on; box on; xlabel('\sigma (standard deviation of ln V)'); ylabel('median log_{10}(predicted / observed coda ratio)');
        if g == 1, legend('Location', 'northwest'); end
        panel_label(gca, g);
    end
    exportgraphics(f2, fullfile(fd, 'fig_sigma_eff.png'), 'Resolution', 140); close(f2);
end
% ---------------- Phi bootstrap and holdout ----------------
fn = fullfile(out, 'table_phi_bootstrap.csv'); fh = fullfile(out, 'table_holdout.csv');
if exist(fn, 'file') && exist(fh, 'file')
    T = readtable(fn); T = sortrows(T, 'PHI');
    fcol = containers.Map({'smooth','deterministic','pink','pink-VpVs','exp15','exp50','gau15','gau50'}, ...
        {[0.45 0.45 0.45],[0.55 0.25 0.65],[0.85 0.45 0.1],[0.95 0.65 0.3],[0.2 0.45 0.8],[0.55 0.7 0.95],[0.2 0.65 0.35],[0.6 0.85 0.6]});
    f3 = figure('Position',[40 40 1500 640], 'Visible','off');
    subplot(1,3,[1 2]); hold on
    for m = 1:height(T)
        c = fcol(T.family{m});
        plot([m m], [T.PHI_lo95(m) T.PHI_hi95(m)], '-', 'Color', c, 'LineWidth', 1.5);
        plot(m, T.PHI(m), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
    end
    nm = regexprep(T.model, '^X_exp_', 'E '); nm = regexprep(nm, '^X_gau_', 'G '); nm = regexprep(nm, '^X_pink_', 'P ');
    nm = strrep(nm, '_', '/');
    set(gca, 'XTick', 1:height(T), 'XTickLabel', nm, 'XTickLabelRotation', 70, 'FontSize', 7);
    xlim([0 height(T)+1]); grid on; box on; ylabel('\Phi (event-bootstrap 95 % interval)');
    set(gca, 'Position', [0.05 0.24 0.62 0.72]);
    panel_label(gca, 1);
    subplot(1,3,3); H = readtable(fh); fams = {'smooth','deterministic','pink','exp15','exp50','gau15','gau50','pink-VpVs'}; hold on
    for f = 1:numel(fams)
        v = H.PHI_test(strcmp(H.family, fams{f}));
        if isempty(v), continue; end
        q = prctile(v, [0 25 50 75 100]); c = fcol(fams{f});
        plot([f f], q([1 5]), '-', 'Color', c);
        patch(f + [-0.3 0.3 0.3 -0.3], q([2 2 4 4]), c, 'FaceAlpha', 0.35, 'EdgeColor', c);
        plot(f, q(3), 'd', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    end
    set(gca, 'XTick', 1:numel(fams), 'XTickLabel', fams, 'XTickLabelRotation', 40); xlim([0.4 numel(fams)+0.6]);
    grid on; box on; ylabel('\Phi on held-out events (best training member)');
    set(gca, 'Position', [0.74 0.24 0.24 0.72]);
    panel_label(gca, 2);
    exportgraphics(f3, fullfile(fd, 'fig_phi_family.png'), 'Resolution', 140); close(f3);
end
end
