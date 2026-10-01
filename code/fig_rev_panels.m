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
f1 = figure('Visible','off');
for g = 1:2
    subplot(1,2,g); hold on
    for m = 1:size(grp,1)
        S = load(fullfile(out, ['s19_' grp{m,1} '.mat']), 'lr'); v = S.lr(:,g,1,1); v = v(isfinite(v));
        q = prctile(v, [10 25 50 75 90]); c = col(grp{m,2});
        plot([m m], q([1 5]), '-', 'Color', c);
        patch(m + [-0.3 0.3 0.3 -0.3], q([2 2 4 4]), c, 'FaceAlpha', 0.35, 'EdgeColor', c);
        plot(m, q(3), 'd', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
    end
    yline(0, 'k-'); yline([-x3 x3], 'k:');
    set(gca, 'XTick', 1:size(grp,1), 'XTickLabel', cellfun(@(n) lab(n), grp(:,1), 'UniformOutput', false), 'XTickLabelRotation', 50);
    xlim([0.4 size(grp,1)+0.6]); ylim([-3.2 2.2]); grid on; box on
    if g == 1, ylabel('log_{10}(predicted / observed late energy)'); else, set(gca, 'YTickLabel', []); end
    set(gca, 'Position', [0.07 + 0.475*(g-1) 0.25 0.445 0.72]);
    panel_label(gca, g);
end
pub_export(f1, fullfile(fd, 'fig_coda_ratio_rev.png'), 17.4, 8.5, 7); close(f1);
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
    f2 = figure('Visible','off');
    for g = 1:2
        subplot(1,2,g); hold on
        yl = [min(sv(:,g)) max(sv(:,g))];
        patch([0.015 0.2 0.2 0.015], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor','none', 'DisplayName', 'smooth media (range)');
        for k = 1:numel(fams)
            r = C(strcmp(C.family, fams{k}), :); y = r{:, 3+g}; ok = isfinite(y);
            plot(r.sigma(ok), y(ok), '-o', 'Color', fc(fams{k}), 'MarkerFaceColor', fc(fams{k}), 'LineWidth', 1, 'MarkerSize', 3.5, 'DisplayName', fl(fams{k}));
        end
        dl = {'L6','6-m layering (L6)','-.'; 'FZ','fracture zones (FZ)','--'};
        for q = 1:2
            f = fullfile(out, ['s19_' dl{q,1} '.mat']);
            if exist(f,'file'), S = load(f,'lr'); yline(median(S.lr(:,g,1,1),'omitnan'), dl{q,3}, 'Color', [0.55 0.25 0.65], 'LineWidth', 1, 'DisplayName', dl{q,2}); end
        end
        yline(0, 'k-', 'HandleVisibility','off'); yline([-x3 x3], 'k:', 'HandleVisibility','off');
        set(gca, 'XScale', 'log'); xlim([0.015 0.2]); set(gca, 'XTick', [0.02 0.045 0.09 0.13 0.18]);
        grid on; box on; xlabel('\sigma (standard deviation of ln V)'); ylabel('median log_{10}(predicted / observed)');
        if g == 1, lg = legend('Orientation', 'horizontal'); lg.NumColumns = 3; lg.Box = 'off'; lg.ItemTokenSize = [14 8]; end
        set(gca, 'Position', [0.08 + 0.5*(g-1) 0.31 0.40 0.66]);
        panel_label(gca, g);
    end
    lg.Position = [0.06 0.0 0.90 0.14];
    pub_export(f2, fullfile(fd, 'fig_sigma_eff.png'), 17.4, 9); close(f2);
end
% ---------------- Phi bootstrap and holdout ----------------
fn = fullfile(out, 'table_phi_bootstrap.csv'); fh = fullfile(out, 'table_holdout.csv');
if exist(fn, 'file') && exist(fh, 'file')
    T = readtable(fn); T = sortrows(T, 'PHI');
    fcol = containers.Map({'smooth','deterministic','pink','pink-VpVs','exp15','exp50','gau15','gau50'}, ...
        {[0.45 0.45 0.45],[0.55 0.25 0.65],[0.85 0.45 0.1],[0.95 0.65 0.3],[0.2 0.45 0.8],[0.55 0.7 0.95],[0.2 0.65 0.35],[0.6 0.85 0.6]});
    f3 = figure('Visible','off');
    ax = subplot(1,2,1); hold on; n = height(T);
    for m = 1:n
        c = fcol(T.family{m}); y = n + 1 - m;
        plot([T.PHI_lo95(m) T.PHI_hi95(m)], [y y], '-', 'Color', c, 'LineWidth', 1.2);
        plot(T.PHI(m), y, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', 3.5, 'LineWidth', 0.4);
    end
    nm = cellfun(@(m, fam, sg, sh) nicename(m, fam, sg, sh), T.model, T.family, num2cell(T.sigma), num2cell(T.shape), 'UniformOutput', false);
    set(ax, 'YTick', 1:n, 'YTickLabel', flipud(nm), 'YLim', [0.3 n+0.7], 'TickLabelInterpreter', 'none');
    grid on; box on; xlabel('\Phi (event-bootstrap 95 % interval)');
    set(ax, 'Position', [0.17 0.08 0.36 0.90]);
    panel_label(ax, 1, 'tr');
    ax2 = subplot(1,2,2); H = readtable(fh); hold on
    fams = {'smooth','deterministic','pink-VpVs','gau15','gau50','exp15','exp50','pink'};
    flab = {'smooth','deterministic','pink, var. Vp/Vs','Gaussian a 15 m','Gaussian a 50 m','exponential a 15 m','exponential a 50 m','pink'};
    for f = 1:numel(fams)
        v = H.PHI_test(strcmp(H.family, fams{f}));
        if isempty(v), continue; end
        q = prctile(v, [0 25 50 75 100]); c = fcol(fams{f}); y = numel(fams) + 1 - f;
        plot(q([1 5]), [y y], '-', 'Color', c);
        patch(q([2 4 4 2]), y + [-0.3 -0.3 0.3 0.3], c, 'FaceAlpha', 0.35, 'EdgeColor', c);
        plot(q(3), y, 'd', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
    end
    set(ax2, 'YTick', 1:numel(fams), 'YTickLabel', fliplr(flab), 'YLim', [0.4 numel(fams)+0.6]);
    grid on; box on; xlabel('\Phi on held-out events');
    set(ax2, 'Position', [0.74 0.08 0.24 0.90]);
    panel_label(ax2, 2, 'br');
    pub_export(f3, fullfile(fd, 'fig_phi_family.png'), 17.4, 15, 7); close(f3);
end
end

function s = nicename(m, fam, sg, sh)
% readable label: family, sigma, shape parameter and seed (seed 11 unless stated)
fixed = containers.Map({'H','L','A','GRAD3','L6','FZ','VSD','VSO'}, {'homogeneous','layered 50 m','VTI', ...
    '3-D gradients','layered 6 m','fracture zones','pink .13, indep. Vs','pink, Vs-dominated'});
if isKey(fixed, m), s = fixed(m); return; end
if startsWith(m, 'X_'), t = strsplit(m, '_'); seed = t{end};
elseif endsWith(m, 'b') || startsWith(m, 'R'), seed = '33';
else, seed = '11';
end
switch fam
    case 'pink',  s = sprintf('pink %.2f, p %.1f, s%s', sg, sh, seed);
    case {'exp15','exp50'}, s = sprintf('exp %.2f, a %d, s%s', sg, round(sh), seed);
    otherwise,    s = sprintf('Gauss %.2f, a %d, s%s', sg, round(sh), seed);
end
s = regexprep(s, ' 0\.', ' .');
end
