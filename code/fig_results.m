% Final result figures: (a) misfit ranking for synthetic tests and FORGE data;
% (b) median S-coda envelopes data vs model classes; (c) sigma-p misfit map.
maxNumCompThreads(2);
cmpf = @(o) load(sprintf('../results/compare_%s_forge.mat', o));
obs = {'T1','T2','data'}; ttl = {'Synthetic truth T1: pink \sigma=0.045, p=1.5', ...
       'Synthetic truth T2: layered + VTI', 'Utah FORGE 2022 (56-32 sensors)'};
cls = @(n) 1 + ismember(n(1),'PTR') + 2*ismember(n(1),'EG');   % 1 deterministic, 2 pink, 3 competitor
col = [0.35 0.35 0.35; 0.85 0.33 0.10; 0.2 0.45 0.75];
figure('Position',[50 50 1500 620]);
for i = 1:3
    C = cmpf(obs{i}); R = C.RES; [~,o] = sort([R.PHI]); R = R(o);
    subplot(1,3,i); hold on
    for m = 1:numel(R)
        barh(m, R(m).PHI, 'FaceColor', col(cls(R(m).name),:));
        text(R(m).PHI+0.01, m, sprintf('%.2f', R(m).PHI), 'FontSize',8);
    end
    set(gca,'YTick',1:numel(R),'YTickLabel',{R.name},'YDir','reverse'); xlabel('\Phi = KS + ENV misfit');
    box on; xlim([0 max([R.PHI])*1.2]); panel_label(gca, i, 'br');
end
exportgraphics(gcf,'../figs/fig_misfit_ranking.png','Resolution',150);

% (b) envelopes
Df = load('../data/forge_features.mat');
mods = {'H','A','E50','G50','P3','P9'}; lag = (-0.08:1e-3:0.12)';
figure('Position',[50 50 1300 420]); pan = [1 1; 1 2; 2 2]; k = 0;
for q = 1:3
    g = pan(q,1); ib = pan(q,2);
    subplot(1,3,q); E = [Df.D(g,ib).F.envS];
    semilogy(lag, median(E,2,'omitnan'), 'k', 'LineWidth', 2.5); hold on
    for m = 1:numel(mods)
        f = sprintf('../results/feat_%s_forge.mat', mods{m}); if ~exist(f,'file'), continue; end
        L = load(f); [~,im] = intersect(L.M(g,ib).k, Df.D(g,ib).k);
        semilogy(lag, median([L.M(g,ib).F(im).envS],2,'omitnan'), 'LineWidth', 1.3);
    end
    legend(['FORGE data' mods], 'Location','northeast'); grid on; ylim([1e-3 1.5]); xlim([-0.06 0.12]);
    xlabel('time after S-envelope peak (s)'); ylabel('median E / E_{peak}');
    panel_label(gca, q, 'bl');
end
exportgraphics(gcf,'../figs/fig_envelopes.png','Resolution',150);

% (c) sigma-p misfit map (data); repeated (sigma,p) = realisations -> mean, range
C = cmpf('data'); R = C.RES;
sp = struct('P1',[0.045 1.5],'P2',[0.045 1.8],'P3',[0.09 1.5],'P4',[0.02 1.5],'P5',[0.045 1.2], ...
            'P6',[0.13 1.5],'P7',[0.09 1.2],'P8',[0.09 1.8],'P9',[0.18 1.5],'T1',[0.045 1.5], ...
            'R1',[0.045 1.5],'R2',[0.045 1.5],'R3',[0.09 1.5], ...
            'P10',[0.25 1.5],'P11',[0.18 1.8],'P12',[0.22 1.5],'R5',[0.18 1.5]);
isp = arrayfun(@(r) isfield(sp, r.name), R);
X = cell2mat(arrayfun(@(r) sp.(r.name), R(isp), 'uni', 0)'); phi = [R(isp).PHI]';
[U,~,j] = unique(X,'rows');
figure('Position',[50 50 600 460]); hold on
for u = 1:size(U,1)
    v = phi(j==u);
    scatter(U(u,1), U(u,2), 420, mean(v), 'filled', 'MarkerEdgeColor','k');
    if numel(v) > 1
        text(U(u,1), U(u,2)-0.13, {sprintf('%.2f', mean(v)), sprintf('(%.2f-%.2f, n=%d)', min(v), max(v), numel(v))}, ...
            'FontSize',8, 'HorizontalAlignment','center');
    else
        text(U(u,1), U(u,2)+0.09, sprintf('%.2f', v), 'FontSize',8, 'HorizontalAlignment','center');
    end
end
det = R(~isp & ~arrayfun(@(r) any(r.name(1)=='EG'), R)); [dm,k] = min([det.PHI]);
cmp = R(arrayfun(@(r) any(r.name(1)=='EG'), R)); [cm,kc] = min([cmp.PHI]);
cb = colorbar; cb.Label.String = '\Phi (FORGE)'; colormap(flipud(parula));
xlabel('\sigma (std of ln V at 5-m grid scale)'); ylabel('3D filter exponent p'); box on
fprintf('sigma-p map: best smooth %s = %.2f, best competitor %s = %.2f\n', det(k).name, dm, cmp(kc).name, cm);xlim([0 0.28]); ylim([1.0 2.0]);
exportgraphics(gcf,'../figs/fig_sigma_p_map.png','Resolution',150);
