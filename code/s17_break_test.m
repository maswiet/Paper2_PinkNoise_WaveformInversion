% S17  Qc^-1 break test over 25-50-100 Hz: slope1 = 25->50 Hz (sensor A, 5-m grid),
% slope2 = 50->100 Hz (sensor B, 2.5-m grid); curvature = slope2 - slope1.
% A correlation length a gives a change of slope near k a ~ 1; a scale-free medium none.
% Also the coda level at 100 Hz (sensor B). Output: results/table_break_test.csv, figs/fig_break_test.png
maxNumCompThreads(4);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
mods = {'H','P6','R6','E15','E15b','G50','G50b','EX5','EX9','EX17','EX34','EX67','EX134'};
av   = [NaN NaN NaN 15 15 50 50 5 9 17 34 67 134];
fam  = {'smooth','pink','pink','exponential','exponential','Gaussian','Gaussian', ...
        'exponential','exponential','exponential','exponential','exponential','exponential'};
for i = 1:numel(mods)                        % make sure the per-model measures exist
    for gr = {'5m','hf'}
        f = fullfile(out, sprintf('codamo_%s_%s.mat', mods{i}, gr{1}));
        if strcmp(mods{i},'EX5') && strcmp(gr{1},'5m'), continue; end   % a = 5 m not resolved on 5-m grid
        if ~exist(f,'file'), s15_model_coda_multioctave(mods{i}, gr{1}); end
    end
end
rows = {};
for i = 1:numel(mods)
    q25 = NaN; q50a = NaN;
    f5 = fullfile(out, sprintf('codamo_%s_5m.mat', mods{i}));
    if exist(f5,'file'), Z = load(f5); m = median(Z.QI{1},'omitnan'); q25 = m(1); q50a = m(2); end
    Z = load(fullfile(out, sprintf('codamo_%s_hf.mat', mods{i})));
    mq = median(Z.QI{2},'omitnan'); mc = median(Z.CL{2},'omitnan');
    s1 = log2(q50a/q25); s2 = log2(mq(2)/mq(1));
    rows(end+1,:) = {mods{i}, fam{i}, av(i), q25, q50a, mq(1), mq(2), s1, s2, s2 - s1, mc(1), mc(2)}; %#ok<SAGROW>
end
% FORGE data with bootstrap over events
D = load(fullfile(out,'coda_multioctave.mat'));
qa = D.res(1).QI(:,1:2); qa(~D.res(1).OK(:,1:2)) = NaN;
qb = D.res(2).QI(:,2:3); qb(~D.res(2).OK(:,2:3)) = NaN;
cb = D.res(2).CL(:,2:3); cb(~D.res(2).OK(:,2:3)) = NaN;
sl = @(q) log2(median(q(:,2),'omitnan')/median(q(:,1),'omitnan'));
rng(9); nb = 1000; B = zeros(nb,3);
for b = 1:nb
    ia = randi(size(qa,1), size(qa,1), 1); ib = randi(size(qb,1), size(qb,1), 1);
    B(b,:) = [sl(qa(ia,:)) sl(qb(ib,:)) sl(qb(ib,:)) - sl(qa(ia,:))];
end
d1 = sl(qa); d2 = sl(qb);
rows(end+1,:) = {'data','FORGE',NaN, median(qa(:,1),'omitnan'), median(qa(:,2),'omitnan'), ...
    median(qb(:,1),'omitnan'), median(qb(:,2),'omitnan'), d1, d2, d2 - d1, median(cb(:,1),'omitnan'), median(cb(:,2),'omitnan')};
T = cell2table(rows, 'VariableNames', {'model','family','a_m','Qc25_A','Qc50_A','Qc50_B','Qc100_B', ...
    'slope_25_50','slope_50_100','curvature','CL50_B','CL100_B'});
ci = prctile(B, [2.5 97.5]);
T.lo95 = NaN(height(T),3); T.hi95 = T.lo95; T.lo95(end,:) = ci(1,:); T.hi95(end,:) = ci(2,:);
T2 = splitvars(T, {'lo95','hi95'}, 'NewVariableNames', {{'s1_lo','s2_lo','curv_lo'},{'s1_hi','s2_hi','curv_hi'}});
disp(T2(:,[1 3 8 9 10 12]));
fprintf('data CI: slope1 [%.2f %.2f], slope2 [%.2f %.2f], curvature [%.2f %.2f]\n', ci(:));
writetable(T2, fullfile(out,'table_break_test.csv'));
% ---------------- figure ----------------
f1 = figure('Position',[40 40 1500 460], 'Visible','off');
ex = startsWith(T2.model,'EX'); aa = T2.a_m(ex);
subplot(1,3,1); hold on
for i = 1:height(T2)-1
    y = [T2.Qc25_A(i) T2.Qc50_A(i) T2.Qc50_B(i)*T2.Qc50_A(i)/T2.Qc50_B(i) T2.Qc100_B(i)*T2.Qc50_A(i)/T2.Qc50_B(i)];
    y = y / (isfinite(y(2))*y(2) + ~isfinite(y(2))*y(3));   % normalise at 50 Hz, stitch A (25-50) and B (50-100)
    x = [25 50 50 100];
    if startsWith(T2.model{i},'EX'), c = [0.2 0.45 0.75]*(0.4 + 0.6*log(T2.a_m(i))/log(134));
    elseif strcmp(T2.family{i},'pink'), c = [0.85 0.33 0.10];
    elseif strcmp(T2.family{i},'Gaussian'), c = [0.3 0.6 0.3];
    elseif strcmp(T2.family{i},'smooth'), c = [0.4 0.4 0.4];
    else, c = [0.45 0.65 0.9]; end
    plot(x([1 2]), y([1 2]), '-o', 'Color', c, 'MarkerSize', 4, 'HandleVisibility','off');
    plot(x([3 4]), y([3 4]), '-s', 'Color', c, 'MarkerSize', 4, 'HandleVisibility','off');
end
yd = [T2.Qc25_A(end) T2.Qc50_A(end) 1 T2.Qc100_B(end)/T2.Qc50_B(end)*T2.Qc50_A(end)]/T2.Qc50_A(end);
yd(3) = 1;
plot([25 50], yd(1:2), 'k-o', 'LineWidth', 2.5, 'MarkerFaceColor','k'); plot([50 100], yd([3 4]), 'k-s', 'LineWidth', 2.5, 'MarkerFaceColor','k');
set(gca,'XScale','log','YScale','log'); xlim([20 120]); grid on; box on
xlabel('frequency (Hz)'); ylabel('Q_c^{-1} / Q_c^{-1}(50 Hz)'); panel_label(gca, 1, 'tr');
subplot(1,3,2); hold on
patch([3 200 200 3], [ci(1,3) ci(1,3) ci(2,3) ci(2,3)], [0 0 0], 'FaceAlpha', 0.08, 'EdgeColor','none');
yline(T2.curvature(end), 'k-', 'LineWidth', 2);
plot(aa, T2.curvature(ex), 'o-', 'Color', [0.2 0.45 0.75], 'MarkerFaceColor', [0.2 0.45 0.75], 'LineWidth', 1.5);
pk = strcmp(T2.family,'pink'); gs = strcmp(T2.family,'Gaussian'); e15 = ismember(T2.model, {'E15','E15b'});
plot(repmat(150,nnz(pk),1), T2.curvature(pk), 'o', 'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerEdgeColor','k');
plot(repmat(50,nnz(gs),1), T2.curvature(gs), 's', 'MarkerFaceColor', [0.3 0.6 0.3], 'MarkerEdgeColor','k');
plot(repmat(15,nnz(e15),1), T2.curvature(e15), 'd', 'MarkerFaceColor', [0.45 0.65 0.9], 'MarkerEdgeColor','k');
plot(4, T2.curvature(strcmp(T2.model,'H')), 'v', 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor','k');
set(gca,'XScale','log'); xlim([3 200]); grid on; box on
xlabel('correlation length a (m)  [pink at right, smooth at left]'); ylabel('slope(50-100 Hz) - slope(25-50 Hz)');
panel_label(gca, 2, 'tr');
subplot(1,3,3); hold on
yline(T2.CL100_B(end), 'k-', 'LineWidth', 2);
plot(aa, T2.CL100_B(ex), 'o-', 'Color', [0.2 0.45 0.75], 'MarkerFaceColor', [0.2 0.45 0.75], 'LineWidth', 1.5);
plot(repmat(150,nnz(pk),1), T2.CL100_B(pk), 'o', 'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerEdgeColor','k');
plot(repmat(50,nnz(gs),1), T2.CL100_B(gs), 's', 'MarkerFaceColor', [0.3 0.6 0.3], 'MarkerEdgeColor','k');
plot(repmat(15,nnz(e15),1), T2.CL100_B(e15), 'd', 'MarkerFaceColor', [0.45 0.65 0.9], 'MarkerEdgeColor','k');
plot(4, T2.CL100_B(strcmp(T2.model,'H')), 'v', 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor','k');
set(gca,'XScale','log'); xlim([3 200]); grid on; box on
xlabel('correlation length a (m)'); ylabel('coda level at 100 Hz, sensor B (log_{10})'); panel_label(gca, 3, 'tr');
exportgraphics(f1, fullfile(here,'..','figs','fig_break_test.png'), 'Resolution', 150); close(f1);
