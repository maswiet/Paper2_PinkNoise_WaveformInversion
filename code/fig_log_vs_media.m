% 1-D vertical spectra through each random-medium class vs the 56-32 sonic log.
% Independent (non-seismic) test: a correlation-length medium has a spectral
% corner; the log and the pink-noise medium do not.
maxNumCompThreads(2);
G = forge_setup(); L = load('../data/log5632_stats.mat'); LOG = L.LOG;
n = G.n; dz = G.dx;
F = {pinknoise3d(n,1.5,11), pinknoise3d(n,1.8,11), randmedium3d(n,dz,'exp',15,11), ...
     randmedium3d(n,dz,'exp',50,11), randmedium3d(n,dz,'gauss',15,11), randmedium3d(n,dz,'gauss',50,11)};
nm = {'pink p=1.5','pink p=1.8','exponential a=15 m','exponential a=50 m','Gaussian a=15 m','Gaussian a=50 m'};
nz = n(3); k = (1:floor(nz/2)-1)'/(nz*dz);           % cycles/m
kr = [1/200 1/12];                                    % common band (12-200 m)
figure('Position',[50 50 760 560]);
% log: detrended lnV, average into 5-m bins then Welch-like segment average
lk = LOG.k; lp = LOG.P; in = lk >= kr(1) & lk <= kr(2);
lpb = movmean(lp, 15);
ref = @(kk, pp) pp / interp1(kk, pp, 1/50);
loglog(lk(in), ref(lk(in), lpb(in)), 'k', 'LineWidth', 2.5); hold on
cl = [0.85 0.33 0.10; 0.95 0.6 0.2; 0.2 0.45 0.75; 0.45 0.65 0.9; 0.3 0.6 0.3; 0.55 0.8 0.5];
slopes = zeros(1,numel(F)); rmsd = slopes;
c0 = polyfit(log(lk(in)), log(lp(in)), 1);
for i = 1:numel(F)
    f = F{i}; P = zeros(numel(k),1); c = 0;
    for ix = 5:6:n(1)-4, for iy = 5:6:n(2)-4
        l = squeeze(f(ix,iy,:)); l = l - mean(l);
        X = abs(fft(l)).^2; P = P + X(2:numel(k)+1); c = c + 1;
    end, end
    P = P/c; w = k >= kr(1) & k <= kr(2);
    cf = polyfit(log(k(w)), log(P(w)), 1); slopes(i) = -cf(1);
    loglog(k(w), ref(k(w), P(w)), 'Color', cl(i,:), 'LineWidth', 1.4);
    % misfit of normalised spectral shape vs log (log10 rms)
    li = interp1(log(lk(in)), log(ref(lk(in), lpb(in))), log(k(w)), 'linear', NaN);
    rmsd(i) = sqrt(mean((log10(exp(li)) - log10(ref(k(w),P(w)))).^2, 'omitnan'));
end
grid on; ylim([1e-3 30]); xlabel('vertical wavenumber (cycles/m)'); ylabel('normalised power (1 at 50 m)');
legend([{sprintf('56-32 sonic log (1/k^{%.2f})', -c0(1))}, cellfun(@(a,s,r) sprintf('%s (1/k^{%.2f}, rms %.2f)', a, s, r), nm, num2cell(slopes), num2cell(rmsd), 'uni', 0)], ...
    'Location','southwest','FontSize',8);
exportgraphics(gcf,'../figs/fig_log_vs_media.png','Resolution',150);
T = table(nm', slopes', rmsd', 'VariableNames', {'medium','slope_12_200m','rms_log10_vs_log'});
disp(T); fprintf('log slope %.2f\n', -c0(1));
writetable(T, '../results/table_log_vs_media.csv');
