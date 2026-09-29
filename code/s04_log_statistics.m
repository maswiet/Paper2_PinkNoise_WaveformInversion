% S04  56-32 sonic-porosity log -> velocity fluctuation statistics and 1D
% layered model for the sensor/event depth interval (granitoid).
% SPHI is Wyllie sonic porosity; invert with Schlumberger defaults
% (DTma = 55.5 us/ft, DTf = 189 us/ft). Only relative fluctuations are used.
L = readmatrix('../../56-32/A.1008700.02.01_UniversityUtah_Forge-56-32_ThruBit_Main.las.dat','FileType','text','NumHeaderLines',0);
L(L==-999.25) = NaN;
z = L(:,1)*0.3048; sphi = L(:,14); cal = L(:,4);
DT = 55.5 + sphi*(189-55.5);  v = 304800./DT;
sel = z > 1950 & z < 2780 & isfinite(v) & cal < 10.5;   % sensor/event interval, gauge hole
z = z(sel); lv = log(v(sel));
fprintf('interval %.0f-%.0f m, n=%d, mean v=%.0f m/s\n', min(z), max(z), numel(z), exp(mean(lv)));
c = polyfit(z, lv, 1); d = lv - polyval(c, z);
fprintf('log-velocity trend: %.2e /m  (%.1f m/s per 100 m)\n', c(1), c(1)*exp(mean(lv))*100);
dz = median(diff(z));
sig_raw = std(d);
w6 = round(6/dz); d6 = movmean(d, w6); sig6 = std(d6);
w50 = round(50/dz); d50 = movmean(d, w50); sig50 = std(d50);
fprintf('sigma(lnV): raw(0.15 m) %.3f | 6 m avg %.3f | 50 m avg %.3f\n', sig_raw, sig6, sig50);
% power spectrum slope
n = numel(d); D = abs(fft(d - mean(d))).^2; k = (1:floor(n/2)-1)'/(n*dz);
P = D(2:floor(n/2)); in = k > 1/200 & k < 1/1;
cf = polyfit(log(k(in)), log(P(in)), 1);
fprintf('well-log spectrum S(k) ~ 1/k^%.2f over 1-200 m\n', -cf(1));
% synthetic check: sigma of 1D line through 3D pinknoise (p=1.35) at 6 m vs 50 m smoothing
f3 = pinknoise3d([128 128 128], 1.35, 7); l = squeeze(f3(64,64,:));
fprintf('3D pink p=1.35, 1D line: std(8-pt avg)/std(1-pt) = %.2f; log-spectrum %.2f vs log %.2f\n', ...
    std(movmean(l,8))/std(l), NaN, sig50/sig6);
LOG.z = z; LOG.d = d; LOG.d6 = d6; LOG.d50 = d50; LOG.trend = c; LOG.sig6 = sig6; LOG.sig50 = sig50;
LOG.beta = -cf(1); LOG.k = k; LOG.P = P;
save ../data/log5632_stats.mat LOG
figure('Position',[50 50 1100 600]);
subplot(1,3,1); plot(exp(lv),z,'Color',[.7 .7 .7]); hold on; plot(exp(polyval(c,z)+d50),z,'k','LineWidth',1.5);
set(gca,'YDir','reverse'); xlabel('V_P from SPHI (m/s)'); ylabel('depth (m)');
yline([2175 2523],'r--'); panel_label(gca, 1);
subplot(1,3,2); histogram(d6,60,'Normalization','pdf'); xlabel('\delta ln V (6 m avg)'); ylabel('probability density');
panel_label(gca, 2);
subplot(1,3,3); loglog(k,P,'.','Color',[.6 .6 .6]); hold on; loglog(k(in),exp(polyval(cf,log(k(in)))),'r','LineWidth',2);
xlabel('k (1/m)'); ylabel('S(k)'); panel_label(gca, 3, 'tr');
exportgraphics(gcf,'../figs/fig_log5632.png','Resolution',150);
