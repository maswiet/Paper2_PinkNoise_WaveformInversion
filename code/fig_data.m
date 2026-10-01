% Figures of the FORGE 2022 56-32 data: record sections (horizontal energy,
% 40-80 Hz, Ricker-convolved) and P travel-time residual correlation.
maxNumCompThreads(2);
load ../data/forge2022_events.mat; load ../data/forge_features.mat; load ../data/tt_inversion.mat
G = forge_setup(); fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
hch = {[2 3],[5 6]}; [b,a] = butter(3,[40 80]/(fs/2));
tr = (-0.05:1/fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
f = figure('Visible','off');
for g = 1:2
    k = D(g,1).k; R = D(g,1).R; [R,o] = sort(R); k = k(o); tp = D(g,1).tpl(o);
    subplot(2,3,3*g-2); hold on
    sel = round(linspace(1,numel(k),40));
    for q = sel
        x = filtfilt(b,a,conv2(double(ev.W(:,hch{g},k(q))),ric,'same'));
        plot(t - tp(q), x(:,1)/max(abs(x(:)))*8 + R(q), 'k', 'LineWidth',0.3);
    end
    plot(R/G.vs - R/G.vp, R, 'r--'); xlim([-0.02 0.2]); ylim([min(R)-10 max(R)+10]);
    xlabel('time after P (s)'); ylabel('hypocentral distance (m)'); panel_label(gca, 1 + 3*(g-1)); box on
    subplot(2,3,3*g-1);
    E = [D(g,2).F.envS]; lag = (-0.08:4/fs:0.12)';
    semilogy(lag, median(E,2,'omitnan'),'k','LineWidth',1.2); hold on
    semilogy(lag, prctile(E,25,2),'k:'); semilogy(lag, prctile(E,75,2),'k:');
    xlabel('time after S-envelope peak (s)'); ylabel('E/E_{S,max}'); ylim([1e-3 2]); grid on; panel_label(gca, 2 + 3*(g-1));
end
subplot(2,3,[3 6]);
e = TT(1).edges; c = (e(1:end-1)+e(2:end))/2;
plot(c, TT(1).rho{1}, 'o-', c, TT(1).rho{2}, 's-','LineWidth',1,'MarkerSize',4); grid on
set(gca,'XScale','log','XLim',[3 300],'XTick',[5 10 20 50 100 200]); xlabel('inter-event separation (m)'); ylabel('correlation of P residuals');
legend('sensor A','sensor B','Box','off'); panel_label(gca, 3);
pub_export(f, '../figs/fig_data_overview.png', 17.4, 11); close(f);
