% FIG_GEOMETRY_REV  Publication version of the FORGE 2022 geometry figure:
% a) map view and b) east-depth section of the 416 analysed MEQs (colour = stimulation stage; the
% other in-box events in grey), sensors A and B in well 56-32 and the FD interior; c) hypocentral
% distances of the 416 events to both sensors.
here = fileparts(mfilename('fullpath'));
G = forge_setup();
ev = load(fullfile(here,'..','data','forge2022_events.mat')); ev = ev.ev;
S = load(fullfile(here,'..','results','s19_H.mat'), 'ev');
use = ismember(G.ev, S.ev); st = ev.stage(G.ev);
lo = G.x0 + G.nab*G.dx; hi = G.x0 + (G.n-1-G.nab)*G.dx;
c2 = [0.13 0.40 0.75]; c3 = [0.90 0.45 0.10]; cg = [0.75 0.75 0.75];
f = figure('Visible','off');
% ---- a) map view
ax1 = subplot(1,3,1); hold on
plot(G.evxyz(~use,1), G.evxyz(~use,2), '.', 'Color', cg, 'MarkerSize', 4);
m2 = use & st == 2; m3 = use & st == 3;
h3 = plot(G.evxyz(m3,1), G.evxyz(m3,2), '.', 'Color', c3, 'MarkerSize', 5);
h2 = plot(G.evxyz(m2,1), G.evxyz(m2,2), '.', 'Color', c2, 'MarkerSize', 7);
rectangle('Position', [lo(1) lo(2) hi(1)-lo(1) hi(2)-lo(2)], 'LineStyle', ':', 'EdgeColor', 'k');
hs = plot(G.sens(1,1), G.sens(1,2), 'kv', 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerSize', 6);
text(G.sens(1,1)+12, G.sens(1,2)+25, '56-32', 'FontSize', 7);
axis equal; box on; grid on; xlim([lo(1)-60 hi(1)+20]); ylim([lo(2)-20 hi(2)+40]);
xlabel('East (m)'); ylabel('North (m)');
lg = legend([h3 h2 hs], {'stage 3', 'stage 2', 'sensors A, B'}, 'Location', 'southwest'); lg.Box = 'off'; lg.ItemTokenSize = [10 8];
panel_label(ax1, 1);
% ---- b) east-depth section
ax2 = subplot(1,3,2); hold on
plot(G.evxyz(~use,1), G.evxyz(~use,3), '.', 'Color', cg, 'MarkerSize', 4);
plot(G.evxyz(m3,1), G.evxyz(m3,3), '.', 'Color', c3, 'MarkerSize', 5);
plot(G.evxyz(m2,1), G.evxyz(m2,3), '.', 'Color', c2, 'MarkerSize', 7);
rectangle('Position', [lo(1) lo(3) hi(1)-lo(1) hi(3)-lo(3)], 'LineStyle', ':', 'EdgeColor', 'k');
plot(G.sens(:,1), G.sens(:,3), 'kv', 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerSize', 6);
text(G.sens(:,1)+14, G.sens(:,3), {'A','B'}, 'FontSize', 8, 'FontWeight', 'normal');
set(ax2, 'YDir', 'reverse'); axis equal; box on; grid on; xlim([lo(1)-60 hi(1)+20]); ylim([lo(3)-80 hi(3)+20]);
xlabel('East (m)'); ylabel('Depth (m)');
panel_label(ax2, 2);
% ---- c) distances
ax3 = subplot(1,3,3); hold on
RA = sqrt(sum((G.evxyz(use,:) - G.sens(1,:)).^2, 2)); RB = sqrt(sum((G.evxyz(use,:) - G.sens(2,:)).^2, 2));
be = 150:10:560;
histogram(RA, be, 'FaceColor', [0.30 0.30 0.30], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
histogram(RB, be, 'FaceColor', [0.40 0.70 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
box on; grid on; xlim([150 560]); ylim([0 75]); xlabel('hypocentral distance (m)'); ylabel('number of events');
lg = legend({'sensor A', 'sensor B'}, 'Location', 'northeast'); lg.Box = 'off'; lg.ItemTokenSize = [12 8];
panel_label(ax3, 3);
set(ax1, 'Position', [0.07 0.17 0.27 0.78]); set(ax2, 'Position', [0.40 0.17 0.24 0.78]); set(ax3, 'Position', [0.72 0.17 0.26 0.78]);
pub_export(f, fullfile(here,'..','figs','fig_geometry_rev.png'), 17.4, 7.2);
close(f);
fprintf('events used %d (stage 2: %d, stage 3: %d), other in-box %d\n', nnz(use), nnz(m2), nnz(m3), nnz(~use));
