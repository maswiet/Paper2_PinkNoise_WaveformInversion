% Figures: (1) FORGE geometry; (2) vertical slices of candidate models
maxNumCompThreads(2);
G = forge_setup();
load ../data/forge2022_events.mat
figure('Position',[50 50 1000 450]);
subplot(1,2,1);
scatter3(G.evxyz(:,1), G.evxyz(:,2), G.evxyz(:,3), 6, ev.stage(ismember((1:numel(ev.M))',G.ev)), 'filled'); hold on
plot3(G.sens(:,1), G.sens(:,2), G.sens(:,3), 'kv', 'MarkerFaceColor','r','MarkerSize',10);
text(G.sens(:,1)+15, G.sens(:,2), G.sens(:,3), {'A','B'},'FontWeight','bold');
lo = G.x0 + G.nab*G.dx; hi = G.x0 + (G.n-1-G.nab)*G.dx;
plotcube(lo, hi); set(gca,'ZDir','reverse'); axis equal; grid on
xlabel('East (m)'); ylabel('North (m)'); zlabel('Depth (m)'); view(-35,20);
panel_label(gca, 1);
subplot(1,2,2);
RA = sqrt(sum((G.evxyz-G.sens(1,:)).^2,2)); RB = sqrt(sum((G.evxyz-G.sens(2,:)).^2,2));
histogram(RA,30); hold on; histogram(RB,30); legend('sensor A','sensor B'); xlabel('hypocentral distance (m)'); ylabel('number of events');
panel_label(gca, 2);
exportgraphics(gcf,'../figs/fig_geometry.png','Resolution',150);

names = {'H','L','A','P4','P1','P3','P5','P2'};
ttl = {'H homogeneous','L layered (56-32 log)','A VTI (TT-inverted)','P4 pink \sigma=0.02 p=1.5', ...
       'P1 pink \sigma=0.045 p=1.5','P3 pink \sigma=0.09 p=1.5','P5 pink \sigma=0.045 p=1.2','P2 pink \sigma=0.045 p=1.8'};
figure('Position',[50 50 1400 650]);
x = G.x0(1) + (0:G.n(1)-1)*G.dx; z = G.x0(3) + (0:G.n(3)-1)*G.dx;
jy = G.isg(1,2);
for i = 1:numel(names)
    M = model_catalog(names{i}, G);
    subplot(2,4,i);
    v = M.vp; if isscalar(v), v = v*ones(G.n); end
    imagesc(x, z, squeeze(v(:,jy,:))'); axis image; colorbar; caxis([4800 7000]);
    hold on; plot(G.sens(:,1), G.sens(:,3), 'wv','MarkerFaceColor','w');
    title(ttl{i}); xlabel('E (m)'); ylabel('depth (m)');
end
colormap(turbo);
exportgraphics(gcf,'../figs/fig_model_slices.png','Resolution',150);

function plotcube(lo, hi)
[X,Y,Z] = ndgrid([lo(1) hi(1)],[lo(2) hi(2)],[lo(3) hi(3)]);
P = [X(:) Y(:) Z(:)];
E = [1 2;1 3;2 4;3 4;5 6;5 7;6 8;7 8;1 5;2 6;3 7;4 8];
for e = 1:size(E,1), plot3(P(E(e,:),1), P(E(e,:),2), P(E(e,:),3), 'k:'); end
end
