% Validation: reciprocity (SGT) vs direct moment-tensor simulation, and
% heterogeneous (pink) medium to check robustness of the approach.
n = [80 80 90]; dx = 5;
g.n = n; g.dx = dx; g.z0 = 0;
Mh = build_model('pink', struct('vp',5800,'vs',3350,'sig',0.045,'p',1.5,'seed',5), g); Mh.rho = 2650;
Mh.eps = 0.05; Mh.del = 0.02; Mh.gam = 0.05;        % also exercise VTI
xs = [30 35 30]; xe = [52 45 62];                  % sensor, event nodes
str = 30; dip = 70; rk = -80;                      % DC mechanism (Aki & Richards)
Mt = dc_moment(str, dip, rk);
opt.dt = 2.5e-4; opt.nt = 700; opt.nabs = 20;
src = struct('ix',xe(1),'iy',xe(2),'iz',xe(3),'f0',45,'t0',0.027,'type','dc','M',Mt);
tic; d = fd3d_elastic_vti(Mh, src, xs, opt); t1 = toc;
vd = [d.vx d.vy d.vz];
opt.erec = xe; vr = zeros(opt.nt,3);
for k = 1:3
    s = fd3d_elastic_vti(Mh, struct('ix',xs(1),'iy',xs(2),'iz',xs(3),'f0',45,'t0',0.027,'type','force','dir',k), xe, opt);
    vr(:,k) = sgt_synth(s.e(:,:,1), Mt);
end
vr = cumsum(vr)*opt.dt;                             % SGT gives time derivative
t = (0:opt.nt-1)'*opt.dt;
for k = 1:3
    c = vd(:,k)\vr(:,k);
    fprintf('comp %d: scale %.3g, corr %.4f, rel. misfit %.3f\n', k, c, corr(vd(:,k),vr(:,k)), ...
        norm(vd(:,k)-vr(:,k)/c)/norm(vd(:,k)));
end
fprintf('direct run %.1f s for %d cells\n', t1, prod(n));
figure('Position',[50 50 900 650]); cmp = {'v_E','v_N','v_Z'};
for k = 1:3
    subplot(3,1,k); plot(t, vd(:,k)/max(abs(vd(:))), 'k', t, vr(:,k)/max(abs(vr(:))), 'r--', 'LineWidth', 1);
    ylabel(sprintf('%s (normalised)', cmp{k})); panel_label(gca, k); grid on
    if k == 3, xlabel('time (s)'); end
end
legend('direct','reciprocity (SGT)', 'Location','southeast');
exportgraphics(gcf,'../figs/val_reciprocity.png','Resolution',150);
