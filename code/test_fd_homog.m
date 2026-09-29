% Validation: homogeneous explosion, compare FD P arrival & timing
n = 96; dx = 10;
M.vp = 5800*ones(n,n,n,'single'); M.vs = 3350*ones(n,n,n,'single');
M.rho = 2650; M.dx = dx;
src.ix=48; src.iy=48; src.iz=48; src.f0=40; src.t0=0.03; src.type='explosion';
rec = [48+(5:5:40)', 48*ones(8,1), 48*ones(8,1)];
opt.dt = 0.8e-3; opt.nt = 300; opt.nabs = 20;
tic; [s,t] = fd3d_elastic_vti(M,src,rec,opt); el = toc;
fprintf('time %.1f s, %.3f s/step\n', el, el/opt.nt);
r = (5:5:40)*dx;
[~,ipk] = max(abs(s.vx)); tpk = t(ipk)';
% far-field P velocity ~ d/dt of Ricker -> peak |v| near t0 + r/vp (+/- ~quarter period)
disp([r' tpk' (src.t0 + r/5800)']);
amp = max(abs(s.vx))'; disp('amp*r (should be ~const, far field):'); disp((amp.*r')'/amp(end)/r(end));
