function run_convergence(cs, nthreads)
% RUN_CONVERGENCE  Numerical convergence of coda measures (reviewer point 12).
% One band-limited random medium (pink-like power spectrum (k0^2+k^2)^(-p/2), p = 1.5,
% k0 = 1/780 m, low-passed to wavelengths >= 12.5 m so that the 5-m grid resolves it;
% sigma = 0.13, truncated at +-3) is defined on a 2.5-m master grid and sampled onto:
%   'C5'   5-m grid,   test box, 20-node (100 m) sponge
%   'C25'  2.5-m grid, test box, 40-node (100 m) sponge          -> grid resolution
%   'CBIG' 5-m grid,   test box enlarged by 150 m on all sides     -> domain size
%   'CSP'  5-m grid,   test box, 40-node (200 m) sponge            -> absorbing boundary
% Sensor B and the events are snapped to the common 5-m lattice. Horizontal forces only.
% Output: ../conv/<cs>.mat  (E{2,1:2}, event list, geometry)
if nargin > 1, maxNumCompThreads(nthreads); end
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'conv'); if ~exist(out, 'dir'), mkdir(out); end
G0 = forge_setup();
IBlo = [800 -330 2380]; IBhi = [1080 40 2640];            % interior test box (m)
xm0 = [550 -580 2130]; xm1 = [1330 290 2890]; dxm = 2.5;   % master grid (covers all cases)
% ---- master field (built once, cached) ----
ff = fullfile(out, 'master_field.mat');
if exist(ff, 'file'), load(ff, 'F'); else
    nm = round((xm1 - xm0)/dxm) + 1; rng(77);
    kk = @(m) [0:ceil(m/2)-1, -floor(m/2):-1]/(m*dxm);      % cycles/m
    [kx, ky, kz] = ndgrid(kk(nm(1)), kk(nm(2)), kk(nm(3)));
    k = sqrt(kx.^2 + ky.^2 + kz.^2); clear kx ky kz
    A = (1/780^2 + k.^2).^(-1.5/2);
    taper = ones(size(k)); t = k > 0.06 & k < 0.08; taper(t) = 0.5*(1 + cos(pi*(k(t) - 0.06)/0.02)); taper(k >= 0.08) = 0;
    F = real(ifftn(fftn(randn(nm)) .* A .* taper)); clear A taper k
    F = single((F - mean(F(:)))/std(F(:)));
    F = max(min(F, 3), -3);
    save(ff, 'F', '-v7.3');
end
if strcmp(cs, 'BUILDONLY'), return; end
switch cs
    case 'C5',   dx = 5;   nab = 20; lo = IBlo;       hi = IBhi;
    case 'C25',  dx = 2.5; nab = 40; lo = IBlo;       hi = IBhi;
    case 'CBIG', dx = 5;   nab = 20; lo = IBlo - 150; hi = IBhi + 150;
    case 'CSP',  dx = 5;   nab = 40; lo = IBlo;       hi = IBhi;
end
x0 = lo - nab*dx; n = round((hi - lo)/dx) + 2*nab + 1;
% ---- sample onto this case's grid ----
step = round(dx/dxm); o = round((x0 - xm0)/dxm);
ix = o(1) + 1 + (0:n(1)-1)*step; iy = o(2) + 1 + (0:n(2)-1)*step; iz = o(3) + 1 + (0:n(3)-1)*step;
f = F(ix, iy, iz); clear F
M.dx = dx; M.rho = 2650; M.vp = G0.vp*exp(0.13*f); M.vs = G0.vs*exp(0.13*f); M.eps = 0; M.del = 0; M.gam = 0;
% ---- geometry on the common 5-m lattice ----
snap = @(x) round((x - xm0)/5)*5 + xm0;
sB = snap(G0.sens(2,:));
evx = snap(G0.evxyz); in = all(evx > IBlo + 10 & evx < IBhi - 10, 2);
evx = evx(in,:); ev = G0.ev(in);
isg = round((sB - x0)/dx) + 1; ievg = round((evx - x0)/dx) + 1;
dt = 2.5e-4*dx/5; nt = round(0.30/dt); decim = round(2.5e-4/dt);
vmax = max(M.vp(:)); cfl = vmax*dt/dx*sqrt(3)*7/6;
fprintf('%s: grid %s (%.1f M cells), %d events, CFL %.2f, sigma(lnV) %.3f\n', cs, mat2str(n), prod(n)/1e6, numel(ev), cfl, std(0.13*double(f(:))));
opt.dt = dt; opt.nt = nt; opt.nabs = nab; opt.erec = ievg;
E = cell(2,2);
for nf = 1:2
    src = struct('ix',isg(1),'iy',isg(2),'iz',isg(3),'f0',45,'t0',0.027,'type','force','dir',nf);
    tic; s = fd3d_elastic_vti(M, src, isg, opt);
    E{2,nf} = s.e(1:decim:end,:,:);
    fprintf('%s force %d: %.0f s\n', cs, nf, toc);
end
C.sens = [NaN NaN NaN; sB]; C.evxyz = evx; C.ev = ev; C.dt = 2.5e-4; C.t0 = 0.027; C.vp = G0.vp; C.vs = G0.vs;
save(fullfile(out, [cs '.mat']), 'E', 'C', '-v7.3');
end
