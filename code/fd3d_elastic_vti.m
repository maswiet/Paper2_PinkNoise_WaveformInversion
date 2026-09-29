function [seis, t, snap] = fd3d_elastic_vti(M, src, rec, opt)
% FD3D_ELASTIC_VTI  3D velocity-stress staggered-grid elastic FD (O(2,4))
% with vertical transverse isotropy (VTI, Thomsen parameters).
% Staggered-grid treatment of anisotropy after Igel, Mora & Riollet (1995,
% Geophysics 60:1203-1216); developed from H. Igel's acoustic FD MATLAB codes
% (Igel 2016, Computational Seismology, Oxford University Press).
%
%   M   : struct with fields (arrays size [nx ny nz] or scalars)
%         vp, vs, rho            - vertical P/S velocity (m/s), density (kg/m3)
%         eps, del, gam          - Thomsen parameters (optional, default 0)
%         dx                     - grid spacing (m), isotropic grid
%   src : struct  .ix .iy .iz (grid indices), .f0 (Hz), .t0 (s),
%                 .type 'explosion' | 'dc' ; for 'dc' give .M (3x3 moment tensor)
%   rec : [nr x 3] receiver grid indices (ix iy iz)
%   opt : .nt, .dt, .nabs (sponge width), .snapit (optional step for snapshot)
%
%   seis: struct .vx .vy .vz  [nt x nr] particle velocity
%
% Grid: vx(i+1/2,j,k) vy(i,j+1/2,k) vz(i,j,k+1/2); sxx,syy,szz(i,j,k);
%       sxy(i+1/2,j+1/2,k) sxz(i+1/2,j,k+1/2) syz(i,j+1/2,k+1/2).
% Media parameters are taken at the integer node (no averaging) - adequate
% for the smooth band-limited random media used here.

[nx,ny,nz] = size(M.vp);
if isscalar(M.vp), error('M.vp must be a 3D array'); end
dx = M.dx; dt = opt.dt; nt = opt.nt;
if ~isfield(M,'eps'), M.eps = 0; end
if ~isfield(M,'del'), M.del = 0; end
if ~isfield(M,'gam'), M.gam = 0; end
cls = 'single';

rho = cast(M.rho .* ones(nx,ny,nz), cls);
C33 = rho .* cast(M.vp,cls).^2;
C44 = rho .* cast(M.vs,cls).^2;
C11 = C33 .* (1 + 2*cast(M.eps,cls));
C66 = C44 .* (1 + 2*cast(M.gam,cls));
C13 = sqrt(max((C33-C44).*(C33.*(1+2*cast(M.del,cls))-C44),0)) - C44;
C12 = C11 - 2*C66;
b   = 1 ./ rho;

% stability check (Courant, O(4) staggered: sum|c| = 7/6)
vmax = max(sqrt(C11(:)./rho(:)));
cfl = vmax*dt/dx*sqrt(3)*7/6;
if cfl > 1, error('Unstable: CFL=%.3f', cfl); end

% sponge (Cerjan) damping
nab = opt.nabs; a = 0.3/nab;           % edge factor exp(-0.09) per step
d = @(n) cast(exp(-(a*max([nab:-1:1, zeros(1,n-2*nab), 1:nab]',0)).^2), cls);
g = reshape(d(nx),[],1,1) .* reshape(d(ny),1,[],1) .* reshape(d(nz),1,1,[]);

z0 = zeros(nx,ny,nz,cls);
vx=z0; vy=z0; vz=z0; sxx=z0; syy=z0; szz=z0; sxy=z0; sxz=z0; syz=z0;

c1 = cast(9/8,cls); c2 = cast(-1/24,cls); rdx = cast(1/dx,cls); dtc = cast(dt,cls);

% source time function: Ricker, moment-rate
t = (0:nt-1)'*dt;
arg = (pi*src.f0*(t-src.t0)).^2;
stf = cast((1-2*arg).*exp(-arg), cls);
sc = cast(1/dx^3, cls);
if ~isfield(src,'type'), src.type = 'explosion'; end
isforce = strcmpi(src.type,'force');
if strcmpi(src.type,'explosion'), Mt = eye(3); elseif ~isforce, Mt = src.M; end

nr = size(rec,1);
ir = sub2ind([nx ny nz], rec(:,1), rec(:,2), rec(:,3));
seis.vx = zeros(nt,nr,cls); seis.vy = seis.vx; seis.vz = seis.vx;
% optional strain-rate recording (for reciprocity / strain Green's tensor)
dostr = isfield(opt,'erec') && ~isempty(opt.erec);
if dostr
    ie = sub2ind([nx ny nz], opt.erec(:,1), opt.erec(:,2), opt.erec(:,3));
    seis.e = zeros(nt, 6, numel(ie), cls);   % exx eyy ezz exy exz eyz (rates)
end
snap = {};
si = src.ix; sj = src.iy; sk = src.iz;

for it = 1:nt
    % ---- velocity update ----
    vx = vx + dtc*b.*(Dxf(sxx) + Dyb(sxy) + Dzb(sxz));
    vy = vy + dtc*b.*(Dxb(sxy) + Dyf(syy) + Dzb(syz));
    vz = vz + dtc*b.*(Dxb(sxz) + Dyb(syz) + Dzf(szz));
    if isforce   % point force (unit impulse-response x Ricker), on staggered v node
        fs_ = dtc*stf(it)*sc/rho(si,sj,sk);
        switch src.dir
            case 1, vx(si,sj,sk) = vx(si,sj,sk) + fs_;
            case 2, vy(si,sj,sk) = vy(si,sj,sk) + fs_;
            case 3, vz(si,sj,sk) = vz(si,sj,sk) + fs_;
        end
    end
    vx = vx.*g; vy = vy.*g; vz = vz.*g;

    % ---- stress update ----
    exx = Dxb(vx); eyy = Dyb(vy); ezz = Dzb(vz);
    sxx = sxx + dtc*(C11.*exx + C12.*eyy + C13.*ezz);
    syy = syy + dtc*(C12.*exx + C11.*eyy + C13.*ezz);
    szz = szz + dtc*(C13.*exx + C13.*eyy + C33.*ezz);
    exy = Dyf(vx) + Dxf(vy); exz = Dzf(vx) + Dxf(vz); eyz = Dzf(vy) + Dyf(vz);
    sxy = sxy + dtc*C66.*exy;
    sxz = sxz + dtc*C44.*exz;
    syz = syz + dtc*C44.*eyz;
    if dostr   % shear rates sampled at nearest staggered node (half-cell offset)
        seis.e(it,:,:) = [exx(ie) eyy(ie) ezz(ie) 0.5*exy(ie) 0.5*exz(ie) 0.5*eyz(ie)]';
    end

    % ---- moment-tensor source (stress injection) ----
    if ~isforce
    s = dtc*stf(it)*sc;
    sxx(si,sj,sk) = sxx(si,sj,sk) - s*Mt(1,1);
    syy(si,sj,sk) = syy(si,sj,sk) - s*Mt(2,2);
    szz(si,sj,sk) = szz(si,sj,sk) - s*Mt(3,3);
    sxy(si,sj,sk) = sxy(si,sj,sk) - s*Mt(1,2);
    sxz(si,sj,sk) = sxz(si,sj,sk) - s*Mt(1,3);
    syz(si,sj,sk) = syz(si,sj,sk) - s*Mt(2,3);
    end
    sxx=sxx.*g; syy=syy.*g; szz=szz.*g; sxy=sxy.*g; sxz=sxz.*g; syz=syz.*g;

    seis.vx(it,:) = vx(ir); seis.vy(it,:) = vy(ir); seis.vz(it,:) = vz(ir);
    if isfield(opt,'snapit') && mod(it,opt.snapit)==0
        snap{end+1} = sqrt(vx.^2+vy.^2+vz.^2); %#ok<AGROW>
    end
end

% ---------- 4th-order staggered derivatives ----------
    function D = Dxf(f)   % forward: result at i+1/2
        D = zeros(size(f),cls);
        D(2:end-2,:,:) = (c1*(f(3:end-1,:,:)-f(2:end-2,:,:)) + c2*(f(4:end,:,:)-f(1:end-3,:,:)))*rdx;
    end
    function D = Dxb(f)   % backward: result at i
        D = zeros(size(f),cls);
        D(3:end-1,:,:) = (c1*(f(3:end-1,:,:)-f(2:end-2,:,:)) + c2*(f(4:end,:,:)-f(1:end-3,:,:)))*rdx;
    end
    function D = Dyf(f)
        D = zeros(size(f),cls);
        D(:,2:end-2,:) = (c1*(f(:,3:end-1,:)-f(:,2:end-2,:)) + c2*(f(:,4:end,:)-f(:,1:end-3,:)))*rdx;
    end
    function D = Dyb(f)
        D = zeros(size(f),cls);
        D(:,3:end-1,:) = (c1*(f(:,3:end-1,:)-f(:,2:end-2,:)) + c2*(f(:,4:end,:)-f(:,1:end-3,:)))*rdx;
    end
    function D = Dzf(f)
        D = zeros(size(f),cls);
        D(:,:,2:end-2) = (c1*(f(:,:,3:end-1)-f(:,:,2:end-2)) + c2*(f(:,:,4:end)-f(:,:,1:end-3)))*rdx;
    end
    function D = Dzb(f)
        D = zeros(size(f),cls);
        D(:,:,3:end-1) = (c1*(f(:,:,3:end-1)-f(:,:,2:end-2)) + c2*(f(:,:,4:end)-f(:,:,1:end-3)))*rdx;
    end
end
