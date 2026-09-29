function M = build_model(cls, P, g)
% BUILD_MODEL  Crustal model classes compared in Paper 2.
%   cls : 'homog' | 'layered' | 'vti' | 'pink'
%   P   : parameter struct (see below)
%   g   : grid struct .n=[nx ny nz] .dx  .z0 (depth of top of grid, m)
%
% homog   : P.vp, P.vs                            (uniform crust)
% layered : P.vp, P.vs at top, P.gz (1/s) vertical gradient, plus optional
%           P.zi (interface depths, m) and P.dv (fractional velocity jumps)
% vti     : homog/layered + P.eps, P.del, P.gam  (Thomsen, uniform)
% pink    : P.vp, P.vs mean + P.sig (fractional std of velocity),
%           P.p (3D filter exponent, 1.35 ~ 1/k well-log scaling), P.seed.
%           Velocity perturbation follows GFI: log-velocity ~ pink noise,
%           i.e. v = v0*exp(sig*f) (lognormal, like k ~ exp(alpha*phi)).
%           P.vpvs_coupled (default true): vp,vs share the same f.
n = g.n; dx = g.dx;
z = g.z0 + ((1:n(3))-1)*dx;              % depth of each z-slice
Z = reshape(z,1,1,[]);
one = ones(n,'single');

M.dx = dx; M.rho = 2650;
M.eps = 0; M.del = 0; M.gam = 0;
if ~isfield(P,'gz'), P.gz = 0; end

switch lower(cls)
    case 'homog'
        M.vp = P.vp*one; M.vs = P.vs*one;
    case {'layered','vti'}
        fz = 1 + P.gz*(Z - g.z0)/P.vp;       % linear gradient (same fraction for vs)
        if isfield(P,'zi')
            for i = 1:numel(P.zi), fz = fz .* (1 + P.dv(i)*(Z >= P.zi(i))); end
        end
        M.vp = single(P.vp*fz).*one; M.vs = single(P.vs*fz).*one;
        if strcmpi(cls,'vti')
            M.eps = P.eps; M.del = P.del; M.gam = P.gam;
        end
    case 'pink'
        f = single(pinknoise3d(n, P.p, P.seed));
        if isfield(P,'clip'), f = max(min(f, P.clip), -P.clip); end   % truncate tails (FD stability)
        M.vp = P.vp*exp(P.sig*f);
        if ~isfield(P,'vpvs_coupled') || P.vpvs_coupled
            M.vs = P.vs*exp(P.sig*f);
        else
            M.vs = P.vs*exp(P.sig*single(pinknoise3d(n, P.p, P.seed+1)));
        end
        if isfield(P,'gz') && P.gz ~= 0
            fz = single(1 + P.gz*(Z - g.z0)/P.vp);
            M.vp = M.vp.*fz; M.vs = M.vs.*fz;
        end
    otherwise
        error('unknown class %s', cls);
end
end
