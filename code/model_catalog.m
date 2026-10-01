function [M, info] = model_catalog(name, G)
% MODEL_CATALOG  Named earth models for Paper 2 (see README for rationale).
g.n = G.n; g.dx = G.dx; g.z0 = G.z0;
vp = G.vp; vs = G.vs;
L = load(fullfile(fileparts(mfilename('fullpath')),'..','data','log5632_stats.mat'));
L = L.LOG;
% 1D layered fluctuation from 56-32 sonic log, 50-m running blocks, on grid depths
zg = G.z0 + ((1:G.n(3))-1)*G.dx;
dl = interp1(L.z, L.d50, zg, 'linear', 'extrap');
dl = reshape(single(dl - mean(dl)),1,1,[]);
switch name
    case 'H'      % homogeneous isotropic
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
    case 'L'      % layered: 56-32 log 50-m blocks
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        M.vp = M.vp.*exp(dl); M.vs = M.vs.*exp(dl);
    case 'A'      % VTI homogeneous, travel-time inverted Thomsen eps, del (gam = 0)
        load(fullfile(fileparts(mfilename('fullpath')),'..','data','tt_inversion.mat'),'TT');
        p = TT(3).par;
        M = build_model('homog', struct('vp',p(1),'vs',p(1)/1.73), g);
        M.eps = p(2); M.del = p(3); M.gam = 0;
    case {'P1','P2','P3','P4','P5','T1','P6','R1','R2','R3','P7','P8','P9','P10','P11','P12','R5','R6','R6c','R6d'}
        par = struct('P1',[0.045 1.5 11], 'P2',[0.045 1.8 11], 'P3',[0.09 1.5 11], ...
                     'P4',[0.02 1.5 11], 'P5',[0.045 1.2 11], 'T1',[0.045 1.5 22], ...
                     'P6',[0.13 1.5 11], 'R1',[0.045 1.5 33], 'R2',[0.045 1.5 44], ...
                     'R3',[0.09 1.5 33], 'P7',[0.09 1.2 11], 'P8',[0.09 1.8 11], 'P9',[0.18 1.5 11], ...
                     'P10',[0.25 1.5 11], 'P11',[0.18 1.8 11], 'P12',[0.22 1.5 11], 'R5',[0.18 1.5 33], 'R6',[0.13 1.5 33], 'R6c',[0.13 1.5 44], 'R6d',[0.13 1.5 55]);
        a = par.(name);
        P = struct('vp',vp,'vs',vs,'sig',a(1),'p',a(2),'seed',a(3));
        if a(1) > 0.1, P.clip = 3; end
        if a(1) > 0.15, P.clip = 2.5; end
        M = build_model('pink', P, g);
    case {'E15','E50','G15','G50','E15b','G50b','E15c','G50c','E15d','G50d'}   % correlation-length media, sigma as P6; suffix b/c/d = seed 33/44/55
        acf = struct('E','exp','G','gauss'); a = str2double(name(2:3));
        sd = 11; if numel(name) > 3, sd = struct('b',33,'c',44,'d',55); sd = sd.(name(4)); end
        f = single(randmedium3d(G.n, G.dx, acf.(name(1)), a, sd));
        f = max(min(f, 3), -3);                    % same tail truncation as P6
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        M.vp = vp*exp(0.13*f); M.vs = vs*exp(0.13*f);
    case {'EX5','EX9','EX17','EX34','EX67','EX134'}   % exponential ACF, correlation-length series for the Qc break test
        a = str2double(name(3:end));
        f = single(randmedium3d(G.n, G.dx, 'exp', a, 11));
        f = max(min(f, 3), -3);
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        M.vp = vp*exp(0.13*f); M.vs = vs*exp(0.13*f);
    case 'L6'     % layered: full-resolution 56-32 log (6-m average) sampled on the 5-m grid (revision)
        d6 = interp1(L.z, L.d6, zg, 'linear', 'extrap'); d6 = reshape(single(d6 - mean(d6)),1,1,[]);
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        M.vp = M.vp.*exp(d6); M.vs = M.vs.*exp(d6);
    case 'GRAD3'  % smooth 3-D gradients: +-3 % vertical, +-2 % in E and N across the box (revision)
        [x, y, z] = ndgrid(linspace(-1,1,G.n(1)), linspace(-1,1,G.n(2)), linspace(-1,1,G.n(3)));
        fz = single(1 + 0.03*z + 0.02*x - 0.02*y);
        M = build_model('homog', struct('vp',vp,'vs',vs), g); M.vp = M.vp.*fz; M.vs = M.vs.*fz;
    case 'FZ'     % deterministic fracture zones through the event volume (revision):
        % 12 planar zones, strike N25E, dip 75 deg, one cell thick, Vp -10 %, Vs -20 %
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        [ix, iy, iz] = ndgrid(1:G.n(1), 1:G.n(2), 1:G.n(3));
        X = G.x0(1) + (ix-1)*G.dx; Y = G.x0(2) + (iy-1)*G.dx; Z3 = G.x0(3) + (iz-1)*G.dx;
        stk = deg2rad(25); dp = deg2rad(75);
        nrm = [cos(stk)*sin(dp), -sin(stk)*sin(dp), cos(dp)];   % normal of plane striking N25E, dipping 75 deg
        rng(11); c0 = [900 -400 2300] + rand(12,3).*[220 400 350];
        for q = 1:12
            dist = abs((X - c0(q,1))*nrm(1) + (Y - c0(q,2))*nrm(2) + (Z3 - c0(q,3))*nrm(3));
            in = dist <= G.dx/2 & abs(Z3 - c0(q,3)) < 120;         % ~240-m tall zones
            M.vp(in) = 0.90*vp; M.vs(in) = 0.80*vs;
        end
    case 'VSD'    % pink sigma .13 p 1.5, independent Vp and Vs fields (Vp/Vs varies in space) (revision)
        P = struct('vp',vp,'vs',vs,'sig',0.13,'p',1.5,'seed',11,'clip',3,'vpvs_coupled',false);
        M = build_model('pink', P, g);
        M.vs = max(min(M.vs, vs*exp(0.39)), vs*exp(-0.39));   % same truncation as Vp
    case 'VSO'    % pink, Vs-dominant: sigma_Vs .13, sigma_Vp .065, same field (revision)
        f = single(pinknoise3d(G.n, 1.5, 11)); f = max(min(f, 3), -3);
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        M.vp = vp*exp(0.065*f); M.vs = vs*exp(0.13*f);
    case 'T2'     % null-test truth: log layering + VTI (eps .06, del .03, gam .06)
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        M.vp = M.vp.*exp(dl); M.vs = M.vs.*exp(dl);
        M.eps = 0.06; M.del = 0.03; M.gam = 0.06;
    otherwise
        if startsWith(name, 'X_')   % fair family grid (revision): X_<pink|exp|gau>_<sigma>_<p or a>_<seed>
            t = strsplit(name, '_'); fam = t{2}; sg = str2double(t{3}); pr = str2double(t{4}); sd = str2double(t{5});
            clip = 3; if sg > 0.15, clip = 2.5; end
            switch fam
                case 'pink', f = single(pinknoise3d(G.n, pr, sd));
                case 'exp',  f = single(randmedium3d(G.n, G.dx, 'exp', pr, sd));
                case 'gau',  f = single(randmedium3d(G.n, G.dx, 'gauss', pr, sd));
                otherwise, error('unknown family %s', fam);
            end
            f = max(min(f, clip), -clip);
            M = build_model('homog', struct('vp',vp,'vs',vs), g);
            M.vp = vp*exp(sg*f); M.vs = vs*exp(sg*f);
        else
            error('unknown model %s', name);
        end
end
M.rho = G.rho;
info = sprintf('%s: vp %.0f-%.0f, vs %.0f-%.0f, eps %.3f del %.3f gam %.3f', name, ...
    min(M.vp(:)), max(M.vp(:)), min(M.vs(:)), max(M.vs(:)), M.eps, M.del, M.gam);
end
