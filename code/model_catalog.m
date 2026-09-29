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
    case 'T2'     % null-test truth: log layering + VTI (eps .06, del .03, gam .06)
        M = build_model('homog', struct('vp',vp,'vs',vs), g);
        M.vp = M.vp.*exp(dl); M.vs = M.vs.*exp(dl);
        M.eps = 0.06; M.del = 0.03; M.gam = 0.06;
    otherwise
        error('unknown model %s', name);
end
M.rho = G.rho;
info = sprintf('%s: vp %.0f-%.0f, vs %.0f-%.0f, eps %.3f del %.3f gam %.3f', name, ...
    min(M.vp(:)), max(M.vp(:)), min(M.vs(:)), max(M.vs(:)), M.eps, M.del, M.gam);
end
