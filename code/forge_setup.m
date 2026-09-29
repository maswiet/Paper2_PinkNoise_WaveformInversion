function G = forge_setup()
% FORGE_SETUP  Common FD geometry for the 56-32 sensor / 16A-stage MEQ volume.
% Local frame: E, N, depth (m) as in the GES catalogue (origin 16A wellhead).
load(fullfile(fileparts(mfilename('fullpath')),'..','data','tt_inversion.mat'),'TT','obs');
load(fullfile(fileparts(mfilename('fullpath')),'..','data','forge2022_events.mat'),'ev');
q = TT(1).q;                                   % homogeneous-model solution
G.sens = [q(1:3); q(1:2) q(4)];                % sensor A, B  (E N Z)
G.off  = q(5:6);                               % clock offsets (s)
G.vp = TT(1).par(1); G.vs = G.vp/1.73; G.rho = 2650;
G.dx = 5; G.nab = 20;
lo = [780 -420 2120]; hi = [1140 60 2660];     % interior box (m)
G.x0 = lo - G.nab*G.dx;                        % grid origin (node 1)
G.n  = round((hi - lo)/G.dx) + 2*G.nab + 1;
G.dt = 2.5e-4; G.nt = 1300; G.f0 = 45; G.t0 = 0.027;
% events: those with a clean pick on either sensor and inside the box
use = unique(obs(:,1));
x = ev.xyz(use,:);
in = all(x > lo+10 & x < hi-10, 2);
G.ev = use(in); G.evxyz = x(in,:);
G.ievg = round((G.evxyz - G.x0)/G.dx) + 1;     % event grid nodes
G.isg  = round((G.sens - G.x0)/G.dx) + 1;      % sensor grid nodes
G.z0 = G.x0(3);
end
