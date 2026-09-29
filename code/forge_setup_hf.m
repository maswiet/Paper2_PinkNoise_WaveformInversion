function G = forge_setup_hf()
% FORGE_SETUP_HF  High-frequency FD geometry (dx = 2.5 m, 40-160 Hz analysis):
% sensor B only and the core of the event cloud, so the grid stays ~5e6 cells.
G = forge_setup();                              % sensors, velocities, all events
G.dx = 2.5; G.nab = 20;
lo = [790 -390 2350]; hi = [1110 30 2630];      % interior box (m): sensor B + event core
G.x0 = lo - G.nab*G.dx;
G.n  = round((hi - lo)/G.dx) + 2*G.nab + 1;
G.dt = 1.25e-4; G.nt = 2200;                    % 0.275 s
G.f0 = 90; G.t0 = 0.015;
G.decim = 2;                                    % records decimated to 0.25 ms (data sampling)
in = all(G.evxyz > lo+8 & G.evxyz < hi-8, 2);
G.ev = G.ev(in); G.evxyz = G.evxyz(in,:);
G.ievg = round((G.evxyz - G.x0)/G.dx) + 1;
G.isg  = round((G.sens - G.x0)/G.dx) + 1;       % row 1 (sensor A) lies outside and is not used
G.sensors = 2;
G.z0 = G.x0(3);
end
