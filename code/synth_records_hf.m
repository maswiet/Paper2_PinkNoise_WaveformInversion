function S = synth_records_hf(name)
% SYNTH_RECORDS_HF  Sensor-B horizontal synthetics from the high-frequency
% (dx 2.5 m, Ricker 90 Hz) strain Green's tensors; FORGE-like DC mechanisms,
% noise-free. S(2).H [nt x 2 x ne], S(2).t, S(2).k, S(2).R
G = forge_setup_hf();
L = load(fullfile(fileparts(mfilename('fullpath')),'..','sgt_hf',[name '.mat']),'E');
ne = numel(G.ev); dt = G.dt*G.decim; nt = size(L.E{2,1},1);
t = (0:nt-1)'*dt - G.t0;
rng(1234);
str = 50*rand(ne,1); dip = 55 + 30*rand(ne,1); rk = -90 + 30*randn(ne,1);
H = zeros(nt, 2, ne, 'single');
for i = 1:ne
    Mt = dc_moment(str(i), dip(i), rk(i));
    for n = 1:2
        H(:,n,i) = single(cumsum(sgt_synth(L.E{2,n}(:,:,i), Mt))*dt);
    end
end
S(2).H = H; S(2).t = t; S(2).k = G.ev; S(2).R = sqrt(sum((G.evxyz - G.sens(2,:)).^2, 2));
S(1).H = []; S(1).t = t; S(1).k = []; S(1).R = [];
end
