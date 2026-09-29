function S = synth_records(name, mech, addnoise)
% SYNTH_RECORDS  Horizontal-component synthetics at sensors A,B for every FORGE
% event node, from the model's strain Green's tensors (reciprocity).
%   mech : 'forge' (strike U(0,50), dip U(55,85), rake N(-90,30)) | 'random'
%   addnoise : true -> add real FORGE noise at the event's observed SNR
% S(g).H [nt x 2 x ne], S(g).t, S(g).k (catalogue index), S(g).R
if nargin < 2, mech = 'forge'; end
if nargin < 3, addnoise = true; end
G = forge_setup();
L = load(fullfile(fileparts(mfilename('fullpath')),'..','sgt',[name '.mat']),'E');
Df = load(fullfile(fileparts(mfilename('fullpath')),'..','data','forge_features.mat'),'D','NZ');
ne = numel(G.ev); t = (0:G.nt-1)'*G.dt - G.t0;
rng(1234);                                   % identical mechanisms for all models
switch mech
    case 'forge'
        str = 50*rand(ne,1); dip = 55 + 30*rand(ne,1); rk = -90 + 30*randn(ne,1);
    case 'random'
        str = 360*rand(ne,1); dip = acosd(rand(ne,1)); rk = 360*rand(ne,1) - 180;
end
fs = 1/G.dt;
[b,a] = butter(3,[20 80]/(fs/2));
for g = 1:2
    H = zeros(G.nt, 2, ne, 'single');
    for i = 1:ne
        Mt = dc_moment(str(i), dip(i), rk(i));
        for n = 1:2
            H(:,n,i) = single(cumsum(sgt_synth(L.E{g,n}(:,:,i), Mt))*G.dt);
        end
    end
    R = sqrt(sum((G.evxyz - G.sens(g,:)).^2, 2));
    if addnoise   % match each event's observed SNR with a real noise sample
        kd = Df.D(g,1).k; N = Df.NZ(g).N; nn = min(size(N,1), G.nt);
        for i = 1:ne
            j = find(kd == G.ev(i), 1);
            if isempty(j), j = randi(numel(kd)); end
            X = filtfilt(b, a, double(H(:,:,i)));
            tP = R(i)/G.vp;
            pe = max(sum(X(t >= tP-0.004 & t < tP+0.015,:).^2, 2));
            H(1:nn,:,i) = H(1:nn,:,i) + single(N(1:nn,:,j))*sqrt(pe);
        end
    end
    S(g).H = H; S(g).t = t; S(g).k = G.ev; S(g).R = R; %#ok<AGROW>
end
end
