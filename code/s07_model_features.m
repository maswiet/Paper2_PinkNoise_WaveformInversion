function s07_model_features(name, mech)
% S07  Features (same as FORGE data) for synthetic records of one model,
% plus low-band P travel-time residuals. Output ../results/feat_<name>_<mech>.mat
if nargin < 2, mech = 'forge'; end
G = forge_setup(); fs = 1/G.dt;
S = synth_records(name, mech, true);
bands = [20 40; 40 80];
M = struct();
for g = 1:2
    t = S(g).t; ne = size(S(g).H,3);
    tp = zeros(ne,1);
    for i = 1:ne
        tth = S(g).R(i)/G.vp;
        tp(i) = aic_pick(S(g).H(:,:,i), t, [20 100], fs, tth-0.02, tth+0.02);
    end
    dSP = S(g).R*(1/G.vs - 1/G.vp);
    for ib = 1:2
        for i = 1:ne
            M(g,ib).F(i) = wf_features(double(S(g).H(:,:,i)), t, tp(i), dSP(i), bands(ib,:), fs);
        end
        M(g,ib).k = S(g).k; M(g,ib).R = S(g).R; M(g,ib).tp = tp;
    end
end
out = fullfile(fileparts(mfilename('fullpath')),'..','results');
if ~exist(out,'dir'), mkdir(out); end
save(fullfile(out, sprintf('feat_%s_%s.mat', name, mech)), 'M', 'bands');
fprintf('features done: %s (%s)\n', name, mech);
end
