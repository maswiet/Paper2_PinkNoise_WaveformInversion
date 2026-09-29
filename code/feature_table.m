function T = feature_table(F, R, tp, k, xyz)
% FEATURE_TABLE  Collect per-record features + travel-time residual stats
% for one sensor/band from a struct array F (wf_features output).
nm = {'pcoda','spd','cdecay','clevel','sp','clevel2','cdecay2','swidth'};
for j = 1:numel(nm), T.(nm{j}) = [F.(nm{j})]'; end
T.envS = [F.envS]; T.envP = [F.envP];
% P travel-time residual w.r.t. best straight-ray homogeneous fit (tp = a + R/v)
A = [ones(numel(R),1) R]; c = A\tp; r = tp - A*c;
r = r(abs(r - median(r)) < 5*1.4826*mad(r,1));
T.ttres_std = 1.4826*mad(r,1);
T.tt_r = tp - A*c;
% residual correlation vs inter-event separation
if nargin >= 5
    D = squareform(pdist(xyz)); rr = T.tt_r - median(T.tt_r);
    ok = abs(rr) < 5*1.4826*mad(rr,1); rr(~ok) = NaN;
    edges = [0 10 20 40 80 160];
    T.ttcorr = NaN(1,numel(edges)-1);
    for e = 1:numel(edges)-1
        [a,b] = find(triu(D > edges(e) & D <= edges(e+1), 1));
        v = isfinite(rr(a)) & isfinite(rr(b));
        if nnz(v) > 30, T.ttcorr(e) = corr(rr(a(v)), rr(b(v))); end
    end
    T.edges = edges;
end
T.k = k;
end
