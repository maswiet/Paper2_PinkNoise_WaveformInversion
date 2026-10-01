function Z = s11_shape_diagnostics(name, mech)
% S11  Shape-sensitive wavefield diagnostics (beyond coda strength):
%   D1 frequency dependence of the S coda: clevel2(40-80 Hz) - clevel2(20-40 Hz), sensor A, per event
%   D2 inter-event coherence of the S coda vs event separation (both sensors, 40-80 Hz):
%      max normalised 2-component cross-correlation (lag +-8 ms) of S+15..S+100 ms windows
% name : 'data' or a model name;  mech : 'forge' | 'fixed' (models only)
% Output: ../results/shape_<name>_<mech>.mat
if nargin < 2, mech = 'forge'; end
here = fileparts(mfilename('fullpath'));
G = forge_setup(); fs = 1/G.dt;
edges = [0 5 10 20 40 80 160];
rng(7);
if strcmp(name, 'data')
    Df = load(fullfile(here,'..','data','forge_features.mat'));
    ev = load(fullfile(here,'..','data','forge2022_events.mat')); ev = ev.ev;
    tr = (-0.05:1/ev.fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
    td = (0:size(ev.W,1)-1)'/ev.fs - ev.pre; hch = {[2 3],[5 6]};
    for g = 1:2
        k = Df.D(g,1).k; [~, loc] = ismember(k, G.ev);
        H = zeros(numel(td), 2, numel(k));
        for i = 1:numel(k), H(:,:,i) = conv2(double(ev.W(:,hch{g},k(i))), ric, 'same'); end
        REC(g).H = H; REC(g).t = td; REC(g).tp = Df.D(g,1).tpl; REC(g).xyz = G.evxyz(loc,:); %#ok<AGROW>
        REC(g).R = sqrt(sum((REC(g).xyz - G.sens(g,:)).^2,2)); REC(g).k = k; %#ok<AGROW>
    end
    Fa = {Df.D(1,1).F, Df.D(1,2).F}; ka = Df.D(1,1).k;
else
    S = synth_records(name, mech, true);
    f = fullfile(here,'..','results',sprintf('feat_%s_forge.mat', name));
    if strcmp(mech,'forge') && exist(f,'file'), L = load(f); tpA = {L.M(1,1).tp, L.M(2,1).tp}; else, tpA = {[],[]}; end
    for g = 1:2
        REC(g).H = double(S(g).H); REC(g).t = S(g).t; REC(g).xyz = G.evxyz; REC(g).R = S(g).R; REC(g).k = S(g).k; %#ok<AGROW>
        if isempty(tpA{g}), tpA{g} = S(g).R/G.vp; end
        REC(g).tp = tpA{g}; %#ok<AGROW>
    end
    if exist(f,'file') && strcmp(mech,'forge'), Fa = {L.M(1,1).F, L.M(1,2).F}; ka = L.M(1,1).k; else, Fa = {}; end
end
% ---------- D1 ----------
if ~isempty(Fa)
    Z.dlev = [Fa{2}.clevel2]' - [Fa{1}.clevel2]'; Z.dlev_k = ka;
else
    Z.dlev = []; Z.dlev_k = [];
end
% ---------- D2 ----------
[b,a] = butter(3,[40 80]/(fs/2));
Z.edges = edges; Z.coh = NaN(2, numel(edges)-1); Z.coh_q = NaN(2, numel(edges)-1, 2); Z.npair = zeros(2, numel(edges)-1);
Z.pairs = cell(2,1);
lagmax = round(0.008*fs);
for g = 1:2
    t = REC(g).t; n = size(REC(g).H,3);
    tS = REC(g).tp + REC(g).R*(1/G.vs - 1/G.vp);
    nw = round(0.085*fs); W = zeros(nw + 2*lagmax, 2, n); ok = false(n,1);
    for i = 1:n
        X = filtfilt(b, a, REC(g).H(:,:,i));
        i0 = find(t >= tS(i) + 0.015, 1) - lagmax;
        if isempty(i0) || i0 < 1 || i0 + size(W,1) - 1 > numel(t), continue; end
        w = X(i0:i0+size(W,1)-1, :);
        W(:,:,i) = w / (sqrt(mean(w(:).^2)) + realmin); ok(i) = true;   % scale-free (synthetics ~1e-9)
    end
    idx = find(ok); D = squareform(pdist(REC(g).xyz(idx,:)));
    [ii, jj] = find(triu(D > 0 & D <= edges(end), 1));
    if numel(ii) > 20000, s = randperm(numel(ii), 20000); ii = ii(s); jj = jj(s); end
    c = zeros(numel(ii),1);
    for q = 1:numel(ii)
        A = W(lagmax+1:lagmax+nw, :, idx(ii(q)));                % reference window
        B = W(:, :, idx(jj(q)));                                  % padded window
        na = sqrt(sum(A(:).^2)); best = 0;
        for L = 0:2*lagmax
            Bs = B(L+1:L+nw, :);
            v = sum(A(:).*Bs(:)) / (na*sqrt(sum(Bs(:).^2)) + eps);
            if abs(v) > best, best = abs(v); end
        end
        c(q) = best;
    end
    d = D(sub2ind(size(D), ii, jj));
    Z.pairs{g} = [d c];
    Z.pairidx{g} = [REC(g).k(idx(ii)) REC(g).k(idx(jj))];   % catalogue event ids of each pair
    for e = 1:numel(edges)-1
        m = d > edges(e) & d <= edges(e+1);
        Z.npair(g,e) = nnz(m);
        if nnz(m) >= 20
            Z.coh(g,e) = median(c(m)); Z.coh_q(g,e,:) = prctile(c(m), [25 75]);
        end
    end
end
save(fullfile(here,'..','results',sprintf('shape_%s_%s.mat', name, mech)), '-struct', 'Z');
fprintf('%s (%s): D1 median %.2f | coherence A %s | B %s\n', name, mech, median(Z.dlev,'omitnan'), ...
    mat2str(Z.coh(1,:),2), mat2str(Z.coh(2,:),2));
end
