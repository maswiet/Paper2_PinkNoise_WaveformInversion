function W = s08_mt_waveform_inversion(obsname, models, nev)
% S08  Deterministic waveform inversion: per-event moment tensor (6 elements,
% linear) using each candidate model's strain Green's tensors; 4 horizontal
% traces (2 sensors x 2 comps), 20-80 Hz, window P-10 ms ... S+60 ms,
% per-sensor time shift (+-6 ms) for location/clock error.
% For FORGE data the tool azimuth/handedness of each sensor is grid-searched
% with model H. Output: VR per event/model, residual (unexplained) envelopes.
if nargin < 3, nev = 150; end
G = forge_setup(); fs = 1/G.dt; t = (0:G.nt-1)'*G.dt - G.t0;
[b,a] = butter(3,[20 80]/(fs/2));
R = [sqrt(sum((G.evxyz - G.sens(1,:)).^2,2)) sqrt(sum((G.evxyz - G.sens(2,:)).^2,2))];
% ---------------- observed records ----------------
if strcmp(obsname,'data')
    ev = load(fullfile(fileparts(mfilename('fullpath')),'..','data','forge2022_events.mat')); ev = ev.ev;
    Df = load(fullfile(fileparts(mfilename('fullpath')),'..','data','forge_features.mat'));
    tr = (-0.05:1/ev.fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
    td = (0:size(ev.W,1)-1)'/ev.fs - ev.pre;
    hch = {[2 3],[5 6]};
    both = intersect(Df.D(1,1).k, Df.D(2,1).k);
    snr = min(Df.NZ(1).snr(ismember(Df.D(1,1).k,both)), Df.NZ(2).snr(ismember(Df.D(2,1).k,both)));
    [~,o] = sort(snr,'descend'); both = both(o(1:min(nev,numel(o))));
    [~, ie] = ismember(both, G.ev);
    OBS = zeros(G.nt, 4, numel(ie));
    for i = 1:numel(ie)
        for g = 1:2
            x = conv2(double(ev.W(:,hch{g},both(i))), ric, 'same');
            x = interp1(td - G.off(g), x, t, 'linear', 0);    % clock correction, FD time axis
            OBS(:,2*g-1:2*g,i) = x;
        end
    end
    frameunknown = true;
else
    S = synth_records(obsname, 'forge', true);
    ie = (1:min(nev, numel(G.ev)))';
    OBS = zeros(G.nt, 4, numel(ie));
    for i = 1:numel(ie), for g = 1:2, OBS(:,2*g-1:2*g,i) = double(S(g).H(:,:,ie(i))); end, end
    frameunknown = false;
end
OBS = filtfilt(b, a, reshape(OBS, G.nt, [])); OBS = reshape(OBS, G.nt, 4, []);
% windows
win = false(G.nt, 4, numel(ie));
for i = 1:numel(ie), for g = 1:2
    tP = R(ie(i),g)/G.vp; tS = R(ie(i),g)/G.vs;
    win(:,2*g-1:2*g,i) = repmat(t >= tP-0.010 & t <= tS+0.060, 1, 2);
end, end
basis = mt_basis();
W = struct();
for m = 1:numel(models)
    L = load(fullfile(fileparts(mfilename('fullpath')),'..','sgt',[models{m} '.mat']),'E');
    GF = zeros(G.nt, 4, 6, numel(ie));         % Green's functions, FD frame (E,N)
    for i = 1:numel(ie), for g = 1:2, for n = 1:2, for j = 1:6
        GF(:,2*g-2+n,j,i) = cumsum(sgt_synth(L.E{g,n}(:,:,ie(i)), basis{j}))*G.dt;
    end, end, end, end
    GF = reshape(filtfilt(b, a, reshape(GF, G.nt, [])), G.nt, 4, 6, []);
    if frameunknown && m == 1       % tool azimuth & handedness per sensor (with first model)
        best = -inf(2,1); ROT = cell(2,1);
        for g = 1:2
            for hand = [1 -1], for phi = 0:5:355
                Rm = [cosd(phi) sind(phi); -hand*sind(phi) hand*cosd(phi)];
                vr = zeros(min(40,numel(ie)),1);
                for i = 1:numel(vr)
                    [vr(i)] = invert_one(OBS(:,:,i), GF(:,:,:,i), win(:,:,i), g, Rm, fs);
                end
                if median(vr) > best(g), best(g) = median(vr); ROT{g} = Rm; end
            end, end
            fprintf('sensor %c orientation: median single-sensor VR %.2f\n', char('A'+g-1), best(g));
        end
    elseif ~frameunknown
        ROT = {eye(2), eye(2)};
    end
    vr = zeros(numel(ie),1); rc = zeros(numel(ie),1); ro = rc;
    for i = 1:numel(ie)
        [vr(i), res, pred] = invert_one(OBS(:,:,i), GF(:,:,:,i), win(:,:,i), 0, ROT, fs);
        % unexplained energy in S coda window (S+25..S+95 ms) relative to S-window observed
        cw = false(G.nt,4); sw = cw;
        for g = 1:2
            tS = R(ie(i),g)/G.vs;
            cw(:,2*g-1:2*g) = repmat(t >= tS+0.025 & t <= tS+0.095, 1, 2);
            sw(:,2*g-1:2*g) = repmat(t >= tS-0.008 & t <= tS+0.060, 1, 2);
        end
        rc(i) = sum(res(cw).^2)/(sum(pred(sw).^2) + realmin);   % residual coda / predicted S
        ro(i) = sum(pred(cw).^2)/(sum(pred(sw).^2) + realmin);  % predicted coda / predicted S
    end
    W(m).name = models{m}; W(m).vr = vr; W(m).rcoda = rc; W(m).pcoda = ro;
    fprintf('%s vs %-3s: VR median %.3f [%.3f %.3f], predicted coda/S %.3f, unexplained coda/S %.3f\n', ...
        obsname, models{m}, median(vr), prctile(vr,25), prctile(vr,75), median(ro), median(rc));
end
out = fullfile(fileparts(mfilename('fullpath')),'..','results');
save(fullfile(out, sprintf('mtinv_%s.mat', obsname)), 'W', 'ie');
end

function [vr, res, pred] = invert_one(obs, gf, win, gsel, ROT, fs)
% linear MT inversion with per-sensor shift search (+-6 ms)
if gsel > 0, sens = gsel; Rg = {ROT, ROT}; else, sens = 1:2; Rg = ROT; end
if ~iscell(Rg), Rg = {Rg, Rg}; end
cols = []; for g = sens, cols = [cols 2*g-1 2*g]; end %#ok<AGROW>
% rotate Green's functions into tool frame
gr = gf;
for g = sens
    for j = 1:6, gr(:,2*g-1:2*g,j) = gf(:,2*g-1:2*g,j) * Rg{g}'; end
end
sh = round(-0.006*fs):round(0.001*fs):round(0.006*fs);
shift = zeros(1,2);
for it = 1:2
    A = []; d = [];
    for g = sens
        for c = 2*g-1:2*g
            w = win(:,c);
            A = [A; squeeze(circshift(gr(w,c,:),0))]; %#ok<AGROW>
            ww = circshift(w, -shift(g));
            d = [d; obs(ww,c)]; %#ok<AGROW>
        end
    end
    mt = A\d;
    % update shifts by cross-correlation of prediction with observation
    for g = sens
        best = -inf;
        for s = sh
            p = 0;
            for c = 2*g-1:2*g
                w = win(:,c); pr = squeeze(gr(w,c,:))*mt;
                p = p + obs(circshift(w,-s),c)'*pr;
            end
            if p > best, best = p; shift(g) = s; end
        end
    end
end
A = []; d = [];
for g = sens, for c = 2*g-1:2*g
    w = win(:,c); A = [A; squeeze(gr(w,c,:))]; d = [d; obs(circshift(w,-shift(g)),c)]; %#ok<AGROW>
end, end
mt = A\d; vr = 1 - sum((d - A*mt).^2)/sum(d.^2);
% full-trace residuals (for coda analysis)
pred = zeros(size(obs)); res = pred;
for g = sens, for c = 2*g-1:2*g
    pred(:,c) = squeeze(gr(:,c,:))*mt;
    res(:,c) = circshift(obs(:,c), -shift(g)) - pred(:,c);
end, end
end

function B = mt_basis()
B = {[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1], ...
     [0 1 0;1 0 0;0 0 0],[0 0 1;0 0 0;1 0 0],[0 0 0;0 0 1;0 1 0]};
end
