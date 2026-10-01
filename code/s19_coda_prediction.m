function s19_coda_prediction(models, nthreads)
% S19  Per-event S-coda prediction with NON-OVERLAPPING fit and coda windows
% (reviewer points 1, 7, 8).  For every event recorded by both sensors (not only the 150
% highest-SNR events), the moment tensor is fitted to horizontal 40-80 Hz traces in
% [P-10 ms, S+fitEnd] and the coda energy is evaluated in [S+c1, S+c2]:
%   fitEnd = 15 ms (primary), 10, 20 ms (sensitivity), 60 ms (original, overlapping; for comparison)
%   coda   = 40-100 ms (primary), 50-110 ms, 30-100 ms (original)
% At the primary windows the full, deviatoric and double-couple-constrained MTs are all used.
% Output: ../results/s19_<model>.mat
if nargin > 1, maxNumCompThreads(nthreads); end
here = fileparts(mfilename('fullpath'));
G = forge_setup(); fs = 1/G.dt; t = (0:G.nt-1)'*G.dt - G.t0;
[b,a] = butter(3,[40 80]/(fs/2));
fitEnd = [0.015 0.010 0.020 0.060];
coda   = [0.040 0.100; 0.050 0.110; 0.030 0.100];
% ---------------- observed ----------------
ev = load(fullfile(here,'..','data','forge2022_events.mat')); ev = ev.ev;
Df = load(fullfile(here,'..','data','forge_features.mat'));
tr = (-0.05:1/ev.fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
td = (0:size(ev.W,1)-1)'/ev.fs - ev.pre; hch = {[2 3],[5 6]};
both = intersect(Df.D(1,1).k, Df.D(2,1).k);
[~, ie] = ismember(both, G.ev); keep = ie > 0; both = both(keep); ie = ie(keep);
snr = min(Df.NZ(1).snr(ismember(Df.D(1,1).k,both)), Df.NZ(2).snr(ismember(Df.D(2,1).k,both)));
[~,o] = sort(snr,'descend'); top150 = false(numel(both),1); top150(o(1:min(150,end))) = true;
ne = numel(ie);
OBS = zeros(G.nt, 4, ne);
for i = 1:ne, for g = 1:2
    x = conv2(double(ev.W(:,hch{g},both(i))), ric, 'same');
    OBS(:,2*g-1:2*g,i) = interp1(td - G.off(g), x, t, 'linear', 0);
end, end
OBS = reshape(filtfilt(b, a, reshape(OBS, G.nt, [])), G.nt, 4, []);
R = [sqrt(sum((G.evxyz(ie,:) - G.sens(1,:)).^2,2)) sqrt(sum((G.evxyz(ie,:) - G.sens(2,:)).^2,2))];
tP = R/G.vp; tS = R/G.vs;
L = load(fullfile(here,'..','data','sensor_rotation.mat')); ROT = L.ROT;
% unit double couples in basis order [xx yy zz xy xz yz]
[S1, D1, K1] = ndgrid(0:10:350, 10:10:90, -180:15:165);
DCg = zeros(6, numel(S1));
for k = 1:numel(S1)
    Mk = dc_moment(S1(k), D1(k), K1(k)); DCg(:,k) = [Mk(1,1) Mk(2,2) Mk(3,3) Mk(1,2) Mk(1,3) Mk(2,3)]';
end
basis = mt_basis();
env = @(x) movmean(sum(abs(hilbert(x)).^2, 2), 20);
for m = 1:numel(models)
    Lm = load(fullfile(here,'..','sgt',[models{m} '.mat']),'E');
    GF = zeros(G.nt, 4, 6, ne);
    for i = 1:ne, for g = 1:2, for n = 1:2, for j = 1:6
        GF(:,2*g-2+n,j,i) = cumsum(sgt_synth(Lm.E{g,n}(:,:,ie(i)), basis{j}))*G.dt;
    end, end, end, end
    GF = reshape(filtfilt(b, a, reshape(GF, G.nt, [])), G.nt, 4, 6, []);
    lr = NaN(ne, 2, numel(fitEnd), size(coda,1)); vr = NaN(ne, numel(fitEnd));
    lrv = NaN(ne, 2, 3); vrv = NaN(ne, 3);                 % full / dev / dc at primary windows
    cnd = NaN(ne,1); dec = NaN(ne,3); obsr = NaN(ne, 2, size(coda,1));
    for i = 1:ne
        for f = 1:numel(fitEnd)
            win = false(G.nt, 4);
            for g = 1:2, win(:,2*g-1:2*g) = repmat(t >= tP(i,g)-0.010 & t <= tS(i,g)+fitEnd(f), 1, 2); end
            if f == 1, vars = {'full','dev','dc'}; else, vars = {'full'}; end
            Rf = mt_fit_variants(OBS(:,:,i), GF(:,:,:,i), win, ROT, fs, vars, DCg);
            vr(i,f) = Rf.full.vr;
            if f == 1, cnd(i) = Rf.cond; dec(i,:) = mt_decomp(Rf.full.m); end
            for v = 1:numel(vars)
                r = Rf.(vars{v});
                for g = 1:2
                    c = 2*g-1:2*g;
                    eo = env(r.obsS(:,c)); ep = env(r.pred(:,c));
                    iS = t >= tS(i,g)-0.008 & t <= tS(i,g)+0.015;
                    for q = 1:size(coda,1)
                        iC = t >= tS(i,g)+coda(q,1) & t < tS(i,g)+coda(q,2);
                        if t(end) < tS(i,g)+coda(q,2), continue; end
                        ro = mean(eo(iC))/max(eo(iS)); rp = mean(ep(iC))/max(ep(iS));
                        if v == 1, lr(i,g,f,q) = log10(rp/ro); if f == 1, obsr(i,g,q) = ro; end, end
                        if f == 1 && q == 1, lrv(i,g,v) = log10(rp/ro); vrv(i,v) = r.vr; end
                    end
                end
            end
        end
    end
    S = struct('model', models{m}, 'ev', G.ev(ie), 'R', R, 'top150', top150, 'fitEnd', fitEnd, 'coda', coda, ...
        'lr', lr, 'vr', vr, 'lrv', lrv, 'vrv', vrv, 'cond', cnd, 'decomp', dec, 'obsratio', obsr);
    save(fullfile(here,'..','results',sprintf('s19_%s.mat', models{m})), '-struct', 'S');
    fprintf('%-6s n=%d | primary (fit S+15, coda S+40-100): median log10(pred/obs) A %+.2f B %+.2f | within x3 A %.0f%% B %.0f%% | VR %.2f | old windows A %+.2f\n', ...
        models{m}, ne, median(lr(:,1,1,1),'omitnan'), median(lr(:,2,1,1),'omitnan'), ...
        100*mean(abs(lr(:,1,1,1)) < log10(3),'omitnan'), 100*mean(abs(lr(:,2,1,1)) < log10(3),'omitnan'), ...
        median(vr(:,1),'omitnan'), median(lr(:,1,4,3),'omitnan'));
end
end

function d = mt_decomp(m)
% percentages [ISO CLVD DC] (Vavrycuk 2001 style)
M = [m(1) m(4) m(5); m(4) m(2) m(6); m(5) m(6) m(3)];
iso = trace(M)/3; e = eig(M - iso*eye(3)); [~, o] = sort(abs(e), 'descend'); e = e(o);
epsl = -e(3)/abs(e(1));
ci = iso/(abs(iso) + abs(e(1))); cc = 2*epsl*(1 - abs(ci)); cd = 1 - abs(ci) - abs(cc);
d = 100*[ci cc cd];
end
