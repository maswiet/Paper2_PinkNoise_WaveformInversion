% S06  Waveform features of the FORGE 2022 records (horizontal pairs),
% plus noise library for synthetic contamination at matched SNR.
load ../data/forge2022_events.mat; load ../data/tt_inversion.mat
G = forge_setup();
fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
bands = [20 40; 40 80];
hch = {[2 3],[5 6]};
if exist('../data/forge_features.mat','file'), Dold = load('../data/forge_features.mat'); Dold = Dold.D; if ~isfield(Dold,'tpl'), clear Dold; end, end
% convolve data with the FD source wavelet (Ricker, G.f0) so that data and
% synthetics share the same source spectrum (real MEQ corner >> 80 Hz)
tr = (-0.05:1/fs:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
W = zeros(size(ev.W),'single');
for i = 1:size(ev.W,3), W(:,:,i) = single(conv2(double(ev.W(:,:,i)), ric, 'same')); end
ev.W = W; clear W
D = struct(); NZ = struct();
for g = 1:2
    sel = obs(:,2)==g & ismember(obs(:,1), G.ev);
    k = obs(sel,1); tp = obs(sel,3);
    if exist('Dold','var'), tp = Dold(g,1).tpl; end   % low-band AIC picks (same picker as synthetics)
    R = sqrt(sum((ev.xyz(k,:) - G.sens(g,:)).^2,2));
    dSP = R*(1/G.vs - 1/G.vp);
    for ib = 1:2
        for i = 1:numel(k)
            H = double(ev.W(:,hch{g},k(i)));
            F = wf_features(H, t, tp(i), dSP(i), bands(ib,:), fs);
            D(g,ib).F(i) = F;
        end
        D(g,ib).k = k; D(g,ib).R = R; D(g,ib).tp = tp; D(g,ib).dSP = dSP;
    end
    % noise library: late-window samples (0.60-0.95 s) and P-energy scale (20-80 Hz)
    [b,a] = butter(3,[20 80]/(fs/2));
    nw = t >= 0.60 & t < 0.60 + G.nt*G.dt;
    N = zeros(nnz(nw)/(G.dt*fs)*0+ceil(nnz(nw)/(G.dt*fs)), 2, numel(k), 'single');
    snr = zeros(numel(k),1);
    for i = 1:numel(k)
        H = double(ev.W(:,hch{g},k(i)));
        X = filtfilt(b,a,H);
        pe = max(sum(X(t>=tp(i)-0.004 & t<tp(i)+0.015,:).^2,2));
        ne = mean(sum(X(nw,:).^2,2));
        snr(i) = pe/ne;
        Hn = H(nw,:);                                % resample to FD dt (1/4000 = 2.5e-4 already)
        N(:,:,i) = single(Hn(1:size(N,1),:) / sqrt(pe));  % noise normalised by P peak energy
    end
    NZ(g).N = N; NZ(g).snr = snr;
end
for g = 1:2, for ib = 1:2, D(g,ib).tpl = D(g,ib).tp; end, end
save ../data/forge_features.mat D NZ bands
% ---- summary ----
nm = {'pcoda','spd','cdecay','clevel','sp','clevel2','cdecay2','swidth'};
for g = 1:2, for ib = 1:2
    fprintf('sensor %c band %d-%d Hz (n=%d):', char('A'+g-1), bands(ib,:), numel(D(g,ib).F));
    for j = 1:numel(nm)
        v = [D(g,ib).F.(nm{j})];
        fprintf('  %s %.3g [%.3g %.3g]', nm{j}, median(v,'omitnan'), prctile(v,25), prctile(v,75));
    end
    fprintf('\n');
end, end
