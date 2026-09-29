function [CL, QI, SW, BD] = coda_band_measures(H, t, tS, fs, B, fc, noiseE)
% CODA_BAND_MEASURES  Same per-band coda measures as S14, for one record.
%   H [nt x 2] horizontals; t time from origin; tS predicted S time; B [nb x 2] band edges;
%   fc band centres; noiseE (optional) pre-event noise energy per band (NaN = skip test)
nb = size(B,1); CL = NaN(1,nb); QI = CL; SW = CL; BD = CL;
for ib = 1:nb
    [b,a] = butter(4, B(ib,:)/(fs/2));
    sm = max(round(2*fs/fc(ib)), round(0.002*fs));
    X = filtfilt(b, a, double(H));
    E = movmean(sum(abs(hilbert(X)).^2, 2), sm);
    iS = find(t >= tS-0.01 & t < tS+0.03); if isempty(iS), continue; end
    [pk, j] = max(E(iS)); ts = t(iS(j));
    ic = t >= ts+0.03 & t < ts+0.10; id = t >= ts+0.02 & t < ts+0.12;
    if nnz(id) < 20 || t(end) < ts+0.12, continue; end
    if nargin > 6 && isfinite(noiseE(ib)) && mean(E(ic)) < 5*noiseE(ib), continue; end
    CL(ib) = log10(mean(E(ic))/pk);
    c1 = polyfit(t(id), log(E(id)), 1); BD(ib) = -c1(1);
    tl = t(id); c2 = polyfit(tl, log(E(id).*tl.^2), 1); QI(ib) = -c2(1)/(2*pi*fc(ib));
    iw = t >= ts-0.03 & t < ts+0.03; w = E(iw)/sum(E(iw)); tw = t(iw);
    mu = sum(w.*tw); SW(ib) = sqrt(sum(w.*(tw-mu).^2));
end
end
