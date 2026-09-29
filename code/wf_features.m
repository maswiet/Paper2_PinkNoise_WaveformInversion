function F = wf_features(H, t, tP, dSP, band, fs)
% WF_FEATURES  Scattering-sensitive waveform/envelope features of one
% event-sensor record, applied identically to FORGE data and FD synthetics.
%   H    : [nt x 2] horizontal components (rotation invariant energy used)
%   t    : time vector (s); tP: P onset (s); dSP: S-P time (s)
%   band : [f1 f2] Hz
% Returns struct: env (envelope on common lag axis rel. P), pcoda, spd,
% cdecay, clevel, sp (log10 energy ratios / s).
[b,a] = butter(3, band/(fs/2));
X = filtfilt(b, a, double(H));
E = sum(abs(hilbert(X)).^2, 2);
w = max(1, round(fs/band(1)/2));                 % half-period smoothing
E = movmean(E, w);
tS = tP + dSP;
iP  = t >= tP-0.004 & t < tP+0.015;
iPc = t >= tP+0.020 & t < tS-0.008;
iS  = t >= tS-0.008 & t < tS+0.060;
iC  = t >= tS+0.025 & t < tS+0.095;
F.pP = max(E(iP)); [F.pS, j] = max(E(iS)); ts = t(iS);
F.sp = log10(F.pS/F.pP);
F.spd = ts(j) - tS;                               % S peak delay
if nnz(iPc) > 10, F.pcoda = log10(mean(E(iPc))/F.pP); else, F.pcoda = NaN; end
F.clevel = log10(mean(E(iC))/F.pS);
c = polyfit(t(iC), log(E(iC)), 1); F.cdecay = -c(1);  % 1/s (energy)
% envelope on common lag axis relative to the S-envelope PEAK (robust to
% picking/Vp-Vs differences), normalised by S peak
lag = (-0.08:1/fs*4:0.12)';
F.envS = interp1(t - ts(j), E/F.pS, lag, 'linear', NaN);
% S-normalised coda measured from the S peak (mechanism-robust)
iC2 = t >= ts(j)+0.030 & t < ts(j)+0.100;
F.clevel2 = log10(mean(E(iC2))/F.pS);
c2 = polyfit(t(iC2), log(E(iC2)), 1); F.cdecay2 = -c2(1);
iW = t >= ts(j)-0.05 & t < ts(j)+0.05;            % S pulse rms width
w = E(iW)/sum(E(iW)); tw = t(iW); mu = sum(w.*tw); F.swidth = sqrt(sum(w.*(tw-mu).^2));
lagP = (-0.01:1/fs*4:0.05)';
F.envP = interp1(t - tP, E/F.pP, lagP, 'linear', NaN);
end
