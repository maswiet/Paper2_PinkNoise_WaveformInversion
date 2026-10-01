% S26  Numerical convergence of the coda measures (reviewer point 12): compare C25 (2.5-m grid),
% CBIG (domain +150 m) and CSP (200-m sponge) with the 5-m reference C5 for the same medium,
% sensor (B), events, source (Ricker 45 Hz) and double-couple mechanisms.
% Measures (horizontal energy, 20-40 and 40-80 Hz): coda/S ratio S+40..100 ms (primary window),
% coda level 30-100 ms after the S peak, energy decay rate, S-pulse width, Qc^-1.
here = fileparts(mfilename('fullpath')); cdir = fullfile(here,'..','conv'); out = fullfile(here,'..','results');
cases = {'C5','C25','CBIG','CSP'}; basisM = [];
B = [20 40; 40 80]; fc = [28 57];
rows = {}; V = struct();
for k = 1:numel(cases)
    L = load(fullfile(cdir, [cases{k} '.mat'])); C = L.C; E = L.E;
    ne = numel(C.ev); nt = size(E{2,1},1); t = (0:nt-1)'*C.dt - C.t0; fs = 1/C.dt;
    rng(1234); str = 50*rand(ne,1); dip = 55 + 30*rand(ne,1); rk = -90 + 30*randn(ne,1);
    R = sqrt(sum((C.evxyz - C.sens(2,:)).^2, 2)); tS = R/C.vs;
    cr = NaN(ne,2); cl = cr; qi = cr; sw = cr; bd = cr;
    for i = 1:ne
        Mt = dc_moment(str(i), dip(i), rk(i)); H = zeros(nt,2);
        for n = 1:2, H(:,n) = cumsum(sgt_synth(E{2,n}(:,:,i), Mt))*C.dt; end
        [cl(i,:), qi(i,:), sw(i,:), bd(i,:)] = coda_band_measures(H, t, tS(i), fs, B, fc);
        for ib = 1:2
            [b,a] = butter(3, B(ib,:)/(fs/2)); e = movmean(sum(abs(hilbert(filtfilt(b,a,H))).^2,2), 20);
            iS = t >= tS(i)-0.008 & t <= tS(i)+0.015; iC = t >= tS(i)+0.040 & t < tS(i)+0.100;
            if t(end) >= tS(i)+0.100, cr(i,ib) = log10(mean(e(iC))/max(e(iS))); end
        end
    end
    V.(cases{k}) = struct('cr', cr, 'cl', cl, 'qi', log10(qi), 'sw', 1e3*sw, 'bd', bd, 'ev', C.ev);
end
fld = {'cr','cl','qi','sw','bd'}; lab = {'log10 coda/S (S+40-100 ms)','coda level (log10)','log10 Qc^-1','S width (ms)','decay rate (1/s)'};
for k = 2:numel(cases)
    for j = 1:numel(fld)
        for ib = 1:2
            a = double(V.C5.(fld{j})(:,ib)); b = double(V.(cases{k}).(fld{j})(:,ib)); ok = isfinite(a) & isfinite(b);
            dd = b(ok) - a(ok);                 % explicit MAD (R2026b mad -> nanmedian -> prctile failed on these inputs)
            rows(end+1,:) = {cases{k}, lab{j}, sprintf('%d-%d Hz', B(ib,:)), median(a(ok)), median(b(ok)), ...
                median(dd), 1.4826*median(abs(dd - median(dd))), nnz(ok)}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'case','measure','band','median_C5','median_case','median_diff','robust_std_diff','n'});
disp(T); writetable(T, fullfile(out, 'table_convergence.csv'));
