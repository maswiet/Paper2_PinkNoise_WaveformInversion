% S21  Can borehole/tool ringing produce the observed 40-80 Hz S coda? (reviewer point 2)
% (1) Median coda/S spectral ratio (raw 4-kHz data): a resonance that rings produces a narrow
%     peak of coda/S at its frequency; volume scattering gives a smooth curve.
% (2) Coda spectral ratio sensor A / sensor B for the same events.
% (3) Null-ringing test: homogeneous-crust synthetics convolved with a damped resonance
%     h(t) = delta(t) + A exp(-pi f_r t/Q) sin(2 pi f_r t) (f_r 40-80 Hz, Q 5-40); the amplitude A
%     is set so that the resonance peak in |H(f)| does not exceed the largest spectral peak that
%     the observed S spectra show in 30-120 Hz; the resulting coda/S ratio is compared with FORGE.
maxNumCompThreads(2);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
load ../data/forge2022_events.mat; load ../data/forge_picks_sensors.mat
G = forge_setup(); fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
hch = {[2 3],[5 6]}; nw = round(0.07*fs); nfft = 4096; f = (0:nfft/2)'*fs/nfft;
[tap, ~] = dpss(nw, 2, 3);
mts = @(x) mean(abs(fft(tap .* x, nfft)).^2, 2);            % multitaper power, one component
[bh, ah] = butter(3, 5/(fs/2), 'high');
lsm = @(y) movmedian(y, 31);                                   % smooth baseline on log-frequency grid
fl = logspace(log10(15), log10(1800), 400)';
RS = cell(2,1); RC = cell(2,1); kk = cell(2,1);
for g = 1:2
    k = find(good(:,g)); k = k(ismember(k, G.ev));
    Rg = sqrt(sum((ev.xyz(k,:) - G.sens(g,:)).^2,2)); tS = tp(k,g) + Rg*(1/G.vs - 1/G.vp);
    Ps = NaN(numel(fl), numel(k)); Pc = Ps;
    for i = 1:numel(k)
        X = filtfilt(bh, ah, double(ev.W(:,hch{g},k(i))));
        is = find(t >= tS(i)-0.010, 1); ic = find(t >= tS(i)+0.040, 1); in = find(t >= -0.095, 1);
        if isempty(ic) || ic + nw - 1 > numel(t), continue; end
        ps = 0; pc = 0; pn = 0;
        for c = 1:2
            ps = ps + mts(X(is:is+nw-1,c)); pc = pc + mts(X(ic:ic+nw-1,c)); pn = pn + mts(X(in:in+nw-1,c));
        end
        ps = interp1(f, ps(1:numel(f)), fl); pc = interp1(f, pc(1:numel(f)), fl); pn = interp1(f, pn(1:numel(f)), fl);
        ok = pc > 3*pn;                                          % coda above noise
        Ps(:,i) = ps; pcr = pc ./ ps; pcr(~ok) = NaN; Pc(:,i) = pcr;
    end
    RS{g} = Ps; RC{g} = Pc; kk{g} = k;
end
% ---- (1) peak prominence of median spectra ----
rows = {};
for g = 1:2
    ls = log10(median(RS{g} ./ median(RS{g}(fl>30 & fl<150,:),1), 2, 'omitnan'));
    lc = log10(median(RC{g}, 2, 'omitnan'));
    ps = ls - lsm(ls); pc = lc - lsm(lc);
    for bnd = {[30 120], [280 380], [420 560], [850 1000]}
        m = fl >= bnd{1}(1) & fl <= bnd{1}(2);
        rows(end+1,:) = {char('A'+g-1), sprintf('%d-%d', bnd{1}), 10*max(ps(m)), 10*max(pc(m))}; %#ok<SAGROW>
    end
    SP{g} = ls; CP{g} = lc; %#ok<SAGROW>
end
T1 = cell2table(rows, 'VariableNames', {'sensor','band_Hz','S_spectrum_peak_prominence_dB','coda_over_S_peak_prominence_dB'});
disp(T1); writetable(T1, fullfile(out,'table_ringing_prominence.csv'));
% ---- (2) A/B coda spectral ratio for common events ----
[cm, ia, ib] = intersect(kk{1}, kk{2});
rab = log10(median(RC{1}(:,ia) ./ RC{2}(:,ib), 2, 'omitnan'));
m = fl > 30 & fl < 150; pr_ab = 10*max(rab(m) - lsm(rab(m)));
fprintf('A/B coda/S ratio: median over 30-150 Hz %.2f (log10), max prominence %.1f dB\n', median(rab(m),'omitnan'), pr_ab);
% ---- (3) null-ringing test with homogeneous synthetics ----
S = synth_records('H', 'forge', false); fs2 = 1/G.dt; t2 = S(1).t;
[b,a] = butter(3,[40 80]/(fs2/2)); env = @(x) movmean(sum(abs(hilbert(x)).^2,2), 20);
Dq = load(fullfile(out,'s19_H.mat'));                     % observed coda/S at the primary windows
obsA = median(Dq.obsratio(:,1,1),'omitnan'); obsB = median(Dq.obsratio(:,2,1),'omitnan');
pmax = max(T1.S_spectrum_peak_prominence_dB(strcmp(T1.band_Hz,'30-120')));   % dB allowed
th = (0:round(0.3*fs2))'/fs2; rows = {};
for fr = 40:10:80
    for Q = [5 10 20 40]
        r = exp(-pi*fr*th/Q) .* sin(2*pi*fr*th);
        % scale A so that the peak of 20log10|1+R(f)| equals pmax dB
        Rf = fft(r, 8192)/fs2; ff = (0:4095)'*fs2/8192;
        Afun = @(A) max(20*log10(abs(1 + A*Rf(ff>20 & ff<150)))) - pmax;
        A = fzero(Afun, [0 1e4]);
        h = [1; zeros(numel(th)-1,1)] + A*r/fs2;             % discrete: delta + A r dt
        lr = NaN(0,2);
        for g = 1:2
            ne = size(S(g).H,3); v = NaN(ne,1);
            for i = 1:ne
                x = filter(h, 1, double(S(g).H(:,:,i)));
                x = filtfilt(b, a, x); e = env(x);
                tS = S(g).R(i)/G.vs;
                iS = t2 >= tS-0.008 & t2 <= tS+0.015; iC = t2 >= tS+0.040 & t2 < tS+0.100;
                if t2(end) < tS+0.100, continue; end
                v(i) = mean(e(iC))/max(e(iS));
            end
            lr(g) = log10(median(v,'omitnan')); %#ok<SAGROW>
        end
        rows(end+1,:) = {fr, Q, A, pmax, lr(1), lr(2), lr(1) - log10(obsA), lr(2) - log10(obsB)}; %#ok<SAGROW>
    end
end
% baseline without ringing
T3 = cell2table(rows, 'VariableNames', {'f_r_Hz','Q','A','peak_dB','log10_coda_S_A','log10_coda_S_B','deficit_A','deficit_B'});
disp(T3); writetable(T3, fullfile(out,'table_null_ringing.csv'));
fprintf('observed median log10(coda/S): A %.2f, B %.2f\n', log10(obsA), log10(obsB));
% ---- figure ----
f1 = figure('Visible','off');
for g = 1:2
    subplot(1,3,g); semilogx(fl, 10*SP{g}, 'k', fl, 10*CP{g} - median(10*CP{g}(fl>30 & fl<150),'omitnan'), 'r', 'LineWidth', 1.2);
    xline([40 80], ':'); xline([330 480 920], 'Color', [0.6 0.6 0.6]); grid on; xlim([15 1800]);
    set(gca, 'XTick', [20 50 100 200 500 1000], 'XTickLabel', {'20','50','100','200','500','1000'});
    xlabel('frequency (Hz)'); if g == 1, ylabel('relative level (dB)'); end
    lg = legend('S spectrum', 'coda/S ratio', 'Location','southwest'); lg.Box = 'off';
    panel_label(gca, g, 'tl'); set(gca, 'Position', [0.07 + 0.315*(g-1) 0.17 0.27 0.78]);
end
subplot(1,3,3); hold on
for Q = [5 10 20 40]
    s3 = T3(T3.Q == Q,:); plot(s3.f_r_Hz, s3.log10_coda_S_A, '-o', 'MarkerSize', 3, 'DisplayName', sprintf('H + resonance, Q = %d', Q));
end
yline(log10(obsA), 'k-', 'LineWidth', 1.5, 'DisplayName', 'FORGE, sensor A'); grid on; box on
xlabel('resonance frequency f_r (Hz)'); ylabel('log_{10}(coda/S)'); lg = legend('Location','east'); lg.Box = 'off';
panel_label(gca, 3, 'tl'); set(gca, 'Position', [0.74 0.17 0.24 0.78]);
pub_export(f1, fullfile(here,'..','figs','fig_ringing_test.png'), 17.4, 6.5, 7); close(f1);
