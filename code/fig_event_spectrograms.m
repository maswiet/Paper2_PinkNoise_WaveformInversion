% Example FORGE 2022 MEQs at the 56-32 downhole sensors (raw 4 kHz records):
% waveforms, spectrograms and signal/noise spectra, plus population bandwidth
% statistics. Justifies analysis above the ~100 Hz Nyquist of 200-sps networks.
maxNumCompThreads(2);
load ../data/forge2022_events.mat; load ../data/forge_picks_sensors.mat; Df = load('../data/forge_features.mat');
G = forge_setup(); fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
hch = {[2 3],[5 6]};
[bh,ah] = butter(3, 5/(fs/2), 'high');
nw = round(0.08*fs); nfft = 2048; f = (0:nfft/2)'*fs/nfft;
hw = hann(nw);
spec = @(x) abs(fft(x.*hw, nfft)); % per component
% ---------------- population statistics ----------------
SNRf = cell(2,1); fr100 = cell(2,1); fpk = cell(2,1); fmax3 = cell(2,1); Ms = cell(2,1);
for g = 1:2
    k = find(good(:,g)); k = k(ismember(k, G.ev));
    R = sqrt(sum((ev.xyz(k,:) - G.sens(g,:)).^2,2)); tS = tp(k,g) + R*(1/G.vs - 1/G.vp);
    S = zeros(numel(f), numel(k)); N = S; ok = false(numel(k),1);
    for i = 1:numel(k)
        X = filtfilt(bh, ah, double(ev.W(:,hch{g},k(i))));
        is = find(t >= tS(i) - 0.01, 1); in = find(t >= -0.095, 1);
        if isempty(is) || is + nw - 1 > numel(t), continue; end
        xs = X(is:is+nw-1,:); xn = X(in:in+nw-1,:);
        s2 = spec(xs(:,1)).^2 + spec(xs(:,2)).^2; n2 = spec(xn(:,1)).^2 + spec(xn(:,2)).^2;
        S(:,i) = s2(1:numel(f)); N(:,i) = n2(1:numel(f)); ok(i) = true;
    end
    S = S(:,ok); N = N(:,ok); Ms{g} = ev.M(k(ok));
    SNRf{g} = sqrt(movmean(S,9) ./ movmean(N,9));          % amplitude SNR per frequency
    band = f >= 5 & f <= 1900;
    for i = 1:size(S,2)
        s = SNRf{g}(:,i); j = find(band & s >= 3, 1, 'last'); fmax3{g}(i) = f(j);
        e = S(:,i) - N(:,i); e(e < 0) = 0;
        fr100{g}(i) = sum(e(f > 100 & band)) / sum(e(band));
        [~, jp] = max(movmean(S(:,i),9) .* band); fpk{g}(i) = f(jp);
    end
end
rows = {};
for g = 1:2
    q = @(v) prctile(v, [25 50 75]);
    a = q(fmax3{g}); b = q(100*fr100{g}); c = q(fpk{g});
    rows(end+1,:) = {char('A'+g-1), numel(fpk{g}), a(2), a(1), a(3), b(2), b(1), b(3), c(2), c(1), c(3), ...
        100*mean(fmax3{g} > 100), 100*mean(fmax3{g} > 500)}; %#ok<SAGROW>
end
T = cell2table(rows, 'VariableNames', {'sensor','n_events','fmax_snr3_med','fmax_q25','fmax_q75', ...
    'pct_S_energy_above100_med','pct_q25','pct_q75','fpeak_med','fpeak_q25','fpeak_q75', ...
    'pct_events_fmax_gt100','pct_events_fmax_gt500'});
disp(T); writetable(T, '../results/table_bandwidth.csv');

% ---------------- example events ----------------
both = intersect(find(good(:,1)), find(good(:,2))); both = both(ismember(both, G.ev));
% keep events whose P picks agree with the calibrated travel-time model (|residual| < 2 ms)
res_ok = true(size(both));
for g = 1:2
    Rg = sqrt(sum((ev.xyz(both,:) - G.sens(g,:)).^2,2));
    r = tp(both,g) - (Rg/G.vp + G.off(g)); r = r - median(r);
    res_ok = res_ok & abs(r) < 0.002;
end
both = both(res_ok);
sn = min(q_snr(ev, both, tp, 1, hch, t, bh, ah), q_snr(ev, both, tp, 2, hch, t, bh, ah));
[~, o] = sort(sn, 'descend'); cand = both(o(1:40));
[~, om] = sort(ev.M(cand)); pick = cand(om(round(linspace(1, numel(om), 4))));   % span magnitudes
f1 = figure('Position',[20 20 1650 1250], 'Visible','off');
for r = 1:4
    k = pick(r); g = 2;                                     % sensor B (closest)
    R = norm(ev.xyz(k,:) - G.sens(g,:)); tS = tp(k,g) + R*(1/G.vs - 1/G.vp);
    X = filtfilt(bh, ah, double(ev.W(:,hch{g},k)));
    t1 = tp(k,g) - 0.03; t2 = tS + 0.15;                  % common time window for the row
    w = t > t1 & t < t2;
    % waveform
    subplot(4,3,3*r-2); plot(t(w), X(w,1)/max(abs(X(w,1))), 'k', 'LineWidth', 0.4); hold on
    xline(tp(k,g),'b-','P','LineWidth',1); xline(tS,'r-','S','LineWidth',1); xlim([t1 t2]); ylim([-1.1 1.1]);
    ylabel(sprintf('M_w %.2f, R = %.0f m', ev.M(k), R)); if r == 4, xlabel('time after origin (s)'); end
    panel_label(gca, 3*r-2, 'bl');
    box on
    % spectrogram (same dB scale and colour bar in every row)
    subplot(4,3,3*r-1);
    i0 = find(t >= t1 - 0.012, 1); i1 = find(t <= t2 + 0.012, 1, 'last');
    [s, fq, tq] = spectrogram(X(i0:i1, 1), hann(96), 88, 512, fs);
    imagesc(tq + t(i0), fq, 20*log10(abs(s)/max(abs(s(:))))); axis xy; caxis([-60 0]);
    xlim([t1 t2]); ylim([0 2000]); hold on
    yline(100,'w--','LineWidth',1); yline(50,'w:','LineWidth',1);
    if r == 1
        text(t2 - 0.005, 180, '200 sps Nyquist', 'Color','w', 'HorizontalAlignment','right', 'FontSize',7);
        text(t2 - 0.005, 20, '100 sps', 'Color','w', 'HorizontalAlignment','right', 'FontSize',7);
    end
    xline(tp(k,g),'w-','LineWidth',1); xline(tS,'w-','LineWidth',1);
    colormap(gca, turbo); ylabel('frequency (Hz)');
    cb = colorbar; cb.Label.String = 'dB re max'; cb.Ticks = -60:20:0;
    panel_label(gca, 3*r-1, 'tl');
    if r == 4, xlabel('time after origin (s)'); end
    % spectra
    subplot(4,3,3*r);
    is = find(t >= tS - 0.01, 1); in = find(t >= -0.095, 1);
    xs = X(is:is+nw-1,:); xn = X(in:in+nw-1,:);
    ss = sqrt(spec(xs(:,1)).^2 + spec(xs(:,2)).^2); nn = sqrt(spec(xn(:,1)).^2 + spec(xn(:,2)).^2);
    ss = movmean(ss(1:numel(f)),5); nn = movmean(nn(1:numel(f)),5);
    loglog(f(2:end), ss(2:end)/max(ss), 'k', 'LineWidth', 1); hold on
    loglog(f(2:end), nn(2:end)/max(ss), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    xline([50 100], 'k:'); patch([20 160 160 20], [1e-4 1e-4 3 3], [0.85 0.33 0.1], 'FaceAlpha', 0.1, 'EdgeColor','none');
    xlim([5 2000]); ylim([1e-4 3]); grid on;
    if r == 1, legend('signal','noise','Location','southwest'); end
    panel_label(gca, 3*r, 'tl'); ylabel('normalised amplitude');
    if r == 4, xlabel('frequency (Hz)'); end
end
exportgraphics(f1, '../figs/fig_event_spectrograms.png', 'Resolution', 140); close(f1);

% ---------------- population figure ----------------
f2 = figure('Position',[20 20 1300 450], 'Visible','off');
for g = 1:2
    subplot(1,3,g);
    Q = prctile(SNRf{g}, [25 50 75], 2);
    loglog(f(2:end), Q(2:end,2), 'k', 'LineWidth', 2); hold on
    loglog(f(2:end), Q(2:end,1), 'k:', f(2:end), Q(2:end,3), 'k:');
    yline(3, 'r--', 'SNR 3'); xline([50 100], 'b:'); xlim([5 2000]); grid on
    xlabel('frequency (Hz)'); ylabel('amplitude SNR (S window / noise)');
    panel_label(gca, g, 'tl');
end
subplot(1,3,3);
histogram(fpk{1}, logspace(1, log10(2000), 30), 'Normalization','probability'); hold on
histogram(fpk{2}, logspace(1, log10(2000), 30), 'Normalization','probability');
set(gca,'XScale','log'); xline([50 100], 'k--', {'100 sps','200 sps'}); xlim([10 2000]);
legend('sensor A','sensor B','Location','northwest'); xlabel('peak frequency of S-window velocity spectrum (Hz)');
ylabel('fraction of events');
panel_label(gca, 3, 'tr');
exportgraphics(f2, '../figs/fig_bandwidth.png', 'Resolution', 140); close(f2);
disp(ev.M(pick)');

function sn = q_snr(ev, k, tp, g, hch, t, bh, ah)
sn = zeros(numel(k),1);
for i = 1:numel(k)
    X = filtfilt(bh, ah, double(ev.W(:,hch{g},k(i))));
    sn(i) = max(abs(X(t > tp(k(i),g) & t < tp(k(i),g)+0.2, 1))) / (std(X(t < -0.01,1)) + eps);
end
end
