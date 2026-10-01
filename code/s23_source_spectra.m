% S23  Source-spectrum evidence and source-duration sensitivity (reviewer point 7).
% (1) Apparent corner of the S-window displacement spectrum (raw 4-kHz data, both horizontals):
%     first frequency above 50 Hz where the smoothed spectrum falls 6 dB below its 30-50 Hz level
%     (searched to 300 Hz, below the suspected resonances; includes path attenuation, so it is a
%     LOWER bound on the source corner); also the 30-150 Hz displacement-spectrum slope.
% (2) Sensitivity: homogeneous and pink (sigma .13) synthetics convolved with a Brune moment-rate
%     function of corner fc = 60, 100, 150 Hz; change of the 40-80 Hz coda/S ratio and S width.
maxNumCompThreads(2);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
load ../data/forge2022_events.mat; load ../data/forge_picks_sensors.mat
G = forge_setup(); fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
hch = {[2 3],[5 6]}; nw = round(0.07*fs); nfft = 4096; f = (0:nfft/2)'*fs/nfft;
tap = dpss(nw, 2, 3); [bh, ah] = butter(3, 5/(fs/2), 'high');
rows = {};
for g = 1:2
    k = find(good(:,g)); k = k(ismember(k, G.ev));
    Rg = sqrt(sum((ev.xyz(k,:) - G.sens(g,:)).^2,2)); tS = tp(k,g) + Rg*(1/G.vs - 1/G.vp);
    for i = 1:numel(k)
        X = filtfilt(bh, ah, double(ev.W(:,hch{g},k(i))));
        is = find(t >= tS(i)-0.010, 1); in = find(t >= -0.095, 1);
        p = 0; pn = 0;
        for c = 1:2
            p = p + mean(abs(fft(tap .* X(is:is+nw-1,c), nfft)).^2, 2);
            pn = pn + mean(abs(fft(tap .* X(in:in+nw-1,c), nfft)).^2, 2);
        end
        p = p(1:numel(f)); pn = pn(1:numel(f));
        d = sqrt(p) ./ (2*pi*max(f,1));                      % displacement amplitude
        ds = movmean(d, 9); ok = p > 9*pn;                    % amplitude SNR > 3
        L0 = mean(ds(f >= 30 & f <= 50));
        j = find(f > 50 & f <= 300 & ds < 0.5*L0, 1);
        fca = NaN; if ~isempty(j), fca = f(j); end
        m = f >= 30 & f <= 150 & ok;
        sl = NaN; if nnz(m) > 10, c1 = polyfit(log10(f(m)), log10(d(m)), 1); sl = c1(1); end
        rows(end+1,:) = {char('A'+g-1), k(i), ev.M(k(i)), Rg(i), fca, sl}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'sensor','event','Mw','R_m','fc_apparent_Hz','disp_slope_30_150'});
writetable(T, fullfile(out,'table_source_spectra.csv'));
for s = 'AB'
    x = T(strcmp(T.sensor, s), :);
    fprintf('sensor %c: n=%d | apparent corner <= 300 Hz found in %.0f%% | of those median %.0f Hz | fraction with corner < 80 Hz %.1f%% | disp slope 30-150 Hz median %.2f [%.2f %.2f]\n', ...
        s, height(x), 100*mean(isfinite(x.fc_apparent_Hz)), median(x.fc_apparent_Hz,'omitnan'), ...
        100*mean(x.fc_apparent_Hz < 80), median(x.disp_slope_30_150,'omitnan'), prctile(x.disp_slope_30_150,25), prctile(x.disp_slope_30_150,75));
    for mb = [-1 -0.3; -0.3 0.3; 0.3 1.5]'
        y = x(x.Mw >= mb(1) & x.Mw < mb(2), :);
        fprintf('    Mw %.1f..%.1f: n=%d, corner<80 Hz %.1f%%, disp slope median %.2f\n', mb, height(y), 100*mean(y.fc_apparent_Hz < 80), median(y.disp_slope_30_150,'omitnan'));
    end
end
% ---- (2) source-duration sensitivity on synthetics ----
fs2 = 1/G.dt; [b,a] = butter(3,[40 80]/(fs2/2)); env = @(x) movmean(sum(abs(hilbert(x)).^2,2), 20);
tt = (0:round(0.1*fs2))'/fs2; rows = {};
for mdl = {'H','P6'}
    S = synth_records(mdl{1}, 'forge', false); t2 = S(1).t;
    for fc = [Inf 150 100 60]
        if isinf(fc), stf = 1; else, w = 2*pi*fc; stf = w^2*tt.*exp(-w*tt)/fs2; end   % Brune moment rate, unit area
        for g = 1:2
            ne = size(S(g).H,3); v = NaN(ne,1); sw = v;
            for i = 1:ne
                x = filtfilt(b, a, filter(stf, 1, double(S(g).H(:,:,i)))); e = env(x);
                tS = S(g).R(i)/G.vs; iS = t2 >= tS-0.008 & t2 <= tS+0.015; iC = t2 >= tS+0.040 & t2 < tS+0.100;
                if t2(end) < tS+0.100, continue; end
                [~, jp] = max(e .* iS); tp_ = t2(jp); iw = t2 >= tp_-0.03 & t2 < tp_+0.03; ww = e(iw)/sum(e(iw));
                v(i) = mean(e(iC))/max(e(iS)); mu = sum(ww.*t2(iw)); sw(i) = sqrt(sum(ww.*(t2(iw)-mu).^2));
            end
            rows(end+1,:) = {mdl{1}, fc, char('A'+g-1), log10(median(v,'omitnan')), 1e3*median(sw,'omitnan')}; %#ok<SAGROW>
        end
    end
end
T2 = cell2table(rows, 'VariableNames', {'model','brune_fc_Hz','sensor','log10_coda_over_S','S_width_ms'});
disp(T2); writetable(T2, fullfile(out,'table_source_duration_sensitivity.csv'));
