% S14  Multi-octave S-coda analysis of the FORGE 2022 56-32 records (raw 4 kHz).
% Per event, sensor and octave band: coda level (E 30-100 ms after S peak / S peak),
% energy decay rate b, single-scattering Qc^-1 (Aki & Chouet: E t^2 ~ exp(-2 pi f t / Qc)),
% and rms S-pulse width.  All are within-band energy ratios, so a linear instrument /
% coupling response (including resonances) cancels.  Test: single power law in f vs a
% two-segment law with a break (BIC) - a correlation-length medium predicts a break
% near k a ~ 1, a scale-free medium none.
maxNumCompThreads(2);
load ../data/forge2022_events.mat; load ../data/forge_picks_sensors.mat
G = forge_setup(); fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
hch = {[2 3],[5 6]};
fc = [25 50 100 200 400 800]; B = [fc'/sqrt(2) fc'*sqrt(2); 1000 1700]; fc = [fc 1304];
nb = size(B,1);
res = struct();
for g = 1:2
    k = find(good(:,g)); k = k(ismember(k, G.ev)); ne = numel(k);
    R = sqrt(sum((ev.xyz(k,:) - G.sens(g,:)).^2,2));
    tS = tp(k,g) + R*(1/G.vs - 1/G.vp);
    CL = NaN(ne,nb); BD = CL; QI = CL; SW = CL; OK = false(ne,nb);
    for ib = 1:nb
        [b,a] = butter(4, B(ib,:)/(fs/2));
        sm = max(round(2*fs/fc(ib)), round(0.002*fs));
        for i = 1:ne
            X = filtfilt(b, a, double(ev.W(:,hch{g},k(i))));
            E = movmean(sum(abs(hilbert(X)).^2, 2), sm);
            nE = mean(E(t >= -0.095 & t < -0.015));
            iS = find(t >= tS(i)-0.01 & t < tS(i)+0.03); [pk, j] = max(E(iS)); ts = t(iS(j));
            ic = t >= ts+0.03 & t < ts+0.10; id = t >= ts+0.02 & t < ts+0.12;
            if mean(E(ic)) < 5*nE, continue; end
            OK(i,ib) = true;
            CL(i,ib) = log10(mean(E(ic))/pk);
            c1 = polyfit(t(id), log(E(id)), 1); BD(i,ib) = -c1(1);
            tl = t(id);                                   % lapse time from origin
            c2 = polyfit(tl, log(E(id).*tl.^2), 1); QI(i,ib) = -c2(1)/(2*pi*fc(ib));
            iw = t >= ts-0.03 & t < ts+0.03; w = E(iw)/sum(E(iw)); tw = t(iw);
            mu = sum(w.*tw); SW(i,ib) = sqrt(sum(w.*(tw-mu).^2));
        end
    end
    res(g).CL = CL; res(g).BD = BD; res(g).QI = QI; res(g).SW = SW; res(g).OK = OK; res(g).R = R; res(g).k = k;
end
% ---------------- summary + break-point test ----------------
rows = {}; brk = {};
for g = 1:2
    for ib = 1:nb
        o = res(g).OK(:,ib); q = @(v) prctile(v(o), [25 50 75]);
        a = q(res(g).CL(:,ib)); b = q(res(g).BD(:,ib)); c = q(res(g).QI(:,ib)); d = q(1e3*res(g).SW(:,ib));
        rows(end+1,:) = {char('A'+g-1), fc(ib), nnz(o), a(2), a(1), a(3), b(2), c(2), c(1), c(3), d(2)}; %#ok<SAGROW>
    end
    for v = {'CL','QI'}
        Y = res(g).(v{1}); X = repmat(log10(fc), size(Y,1), 1); m = res(g).OK & isfinite(Y);
        if strcmp(v{1},'QI'), m = m & Y > 0; Y = log10(max(Y, eps)); end
        x = X(m); y = Y(m); n = numel(y);
        p1 = polyfit(x, y, 1); r1 = y - polyval(p1, x); bic1 = n*log(mean(r1.^2)) + 2*log(n);
        best = inf; xb = NaN;
        for xc = log10(fc(2:end-1))                      % continuous two-segment (hinge) fit
            A = [ones(n,1) x max(x - xc, 0)]; cc = A\y; r2 = y - A*cc;
            bb = n*log(mean(r2.^2)) + 4*log(n);
            if bb < best, best = bb; xb = xc; cb = cc; end
        end
        brk(end+1,:) = {char('A'+g-1), v{1}, p1(1), bic1, best, best - bic1, 10^xb, cb(2), cb(2)+cb(3), n}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'sensor','fc_Hz','n','coda_level_med','cl_q25','cl_q75', ...
    'decay_per_s_med','Qc_inv_med','qi_q25','qi_q75','S_width_ms'});
disp(T); writetable(T, '../results/table_coda_multioctave.csv');
Tb = cell2table(brk, 'VariableNames', {'sensor','measure','slope_single','BIC_single','BIC_break', ...
    'dBIC_break_minus_single','f_break_Hz','slope_below','slope_above','n'});
disp(Tb); writetable(Tb, '../results/table_coda_break_test.csv');
save ../results/coda_multioctave.mat res fc B

% ---------------- figure ----------------
f1 = figure('Position',[30 30 1500 460], 'Visible','off');
lab = {'coda level log_{10}(E_{coda}/E_{S})', 'Q_c^{-1} (single scattering)', 'rms S-pulse width (ms)'};
fld = {'CL','QI','SW'}; sc = [1 1 1e3];
for p = 1:3
    subplot(1,3,p); hold on
    for g = 1:2
        Y = res(g).(fld{p})*sc(p); Y(~res(g).OK) = NaN;
        Q = prctile(Y, [25 50 75]);
        errorbar(fc*(1 + 0.04*(g-1)), Q(2,:), Q(2,:)-Q(1,:), Q(3,:)-Q(2,:), '-o', 'LineWidth', 1.5, 'MarkerFaceColor','w');
    end
    set(gca,'XScale','log'); if p == 2, set(gca,'YScale','log'); end
    xline([36 11], ':', {'k a = 1 (a = 15 m)','k a = 1 (a = 50 m)'}, 'FontSize', 7);
    xline([330 480 920], 'Color', [0.8 0.8 0.8]);
    xlim([15 2000]); grid on; xlabel('octave-band centre frequency (Hz)'); ylabel(lab{p});
    legend('sensor A','sensor B','Location','best');
end
sgtitle('FORGE 2022 MEQ S coda over two decades of frequency (56-32 sensors; grey: suspected instrument resonances)');
exportgraphics(f1, '../figs/fig_coda_multioctave.png', 'Resolution', 140); close(f1);
