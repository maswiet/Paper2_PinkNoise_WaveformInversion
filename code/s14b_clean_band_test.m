% S14b  Break test restricted to the resonance-free band 25-200 Hz
% (sensor B 25-Hz band excluded: known borehole plateau).
load ../results/coda_multioctave.mat
rows = {};
for g = 1:2
    use = fc <= 200; if g == 2, use = use & fc > 25; end
    for v = {'CL','QI','SW'}
        Y = res(g).(v{1})(:,use); m = res(g).OK(:,use) & isfinite(Y);
        X = repmat(log10(fc(use)), size(Y,1), 1);
        if strcmp(v{1},'QI'), m = m & Y > 0; Y = log10(max(Y, eps)); end
        if strcmp(v{1},'SW'), Y = log10(Y); end
        x = X(m); y = Y(m); n = numel(y);
        p1 = polyfit(x, y, 1); r1 = y - polyval(p1, x); bic1 = n*log(mean(r1.^2)) + 2*log(n);
        best = inf; fb = NaN; cb = [NaN NaN NaN];
        xcs = log10(fc(use)); xcs = xcs(2:end-1);
        for xc = xcs
            A = [ones(n,1) x max(x - xc, 0)]; cc = A\y; r2 = y - A*cc;
            bb = n*log(mean(r2.^2)) + 4*log(n);
            if bb < best, best = bb; fb = 10^xc; cb = cc; end
        end
        % bootstrap slope of single power law (events resampled)
        ev_i = repmat((1:size(m,1))', 1, nnz(use)); ev_i = ev_i(m); ue = unique(ev_i); s = zeros(500,1);
        for bs = 1:500
            pick = ue(randi(numel(ue), numel(ue), 1)); idx = cell2mat(arrayfun(@(e) find(ev_i == e), pick, 'uni', 0));
            pp = polyfit(x(idx), y(idx), 1); s(bs) = pp(1);
        end
        rows(end+1,:) = {char('A'+g-1), v{1}, sprintf('%g-%g', min(fc(use)), max(fc(use))), n, p1(1), prctile(s,2.5), prctile(s,97.5), best - bic1, fb, cb(2), cb(2)+cb(3)}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'sensor','measure','band_Hz','n','slope','slope_lo95','slope_hi95','dBIC_break','f_break','slope_below','slope_above'});
disp(T); writetable(T, '../results/table_coda_cleanband.csv');
