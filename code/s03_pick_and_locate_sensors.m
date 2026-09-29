% S03  Automatic onset picks on both 56-32 3C sensors (chs 1-3 = A, 4-6 = B)
% and joint inversion for sensor position, clock offset and velocity using
% catalogue event locations (sensor self-calibration).
load ../data/forge2022_events.mat
load ../data/qc_snr.mat
fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
[b,a] = butter(4,[30 1500]/(fs/2));
ne = size(ev.W,3); grp = {1:3, 4:6};
tp = NaN(ne,2); ts = NaN(ne,2); q = zeros(ne,2);
for k = 1:ne
    for g = 1:2
        X = filtfilt(b,a,double(ev.W(:,grp{g},k)));
        E = sum(X.^2,2);
        nz = mean(E(t<-0.02)); if nz==0, continue; end
        q(k,g) = max(E(t>0 & t<0.4))/nz;
        % AIC onset on energy trace, window before the first strong peak
        i1 = find(t>-0.02,1); i2 = find(E > 0.2*max(E(t>0&t<0.4)) & t>0, 1);
        if isempty(i2) || i2-i1 < 20, continue; end
        x = sqrt(E(i1:i2)); n = numel(x); aic = inf(n,1);
        for m = 5:n-5
            aic(m) = m*log(var(x(1:m))+eps) + (n-m)*log(var(x(m+1:end))+eps);
        end
        [~,im] = min(aic); tp(k,g) = t(i1+im-1);
        % second arrival: max of horizontal-ish energy after first onset + 10 ms
        j = find(t > tp(k,g)+0.010 & t < tp(k,g)+0.20);
        Es = movmean(E(j), 20); [~,jm] = max(Es);
        ts(k,g) = t(j(jm));
    end
end
good = q > 1e3 & isfinite(tp);
fprintf('good picks: A %d, B %d\n', sum(good));

% ---- joint fit: t = t0 + |x_ev - x_s|/v, per sensor; robust (L1) ----
opts = optimset('MaxFunEvals',2e4,'MaxIter',2e4,'TolX',1e-3);
S = struct();
for g = 1:2
    k = find(good(:,g));
    xe = ev.xyz(k,:); to = tp(k,g);
    f = @(p) sum(abs(to - p(4) - sqrt(sum((xe - p(1:3)).^2,2))/p(5)));
    best = inf;
    for d0 = [1500 1900 2200]            % multistart in sensor depth
        p0 = [ev.rec(1:2) d0 0 5800];
        [p,fv] = fminsearch(f, p0, opts);
        if fv < best, best = fv; pb = p; end
    end
    res = to - pb(4) - sqrt(sum((xe - pb(1:3)).^2,2))/pb(5);
    fprintf('Sensor %s: E=%.0f N=%.0f Z=%.0f m, t0=%.4f s, v=%.0f m/s, MAD res=%.2f ms (n=%d)\n', ...
        char('A'+g-1), pb(1:3), pb(4), pb(5), 1e3*median(abs(res)), numel(k));
    S(g).p = pb; S(g).res = res; S(g).k = k;
end
% also S-P (clock independent) check with fitted positions
for g = 1:2
    k = S(g).k; R = sqrt(sum((ev.xyz(k,:) - S(g).p(1:3)).^2,2));
    sp = ts(k,g) - tp(k,g);
    c = polyfit(R, sp, 1);
    fprintf('Sensor %s: second-arrival minus onset vs R slope %.2e s/m (S-P theory %.2e)\n', ...
        char('A'+g-1), c(1), 1/3350-1/5800);
end
save ../data/forge_picks_sensors.mat tp ts q good S
