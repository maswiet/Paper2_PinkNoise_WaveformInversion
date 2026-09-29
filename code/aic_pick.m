function tp = aic_pick(X, t, band, fs, t1, t2)
% AIC_PICK  Onset pick (Maeda AIC) on multi-component energy within [t1 t2].
[b,a] = butter(3, band/(fs/2));
Y = filtfilt(b, a, double(X));
E = sqrt(sum(Y.^2, 2));
i = find(t >= t1 & t <= t2);
[~, im] = max(E(i)); i = i(1:max(im,12));       % up to the first strong peak
x = E(i); n = numel(x); aic = inf(n,1);
for m = 4:n-4
    aic(m) = m*log(var(x(1:m))+eps) + (n-m)*log(var(x(m+1:end))+eps);
end
[~, j] = min(aic); tp = t(i(j));
end
