% S02  QC: timing alignment of extracted FORGE windows vs catalogue origin times
load ../data/forge2022_events.mat
fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
[b,a] = butter(4, [5 200]/(fs/2));
[~,o] = sort(ev.R);
figure('Position',[50 50 1400 800]);
for c = 1:6
    subplot(2,3,c); hold on
    E = zeros(numel(t), numel(o));
    for k = 1:numel(o)
        x = filtfilt(b,a,double(ev.W(:,c,o(k))));
        E(:,k) = abs(hilbert(x)); E(:,k) = E(:,k)/max(E(:,k));
    end
    imagesc(t, 1:numel(o), E'); axis tight; set(gca,'YDir','normal');
    plot(ev.R(o)/5800, 1:numel(o), 'w.', 'MarkerSize',2);
    plot(ev.R(o)/3350, 1:numel(o), 'r.', 'MarkerSize',2);
    title(sprintf('ch %d  (sorted by R)',c)); xlabel('t - origin (s)'); xlim([-0.1 0.6]);
end
colormap(parula);
exportgraphics(gcf,'../figs/qc_forge_timing.png','Resolution',110);
% stacked envelope per channel
Es = zeros(numel(t),6);
for c = 1:6
    for k = 1:numel(o)
        x = filtfilt(b,a,double(ev.W(:,c,k))); e = abs(hilbert(x)); Es(:,c) = Es(:,c) + e/max(e);
    end
end
[~,ip] = max(Es(t<0.25,:)); disp('peak time of stacked envelope per channel (s):'); disp(t(ip)');
fprintf('median R = %.0f m -> P %.3f s, S %.3f s\n', median(ev.R), median(ev.R)/5800, median(ev.R)/3350);
