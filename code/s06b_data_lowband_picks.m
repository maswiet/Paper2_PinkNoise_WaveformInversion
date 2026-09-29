% S06b  Re-pick FORGE P onsets in the FD band (20-100 Hz), horizontals only,
% so travel-time residual statistics are comparable with the synthetics.
load ../data/forge2022_events.mat; load ../data/forge_features.mat
G = forge_setup(); fs0 = ev.fs; tr = (-0.05:1/fs0:0.05)'; ar = (pi*G.f0*tr).^2; ric = (1-2*ar).*exp(-ar);
for i = 1:size(ev.W,3), ev.W(:,:,i) = single(conv2(double(ev.W(:,:,i)), ric, 'same')); end
fs = ev.fs; t = (0:size(ev.W,1)-1)'/fs - ev.pre;
hch = {[2 3],[5 6]};
for g = 1:2
    k = D(g,1).k; tph = D(g,1).tp; tpl = zeros(size(k));
    for i = 1:numel(k)
        tpl(i) = aic_pick(ev.W(:,hch{g},k(i)), t, [20 100], fs, tph(i)-0.02, tph(i)+0.02);
    end
    D(g,1).tpl = tpl; D(g,2).tpl = tpl;
    fprintf('sensor %c: low-band minus HF pick: median %.2f ms, MAD %.2f ms\n', char('A'+g-1), ...
        1e3*median(tpl-tph), 1e3*median(abs(tpl-tph-median(tpl-tph))));
end
save ../data/forge_features.mat D NZ bands
