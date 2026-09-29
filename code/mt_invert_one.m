function [vr, res, pred, mt, shift] = mt_invert_one(obs, gf, win, gsel, ROT, fs)
% MT_INVERT_ONE  Linear moment-tensor waveform inversion of one event (see S08).
% linear MT inversion with per-sensor shift search (+-6 ms)
if gsel > 0, sens = gsel; Rg = {ROT, ROT}; else, sens = 1:2; Rg = ROT; end
if ~iscell(Rg), Rg = {Rg, Rg}; end
cols = []; for g = sens, cols = [cols 2*g-1 2*g]; end %#ok<AGROW>
% rotate Green's functions into tool frame
gr = gf;
for g = sens
    for j = 1:6, gr(:,2*g-1:2*g,j) = gf(:,2*g-1:2*g,j) * Rg{g}'; end
end
sh = round(-0.006*fs):round(0.001*fs):round(0.006*fs);
shift = zeros(1,2);
for it = 1:2
    A = []; d = [];
    for g = sens
        for c = 2*g-1:2*g
            w = win(:,c);
            A = [A; squeeze(circshift(gr(w,c,:),0))]; %#ok<AGROW>
            ww = circshift(w, -shift(g));
            d = [d; obs(ww,c)]; %#ok<AGROW>
        end
    end
    mt = A\d;
    % update shifts by cross-correlation of prediction with observation
    for g = sens
        best = -inf;
        for s = sh
            p = 0;
            for c = 2*g-1:2*g
                w = win(:,c); pr = squeeze(gr(w,c,:))*mt;
                p = p + obs(circshift(w,-s),c)'*pr;
            end
            if p > best, best = p; shift(g) = s; end
        end
    end
end
A = []; d = [];
for g = sens, for c = 2*g-1:2*g
    w = win(:,c); A = [A; squeeze(gr(w,c,:))]; d = [d; obs(circshift(w,-shift(g)),c)]; %#ok<AGROW>
end, end
mt = A\d; vr = 1 - sum((d - A*mt).^2)/sum(d.^2);
% full-trace residuals (for coda analysis)
pred = zeros(size(obs)); res = pred;
for g = sens, for c = 2*g-1:2*g
    pred(:,c) = squeeze(gr(:,c,:))*mt;
    res(:,c) = circshift(obs(:,c), -shift(g)) - pred(:,c);
end, end
end
