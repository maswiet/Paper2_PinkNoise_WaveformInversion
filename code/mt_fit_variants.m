function R = mt_fit_variants(obs, gf, win, ROT, fs, variants, DCgrid)
% MT_FIT_VARIANTS  Moment-tensor waveform fit of one event with explicit time shifts.
%   obs  [nt x 4]  observed horizontals (sensor A: cols 1-2, B: cols 3-4), model time axis
%   gf   [nt x 4 x 6] Green's functions for the six MT basis tensors (FD frame E,N)
%   win  [nt x 4]  logical fit window per trace
%   ROT  {2}       2x2 tool rotation per sensor
%   variants: cellstr of 'full' | 'dev' | 'dc'
%   DCgrid [6 x K] basis coefficients of unit double couples (for 'dc')
% Per-sensor shifts s (+-6 ms) are applied to the OBSERVED traces: obsS = circshift(obs, s),
% i.e. obsS(i) = obs(i - s); the prediction and obsS share the model time axis.
% Returns R.(variant) with fields pred [nt x 4], obsS [nt x 4], vr, m (6x1), and R.cond,
% R.shift (samples, from the full-MT fit).
gr = gf;
for g = 1:2, for j = 1:6, gr(:,2*g-1:2*g,j) = gf(:,2*g-1:2*g,j) * ROT{g}'; end, end
sh = round(-0.006*fs):round(0.001*fs):round(0.006*fs);
shift = [0 0];
% ---- full MT with iterative shift search ----
for it = 1:3
    [A, d] = stack(gr, obs, win, shift);
    m = A \ d;
    for g = 1:2
        best = -inf;
        for s = sh
            p = 0;
            for c = 2*g-1:2*g
                w = win(:,c); o = circshift(obs(:,c), s);
                p = p + o(w)' * (squeeze(gr(w,c,:)) * m);
            end
            if p > best, best = p; shift(g) = s; end
        end
    end
end
[A, d] = stack(gr, obs, win, shift);
An = A ./ vecnorm(A); R.cond = cond(An); R.shift = shift;
obsS = obs; for g = 1:2, for c = 2*g-1:2*g, obsS(:,c) = circshift(obs(:,c), shift(g)); end, end
G2 = @(mm) full_pred(gr, mm);
for v = variants
    switch v{1}
        case 'full'
            mm = A \ d;
        case 'dev'                                   % trace-free: minimise |Am-d| s.t. m1+m2+m3 = 0
            ml = A \ d; N = (A'*A) \ eye(6); c = [1 1 1 0 0 0]';
            mm = ml - N*c * ((c'*N*c) \ (c'*ml));
        case 'dc'                                    % grid search over unit double couples, LS scalar moment
            P = A * DCgrid; num = (P' * d); den = sum(P.^2, 1)';
            [~, k] = max(num.^2 ./ den); mm = DCgrid(:,k) * (num(k)/den(k));
    end
    r.m = mm; r.vr = 1 - sum((d - A*mm).^2)/sum(d.^2);
    r.pred = G2(mm); r.obsS = obsS;
    R.(v{1}) = r;
end
end

function [A, d] = stack(gr, obs, win, shift)
A = []; d = [];
for g = 1:2
    for c = 2*g-1:2*g
        w = win(:,c); o = circshift(obs(:,c), shift(g));
        A = [A; squeeze(gr(w,c,:))]; d = [d; o(w)]; %#ok<AGROW>
    end
end
end

function p = full_pred(gr, mm)
p = zeros(size(gr,1), 4);
for c = 1:4, p(:,c) = squeeze(gr(:,c,:)) * mm; end
end
