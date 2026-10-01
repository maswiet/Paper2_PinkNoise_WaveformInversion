function s20_family_inference(nthreads)
% S20  Fair model-family comparison with uncertainty (reviewer points 5, 6, 10).
% Phi (explicit definition): for the three sensor-band combinations used (A 20-40, A 40-80,
% B 40-80 Hz; B 20-40 excluded, borehole plateau) and the three S-normalised features
% (coda level, coda decay, S-pulse width):
%   Phi = mean_{3 combos x 3 features} D_KS  +  mean_{3 combos} rms_{lag 0-100 ms}(log10 E_obs - log10 E_mod)
% with E the median S envelope; all terms equally weighted; events equally weighted.
% (a) Event bootstrap (B = 200, same resample for every model): Phi distribution per model.
% (b) Holdout: each family's best (sigma, shape parameter, seed) is chosen on a training half of
%     the events and scored on the other half (20 random splits + near/far split).
if nargin > 0, maxNumCompThreads(nthreads); end
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
G = forge_setup();
% ---------------- model table ----------------
M = {'H','smooth',NaN,NaN; 'L','smooth',NaN,NaN; 'A','smooth',NaN,NaN; 'L6','deterministic',NaN,NaN; ...
     'GRAD3','smooth',NaN,NaN; 'FZ','deterministic',NaN,NaN; 'VSD','pink-VpVs',0.13,1.5; 'VSO','pink-VpVs',0.13,1.5};
pk = {'P3',0.09,1.5; 'R3',0.09,1.5; 'P6',0.13,1.5; 'R6',0.13,1.5; 'P9',0.18,1.5; 'R5',0.18,1.5; ...
      'P8',0.09,1.8; 'X_pink_0.09_1.8_33',0.09,1.8; 'X_pink_0.13_1.8_11',0.13,1.8; 'X_pink_0.13_1.8_33',0.13,1.8; ...
      'P11',0.18,1.8; 'X_pink_0.18_1.8_33',0.18,1.8};
for i = 1:size(pk,1), M(end+1,:) = {pk{i,1},'pink',pk{i,2},pk{i,3}}; end %#ok<AGROW>
for fam = {'exp','gau'}
    for a = [15 50]
        for sg = [0.09 0.13 0.18]
            for sd = [11 33]
                nm = sprintf('X_%s_%.2f_%d_%d', fam{1}, sg, a, sd);
                if sg == 0.13 && sd == 11, nm = sprintf('%s%d', upper(fam{1}(1)), a); end        % E15 E50 G15 G50
                if sg == 0.13 && sd == 33 && ((strcmp(fam{1},'exp') && a == 15) || (strcmp(fam{1},'gau') && a == 50))
                    nm = sprintf('%s%db', upper(fam{1}(1)), a);                                     % E15b G50b
                end
                M(end+1,:) = {nm, [fam{1} num2str(a)], sg, a}; %#ok<AGROW>
            end
        end
    end
end
have = cellfun(@(n) exist(fullfile(out, sprintf('feat_%s_forge.mat', n)), 'file') == 2, M(:,1));
if any(~have), fprintf('missing features (skipped): %s\n', strjoin(M(~have,1)', ' ')); end
M = M(have,:); nmod = size(M,1);
% ---------------- observed and model feature tables ----------------
Df = load(fullfile(here,'..','data','forge_features.mat'));
use = {[1 1],[1 2],[2 2]}; fn = {'clevel2','cdecay2','swidth'};
O = cell(1,3); ev0 = [];
for u = 1:3
    g = use{u}(1); ib = use{u}(2); d = Df.D(g,ib);
    O{u} = struct('k', d.k, 'F', feats(d.F, fn), 'E', [d.F.envS]);
    ev0 = union(ev0, d.k);
end
MOD = cell(nmod, 3);
for m = 1:nmod
    L = load(fullfile(out, sprintf('feat_%s_forge.mat', M{m,1})));
    for u = 1:3
        g = use{u}(1); ib = use{u}(2); d = L.M(g,ib);
        MOD{m,u} = struct('k', d.k, 'F', feats(d.F, fn), 'E', [d.F.envS]);
    end
end
lag = (-0.08:1e-3:0.12)'; wl = lag >= 0 & lag <= 0.10;
phi = @(m, evs) phi_subset(O, MOD(m,:), evs, wl);
% ---------------- (a) full-sample Phi and event bootstrap ----------------
rng(5); B = 200; P0 = zeros(nmod,1); PB = zeros(B, nmod);
for m = 1:nmod, P0(m) = phi(m, ev0); end
for b = 1:B
    evs = ev0(randi(numel(ev0), numel(ev0), 1));
    for m = 1:nmod, PB(b,m) = phi(m, evs); end
end
T = table(M(:,1), M(:,2), cell2mat(M(:,3)), cell2mat(M(:,4)), P0, prctile(PB,2.5)', prctile(PB,97.5)', ...
    'VariableNames', {'model','family','sigma','shape','PHI','PHI_lo95','PHI_hi95'});
T = sortrows(T, 'PHI'); disp(T); writetable(T, fullfile(out, 'table_phi_bootstrap.csv'));
% family-level: best member per bootstrap sample
fams = unique(M(:,2), 'stable'); fb = zeros(B, numel(fams));
for f = 1:numel(fams), fb(:,f) = min(PB(:, strcmp(M(:,2), fams{f})), [], 2); end
[~, win] = min(fb, [], 2);
Tf = table(fams, arrayfun(@(f) min(P0(strcmp(M(:,2),fams{f}))), (1:numel(fams))'), median(fb)', prctile(fb,2.5)', prctile(fb,97.5)', ...
    accumarray(win, 1, [numel(fams) 1])/B*100, 'VariableNames', {'family','best_PHI','boot_median','lo95','hi95','pct_boot_wins'});
disp(Tf); writetable(Tf, fullfile(out, 'table_family_bootstrap.csv'));
% ---------------- (b) holdout ----------------
splits = cell(0,2); rng(8);
for r = 1:20, p = ev0(randperm(numel(ev0))); h = floor(numel(p)/2); splits(end+1,:) = {p(1:h), p(h+1:end)}; end %#ok<AGROW>
[~, loc] = ismember(ev0, G.ev); Rm = sqrt(sum((G.evxyz(loc,:) - mean(G.sens,1)).^2, 2));
splits(end+1,:) = {ev0(Rm <= median(Rm)), ev0(Rm > median(Rm))};
rows = {};
for s = 1:size(splits,1)
    for dirn = 1:2
        tr = splits{s, dirn}; te = splits{s, 3-dirn};
        for f = 1:numel(fams)
            ii = find(strcmp(M(:,2), fams{f})); ptr = arrayfun(@(m) phi(m, tr), ii);
            [~, j] = min(ptr); rows(end+1,:) = {s, dirn, fams{f}, M{ii(j),1}, ptr(j), phi(ii(j), te)}; %#ok<AGROW>
        end
    end
end
H = cell2table(rows, 'VariableNames', {'split','direction','family','chosen','PHI_train','PHI_test'});
writetable(H, fullfile(out, 'table_holdout.csv'));
S = groupsummary(H, 'family', {'median','min','max'}, 'PHI_test'); disp(S);
writetable(S, fullfile(out, 'table_holdout_summary.csv'));
end

function F = feats(Fs, fn)
F = zeros(numel(Fs), numel(fn));
for j = 1:numel(fn), F(:,j) = [Fs.(fn{j})]'; end
end

function p = phi_subset(O, MD, evs, wl)
ks = []; en = [];
for u = 1:3
    [tf, io] = ismember(evs, O{u}.k); [tf2, im] = ismember(evs, MD{u}.k); ok = tf & tf2;
    io = io(ok); im = im(ok);
    for j = 1:3
        a = O{u}.F(io, j); b = MD{u}.F(im, j); a = a(isfinite(a)); b = b(isfinite(b));
        [~,~,d] = kstest2(a, b); ks(end+1) = d; %#ok<AGROW>
    end
    eo = median(O{u}.E(:, io), 2, 'omitnan'); em = median(MD{u}.E(:, im), 2, 'omitnan');
    en(end+1) = sqrt(mean((log10(eo(wl)) - log10(em(wl))).^2, 'omitnan')); %#ok<AGROW>
end
p = mean(ks) + mean(en);
end
