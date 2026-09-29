function RES = s09_compare_models(obsname, models, mech)
% S09  Statistical waveform inversion / model-class comparison.
% obsname : 'data' (FORGE 2022) or a model name used as synthetic truth ('T1','T2')
% models  : candidate model names
% Misfit per feature = two-sample KS distance between observed and model
% feature distributions (sensor x band); envelope misfit = RMS log10 diff of
% median S-envelopes; TT misfit = |log(std ratio)| + |corr diff| (<20 m).
if nargin < 3, mech = 'forge'; end
G = forge_setup();
ev = load(fullfile(fileparts(mfilename('fullpath')),'..','data','forge2022_events.mat')); ev = ev.ev;
nm = {'clevel2','cdecay2','swidth'};           % mechanism-robust, S-normalised
use = [1 1; 1 0];                              % sensor x band used (B 20-40 Hz excluded: borehole resonance)
O = get_obs(obsname, mech, G, ev);
RES = struct();
for m = 1:numel(models)
    L = load(fullfile(fileparts(mfilename('fullpath')),'..','results', sprintf('feat_%s_%s.mat', models{m}, mech)));
    ks = []; envm = []; ttm = [];
    for g = 1:2
        for ib = 1:2
            if ~use(g,ib), continue; end
            Tm = feature_table(L.M(g,ib).F, L.M(g,ib).R, L.M(g,ib).tp, L.M(g,ib).k, G.evxyz);
            To = O(g,ib);
            % restrict to common events
            [~, io, im] = intersect(To.k, Tm.k);
            for j = 1:numel(nm)
                a = To.(nm{j})(io); b = Tm.(nm{j})(im);
                a = a(isfinite(a)); b = b(isfinite(b));
                [~,~,d] = kstest2(a, b); ks(end+1) = d; %#ok<AGROW>
            end
            eo = median(To.envS(:,io),2,'omitnan'); em = median(Tm.envS(:,im),2,'omitnan');
            lag = (-0.08:1e-3:0.12)'; wl = lag >= 0 & lag <= 0.10;   % post-peak coda envelope
            envm(end+1) = sqrt(mean((log10(eo(wl)) - log10(em(wl))).^2,'omitnan')); %#ok<AGROW>
            if ib == 1
                ttm(end+1) = abs(log(Tm.ttres_std/To.ttres_std)) + ...
                    mean(abs(Tm.ttcorr(1:2) - To.ttcorr(1:2)),'omitnan'); %#ok<AGROW>
                RES(m).ttstd(g) = Tm.ttres_std; RES(m).ttcorr(g,:) = Tm.ttcorr;
            end
        end
    end
    RES(m).name = models{m}; RES(m).ks = ks; RES(m).KS = mean(ks);
    RES(m).ENV = mean(envm); RES(m).TT = mean(ttm);
    RES(m).PHI = RES(m).KS + RES(m).ENV;   % TT reported separately (catalogue location error dominates data)
end
fprintf('\nObserved: %s  (mechanisms: %s)   obs TT std A/B = %.2f/%.2f ms, corr<10m A/B = %.2f/%.2f\n', ...
    obsname, mech, 1e3*O(1,1).ttres_std, 1e3*O(2,1).ttres_std, O(1,1).ttcorr(1), O(2,1).ttcorr(1));
fprintf('%-4s %7s %7s %7s %7s   %s\n','mod','KS','ENV','TT','PHI','TT std A/B (ms) | corr<10m A/B');
for m = 1:numel(RES)
    fprintf('%-4s %7.3f %7.3f %7.3f %7.3f   %.2f/%.2f | %.2f/%.2f\n', RES(m).name, RES(m).KS, RES(m).ENV, ...
        RES(m).TT, RES(m).PHI, 1e3*RES(m).ttstd, RES(m).ttcorr(1,1), RES(m).ttcorr(2,1));
end
out = fullfile(fileparts(mfilename('fullpath')),'..','results');
save(fullfile(out, sprintf('compare_%s_%s.mat', obsname, mech)), 'RES', 'O');
end

function O = get_obs(obsname, mech, G, ev)
if strcmp(obsname, 'data')
    Df = load(fullfile(fileparts(mfilename('fullpath')),'..','data','forge_features.mat'));
    for g = 1:2, for ib = 1:2
        d = Df.D(g,ib);
        [~, loc] = ismember(d.k, G.ev);
        O(g,ib) = feature_table(d.F, d.R, d.tpl, d.k, G.evxyz(loc,:)); %#ok<AGROW>
    end, end
else
    L = load(fullfile(fileparts(mfilename('fullpath')),'..','results', sprintf('feat_%s_%s.mat', obsname, mech)));
    for g = 1:2, for ib = 1:2
        O(g,ib) = feature_table(L.M(g,ib).F, L.M(g,ib).R, L.M(g,ib).tp, L.M(g,ib).k, G.evxyz); %#ok<AGROW>
    end, end
end
end
