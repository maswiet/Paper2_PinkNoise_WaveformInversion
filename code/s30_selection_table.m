% S30  Event-selection flow and the distributions at each step (reviewer point 8).
% Catalogue (36 641 events) -> 3556 SEG-2 files (18-24 Apr 2022) -> the 800 largest events with
% M >= -1 whose 1-s window lies inside a file (extraction cap, ordered by magnitude) -> clean P
% pick on at least one sensor -> inside the FD interior box (555) -> clean P picks on both sensors
% (416, per-event coda test) -> 150 highest P-SNR (original subset).
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
ev = load(fullfile(here,'..','data','forge2022_events.mat')); ev = ev.ev;
L = load(fullfile(here,'..','data','tt_inversion.mat'), 'obs'); obs = L.obs;
Df = load(fullfile(here,'..','data','forge_features.mat'));
G = forge_setup(); S = load(fullfile(out, 's19_H.mat'));
n = numel(ev.M); all800 = (1:n)';
picked = unique(obs(:,1));
both = intersect(Df.D(1,1).k, Df.D(2,1).k); both = both(ismember(both, G.ev));
top = S.ev(S.top150);
sets = {'extracted (M >= -1, largest 800)', all800; 'clean P pick, any sensor', picked; ...
        'inside FD box (either sensor)', G.ev(:); 'both sensors (coda test)', both(:); '150 highest P-SNR', top(:)};
Rm = sqrt(sum((ev.xyz - mean(G.sens,1)).^2, 2));
snrS = ev.snrS;
rows = {};
for s = 1:size(sets,1)
    k = sets{s,2};
    rows(end+1,:) = {sets{s,1}, numel(k), median(ev.M(k)), prctile(ev.M(k),10), prctile(ev.M(k),90), ...
        median(Rm(k)), prctile(Rm(k),10), prctile(Rm(k),90), 100*mean(ev.stage(k) == 3), median(snrS(k),'omitnan')}; %#ok<SAGROW>
end
T = cell2table(rows, 'VariableNames', {'step','n','M_med','M_p10','M_p90','R_med_m','R_p10_m','R_p90_m','pct_stage3','cat_snrS_med'});
disp(T); writetable(T, fullfile(out, 'table_event_selection.csv'));
% observed coda ratio: all 416 versus the 150 subset
oc = log10(S.obsratio(:,:,1));
C = table({'both sensors (416)'; '150 highest P-SNR'; 'remaining 266'}, ...
    [median(oc(:,1),'omitnan'); median(oc(S.top150,1),'omitnan'); median(oc(~S.top150,1),'omitnan')], ...
    [median(oc(:,2),'omitnan'); median(oc(S.top150,2),'omitnan'); median(oc(~S.top150,2),'omitnan')], ...
    'VariableNames', {'subset','obs_log10_coda_S_A','obs_log10_coda_S_B'});
disp(C); writetable(C, fullfile(out, 'table_selection_coda.csv'));
