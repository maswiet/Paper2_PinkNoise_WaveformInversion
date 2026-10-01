function post_revision()
% POST_REVISION  After queue_revision.sh: features (S07) and per-event coda test (S19) for all
% new models, then the fair family inference (S20), the convergence comparison (S26), sigma_eff
% (S29), summaries and revised figures. Each step is independent (errors are reported, not fatal);
% script steps run in the base workspace so they cannot overwrite this function's loop variables.
maxNumCompThreads(12);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
d = dir(fullfile(here,'..','sgt','*.mat')); names = erase({d.name}, '.mat');
new = names(startsWith(names, 'X_') | ismember(names, {'L6','GRAD3','FZ','VSD','VSO'}));
for i = 1:numel(new)
    try
        if ~exist(fullfile(out, sprintf('feat_%s_forge.mat', new{i})), 'file'), s07_model_features(new{i}); end
    catch e, fprintf(2, 'S07 %s failed: %s\n', new{i}, e.message);
    end
end
todo = new(~cellfun(@(n) exist(fullfile(out, sprintf('s19_%s.mat', n)), 'file') == 2, new));
try, if ~isempty(todo), s19_coda_prediction(todo, 12); end, catch e, fprintf(2, 'S19 failed: %s\n', e.message); end
steps = {@() s20_family_inference(12), @() evalin('base','s26_convergence'), @() evalin('base','s29_sigma_eff'), ...
         @() evalin('base','s27_summarise_s19'), @() evalin('base','s28_s_timing_check'), ...
         @() fig_seis_rev({'H','FZ','P1','P3'}), @() fig_rev_panels};
for k = 1:numel(steps)
    try, steps{k}(); catch e, fprintf(2, 'post step %d (%s) failed: %s\n', k, func2str(steps{k}), e.message); end
end
disp('post_revision done');
end
