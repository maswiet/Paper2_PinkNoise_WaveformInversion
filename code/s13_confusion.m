% S13  Shape-resolution confusion matrix at equal sigma (0.13).
% Truths: 3 families x 3 seeds (33, 44, 55). Candidates: 2 realisations per
% family (seed 11 and 33); a candidate sharing the truth's seed is excluded
% (then only seed 11 is used for every family, keeping the comparison fair).
% Family score = mean Phi over its candidate realisations; also the same for FORGE data.
maxNumCompThreads(6);
here = fileparts(mfilename('fullpath')); out = fullfile(here,'..','results');
fam  = {'pink','exponential','Gaussian'};
cand = {{'P6','R6'}, {'E15','E15b'}, {'G50','G50b'}};       % seed 11, seed 33
truths = {'R6','R6c','R6d'; 'E15b','E15c','E15d'; 'G50b','G50c','G50d'};
tseed  = [33 44 55];
for t = truths(:)'
    f = fullfile(out, sprintf('feat_%s_forge.mat', t{1}));
    if ~exist(f,'file'), s07_model_features(t{1}); end
end
allc = [cand{:}];
rows = {}; CM = zeros(3); margin = NaN(3,3);
for i = 1:3
    for j = 1:3
        tr = truths{i,j};
        R = s09_compare_models(tr, allc);
        phi = containers.Map({R.name}, num2cell([R.PHI]));
        sc = zeros(1,3);
        for k = 1:3
            use = cand{k}; if tseed(j) == 33, use = use(1); end   % drop same-seed candidate
            sc(k) = mean(cellfun(@(m) phi(m), use));
        end
        [~, win] = min(sc); CM(i,win) = CM(i,win) + 1;
        s = sort(sc); margin(i,j) = sc(i) - min(sc(setdiff(1:3,i)));   % <0: correct family ahead
        rows(end+1,:) = {fam{i}, tr, sc(1), sc(2), sc(3), fam{win}}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'truth_family','truth','score_pink','score_exponential','score_Gaussian','selected'});
disp(T); writetable(T, fullfile(out,'table_confusion_runs.csv'));
C = array2table(CM, 'VariableNames', strcat('sel_', fam), 'RowNames', strcat('true_', fam));
disp(C); writetable(C, fullfile(out,'table_confusion.csv'), 'WriteRowNames', true);
fprintf('overall correct: %d / 9\n', trace(CM));
fprintf('margin (true-family score minus best other; negative = correct):\n'); disp(margin);
% FORGE data, both candidate realisations per family
R = s09_compare_models('data', allc);
phi = containers.Map({R.name}, num2cell([R.PHI]));
D = table(fam', cellfun(@(c) phi(c{1}), cand)', cellfun(@(c) phi(c{2}), cand)', ...
    cellfun(@(c) mean(cellfun(@(m) phi(m), c)), cand)', 'VariableNames', {'family','PHI_seed11','PHI_seed33','mean'});
disp(D); writetable(D, fullfile(out,'table_data_families.csv'));
save(fullfile(out,'confusion.mat'), 'CM', 'margin', 'T', 'D');
