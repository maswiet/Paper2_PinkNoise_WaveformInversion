function run_d2_windows(names)
% RUN_D2_WINDOWS  Inter-event coda coherence (S11/D2) in the late-coda windows S+40..100 ms and
% S+50..110 ms (second review), for the listed media ('data' = FORGE records).
maxNumCompThreads(3);
for i = 1:numel(names)
    for w = {[0.040 0.100], [0.050 0.110]}
        try, s11_shape_diagnostics(names{i}, 'forge', w{1}); catch e, fprintf(2, '%s failed: %s\n', names{i}, e.message); end
    end
end
disp('run_d2_windows done');
end
