% S22  Uncertainty of the sensor self-calibration (reviewer point 11):
% event bootstrap (B = 300) and leave-one-stage-out of the homogeneous travel-time fit
% (common well E,N; depths zA, zB; clock offsets; Vp), same screened picks as S05.
maxNumCompThreads(1);
load ../data/forge2022_events.mat; load ../data/tt_inversion.mat
xe = ev.xyz(obs(:,1),:); sg = obs(:,2); to = obs(:,3); st = ev.stage(obs(:,1));
q0 = TT(1).q; opts = optimset('MaxFunEvals',2e4,'MaxIter',2e4,'TolX',1e-4,'TolFun',1e-7,'Display','off');
fitq = @(ii, q) fminsearch(@(q) sum(abs(to(ii) - tth(q, xe(ii,:), sg(ii)))), q, opts);
rng(3); B = 300; Q = zeros(B, numel(q0));
for b = 1:B
    ii = randi(numel(to), numel(to), 1);
    Q(b,:) = fitq(ii, q0);
end
nm = {'E_m','N_m','zA_m','zB_m','offA_ms','offB_ms','Vp_m_s'}; sc = [1 1 1 1 1e3 1e3 1];
T = table(nm', (q0.*sc)', (std(Q).*sc)', (prctile(Q,2.5).*sc)', (prctile(Q,97.5).*sc)', ...
    'VariableNames', {'param','estimate','bootstrap_std','lo95','hi95'});
disp(T); writetable(T, '../results/table_calib_bootstrap.csv');
C = corrcoef(Q); disp('parameter correlation matrix (bootstrap):'); disp(round(C,2));
writematrix(round(C,3), '../results/table_calib_corr.csv');
rows = {};
for s = unique(st)'
    ii = find(st ~= s); q = fitq(ii, q0);
    rows(end+1,:) = [{sprintf('without stage %d', s)}, num2cell(q.*sc)]; %#ok<SAGROW>
end
T2 = cell2table(rows, 'VariableNames', [{'case'}, nm]); disp(T2); writetable(T2, '../results/table_calib_loso.csv');

function tt = tth(q, xe, sg)
xs = [q(1:3); q(1:2) q(4)]; off = q(5:6);
tt = off(sg)' + sqrt(sum((xe - xs(sg,:)).^2, 2))/q(7); tt = tt(:);
end
