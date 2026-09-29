% S10  Export all numbers quoted in the manuscript to CSV (results/).
maxNumCompThreads(2);
out = '../results/';
sp = struct('H',[0 NaN],'L',[0 NaN],'A',[0 NaN],'P1',[0.045 1.5],'P2',[0.045 1.8],'P3',[0.09 1.5], ...
    'P4',[0.02 1.5],'P5',[0.045 1.2],'P6',[0.13 1.5],'P7',[0.09 1.2],'P8',[0.09 1.8],'P9',[0.18 1.5], ...
    'T1',[0.045 1.5],'R1',[0.045 1.5],'R2',[0.045 1.5],'R3',[0.09 1.5], ...
    'P10',[0.25 1.5],'P11',[0.18 1.8],'P12',[0.22 1.5],'R5',[0.18 1.5], ...
    'E15',[0.13 NaN],'E50',[0.13 NaN],'G15',[0.13 NaN],'G50',[0.13 NaN]);
for o = {'T1_forge','T2_forge','data_forge','data_random'}
    f = [out 'compare_' o{1} '.mat']; if ~exist(f,'file'), continue; end
    C = load(f); R = C.RES; n = numel(R);
    T = table({R.name}', zeros(n,1), zeros(n,1), [R.KS]', [R.ENV]', [R.PHI]', ...
        'VariableNames', {'model','sigma','p','KS','ENV','PHI'});
    T.class = repmat({''}, n, 1);
    for m = 1:n, a = sp.(R(m).name); T.sigma(m) = a(1); T.p(m) = a(2); T.class{m} = cls(R(m).name); end
    writetable(sortrows(T,'PHI'), [out 'table_misfit_' o{1} '.csv']);
end
for o = {'T1','T2','data'}
    f = [out 'mtinv_' o{1} '.mat']; if ~exist(f,'file'), continue; end
    L = load(f); W = L.W; n = numel(W);
    T = table({W.name}', arrayfun(@(w) median(w.vr), W)', arrayfun(@(w) prctile(w.vr,25), W)', ...
        arrayfun(@(w) prctile(w.vr,75), W)', arrayfun(@(w) median(w.pcoda), W)', arrayfun(@(w) median(w.rcoda), W)', ...
        'VariableNames', {'model','VR_med','VR_q25','VR_q75','pred_coda_over_S','unexpl_coda_over_S'});
    writetable(T, [out 'table_mt_' o{1} '.csv']);
end
L = load('../data/tt_inversion.mat'); TT = L.TT;
T = table({TT.cls}', [TT.mad]'*1e3, arrayfun(@(t) mean(abs(t.res))*1e3, TT)', [TT.aic]', [TT.bic]', ...
    arrayfun(@(t) mat2str(t.par,4), TT, 'uni', 0)', 'VariableNames', {'model','MAD_ms','meanabs_ms','AIC','BIC','params'});
writetable(T, [out 'table_tt.csv']);
T = table(TT(1).edges(1:end-1)', TT(1).edges(2:end)', TT(1).rho{1}', TT(1).rho{2}', ...
    'VariableNames', {'d_lo_m','d_hi_m','rho_A','rho_B'});
writetable(T, [out 'table_ttcorr.csv']);
% data feature medians and model medians (Table of physical features)
Df = load('../data/forge_features.mat'); nm = {'clevel2','cdecay2','swidth'};
rows = {}; mods = {'data','H','L','A','P1','P3','P6','P9','E15','E50','G15','G50'};
for g = 1:2, for ib = 1:2
    for m = 1:numel(mods)
        if strcmp(mods{m},'data'), F = Df.D(g,ib).F;
        else
            f = sprintf('%sfeat_%s_forge.mat', out, mods{m}); if ~exist(f,'file'), continue; end
            Lm = load(f); [~,im] = intersect(Lm.M(g,ib).k, Df.D(g,ib).k); F = Lm.M(g,ib).F(im);
        end
        v = cellfun(@(x) median([F.(x)],'omitnan'), nm);
        rows(end+1,:) = {char('A'+g-1), sprintf('%d-%d', Df.bands(ib,:)), mods{m}, v(1), v(2), v(3)*1e3}; %#ok<SAGROW>
    end
end, end
T = cell2table(rows, 'VariableNames', {'sensor','band_Hz','model','coda_level_log10','coda_decay_per_s','S_width_ms'});
writetable(T, [out 'table_features.csv']);
G = forge_setup(); Lg = load('../data/log5632_stats.mat'); LOG = Lg.LOG;
K = table({'n_events_box';'sensorA_depth_m';'sensorB_depth_m';'sensor_E_m';'sensor_N_m';'vp_m_s';'log_sig6';'log_sig50';'log_beta'; ...
    'nx';'ny';'nz';'dx_m';'dt_s';'nt';'f0_Hz'}, ...
    [numel(G.ev); G.sens(1,3); G.sens(2,3); G.sens(1,1); G.sens(1,2); G.vp; LOG.sig6; LOG.sig50; LOG.beta; G.n(:); G.dx; G.dt; G.nt; G.f0], ...
    'VariableNames', {'key','value'});
writetable(K, [out 'key_numbers.csv']);
disp('tables exported');

function c = cls(n)
if any(strcmp(n, {'H','L','A'})), c = 'smooth';
elseif any(n(1) == 'EG'), c = 'competitor';
else, c = 'pink'; end
end
