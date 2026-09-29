function run_sgt_hf(name, nthreads)
% RUN_SGT_HF  High-frequency strain-Green's-tensor runs (dx 2.5 m, Ricker 90 Hz),
% sensor B only, 3 force directions. Output: ../sgt_hf/<name>.mat, E{2,n} at 0.25 ms.
if nargin > 1, maxNumCompThreads(nthreads); end
G = forge_setup_hf();
[M, info] = model_catalog(name, G); disp(info);
out = fullfile(fileparts(mfilename('fullpath')),'..','sgt_hf');
if ~exist(out,'dir'), mkdir(out); end
vmax = max(M.vp(:))*sqrt(1 + 2*max(M.eps(:),0));
cfl = vmax*G.dt/G.dx*sqrt(3)*7/6;
if cfl > 0.97, error('CFL %.2f too high for %s', cfl, name); end
opt.dt = G.dt; opt.nt = G.nt; opt.nabs = G.nab; opt.erec = G.ievg;
fprintf('%s: grid %s, %d events, CFL %.2f, min vs %.0f m/s (%.1f pts/lambda at 160 Hz)\n', name, ...
    mat2str(G.n), numel(G.ev), cfl, min(M.vs(:)), min(M.vs(:))/160/G.dx);
E = cell(2,3);
for n = 1:3
    g = 2;
    src = struct('ix',G.isg(g,1),'iy',G.isg(g,2),'iz',G.isg(g,3),'f0',G.f0,'t0',G.t0,'type','force','dir',n);
    tic; s = fd3d_elastic_vti(M, src, G.isg(g,:), opt);
    E{g,n} = s.e(1:G.decim:end,:,:);
    fprintf('%s sensor B force %d: %.0f s\n', name, n, toc);
end
save(fullfile(out,[name '.mat']), 'E', 'info', '-v7.3');
end
