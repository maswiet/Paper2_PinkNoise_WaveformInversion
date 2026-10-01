function run_sgt(name, nthreads, forces)
% RUN_SGT  Strain-Green's-tensor runs (reciprocity) for one model:
% 2 sensors x force directions (default 1:3; 1:2 = horizontal components only,
% sufficient for all horizontal-component analyses); strain rates at all event nodes.
% Output: ../sgt/<name>.mat  E{g,n} [nt x 6 x ne] single
if nargin > 1 && ~isempty(nthreads), maxNumCompThreads(nthreads); end
if nargin < 3, forces = 1:3; end
G = forge_setup();
[M, info] = model_catalog(name, G); disp(info);
out = fullfile(fileparts(mfilename('fullpath')),'..','sgt');
if ~exist(out,'dir'), mkdir(out); end
% automatic time sub-stepping for strong heterogeneity (CFL, O(4) staggered)
vmax = max(M.vp(:))*sqrt(1 + 2*max(M.eps(:),0));
nsub = max(1, ceil(vmax*G.dt/G.dx*sqrt(3)*7/6 / 0.95));
opt.dt = G.dt/nsub; opt.nt = G.nt*nsub; opt.nabs = G.nab; opt.erec = G.ievg;
fprintf('%s: vmax %.0f m/s, nsub %d\n', name, vmax, nsub);
E = cell(2,3);
for g = 1:2
    for n = forces
        src = struct('ix',G.isg(g,1),'iy',G.isg(g,2),'iz',G.isg(g,3), ...
            'f0',G.f0,'t0',G.t0,'type','force','dir',n);
        tic; s = fd3d_elastic_vti(M, src, G.isg(g,:), opt);
        E{g,n} = s.e(1:nsub:end,:,:);          % back to data sampling (content < 120 Hz)
        fprintf('%s sensor %d force %d: %.0f s\n', name, g, n, toc);
    end
end
save(fullfile(out,[name '.mat']), 'E', 'info', '-v7.3');
end
