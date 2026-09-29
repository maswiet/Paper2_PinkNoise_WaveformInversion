function Z = s15_model_coda_multioctave(name, grid)
% S15  Multi-octave coda measures (as S14) on FD synthetics.
%   grid '5m'  : sgt/<name>.mat, sensors A and B, octave bands 25 and 50 Hz (real noise added)
%   grid 'hf'  : sgt_hf/<name>.mat, sensor B, octave bands 50 and 100 Hz (noise-free)
here = fileparts(mfilename('fullpath'));
if strcmp(grid, '5m')
    G = forge_setup(); S = synth_records(name, 'forge', true); sens = 1:2; fcs = [25 50];
else
    G = forge_setup_hf(); S = synth_records_hf(name); sens = 2; fcs = [50 100];
end
B = [fcs'/sqrt(2) fcs'*sqrt(2)];
if strcmp(grid, 'hf'), fs = 1/(G.dt*G.decim); else, fs = 1/G.dt; end   % both 4 kHz
Z = struct('name', name, 'grid', grid, 'fc', fcs);
for g = sens
    t = S(g).t; ne = size(S(g).H, 3);
    tS = S(g).R/G.vs;
    CL = NaN(ne, numel(fcs)); QI = CL; SW = CL;
    for i = 1:ne
        [CL(i,:), QI(i,:), SW(i,:)] = coda_band_measures(double(S(g).H(:,:,i)), t, tS(i), fs, B, fcs);
    end
    Z.CL{g} = CL; Z.QI{g} = QI; Z.SW{g} = SW; Z.k{g} = S(g).k;
end
save(fullfile(here,'..','results',sprintf('codamo_%s_%s.mat', name, grid)), '-struct', 'Z');
for g = sens
    fprintf('%s %s sensor %c: CL %s | Qc^-1 %s | SW(ms) %s\n', name, grid, char('A'+g-1), ...
        mat2str(median(Z.CL{g},'omitnan'),3), mat2str(median(Z.QI{g},'omitnan'),3), mat2str(1e3*median(Z.SW{g},'omitnan'),3));
end
end
