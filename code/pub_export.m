function pub_export(f, file, wcm, hcm, fs)
% PUB_EXPORT  Export a figure at its print size for publication: width wcm x height hcm (cm;
% 17.4 cm = double column), axis/tick font fs pt (default 8), legends fs-1, panel labels fs+1 bold,
% 300 dpi PNG plus a vector PDF next to it.
if nargin < 5, fs = 8; end
set(f, 'Units', 'centimeters', 'Position', [1 1 wcm hcm], 'PaperUnits', 'centimeters', 'PaperSize', [wcm hcm], 'Color', 'w');
ax = findall(f, 'Type', 'axes');
set(ax, 'FontSize', fs, 'LineWidth', 0.6, 'TickLength', [0.015 0.015], 'FontName', 'Helvetica');
for a = ax'
    set([a.XLabel a.YLabel], 'FontSize', fs);
end
set(findall(f, 'Type', 'legend'), 'FontSize', fs - 1);
set(findall(f, 'Type', 'colorbar'), 'FontSize', fs - 1);
tx = findall(f, 'Type', 'text');
for t = tx'
    if strcmp(t.FontWeight, 'bold') && numel(t.String) <= 3 && endsWith(char(t.String), ')')
        set(t, 'FontSize', fs + 1);                  % panel labels a), b) ...
    elseif t.FontSize > fs
        set(t, 'FontSize', fs);
    end
end
drawnow;
exportgraphics(f, file, 'Resolution', 300);
[p, n] = fileparts(file);
exportgraphics(f, fullfile(p, [n '.pdf']), 'ContentType', 'vector');
end
