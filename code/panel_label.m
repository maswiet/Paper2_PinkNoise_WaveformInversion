function panel_label(ax, k, pos)
% PANEL_LABEL  Put 'a)', 'b)', ... inside the axes frame (top-left by default).
%   ax : axes handle; k : panel index (1 -> 'a)') or a char label; pos : 'tl' | 'tr' | 'bl' | 'br'
if nargin < 3, pos = 'tl'; end
if isnumeric(k), s = [char('a' + k - 1) ')']; else, s = k; end
switch pos
    case 'tl', xy = [0.02 0.97]; ha = 'left';  va = 'top';
    case 'tr', xy = [0.98 0.97]; ha = 'right'; va = 'top';
    case 'bl', xy = [0.02 0.03]; ha = 'left';  va = 'bottom';
    case 'br', xy = [0.98 0.03]; ha = 'right'; va = 'bottom';
end
text(ax, xy(1), xy(2), s, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 11, ...
    'HorizontalAlignment', ha, 'VerticalAlignment', va, 'BackgroundColor', 'w', 'Margin', 1, 'Clipping', 'off');
end
