% S01  Match GES April-2022 catalogue events (16A stages 1-3) to the 56-32 PSS
% continuous SEG-2 files and extract 3C event windows.
% Output: ../data/forge2022_events.mat
root = '../../56-32/';
f = dir([root 'Seg2byGES_56-32_PSS_*.sg2']);
nf = numel(f); t0 = NaN(nf,1);
for i = 1:nf                                   % file start times from file header
    fid = fopen([root f(i).name],'r'); b = fread(fid,600,'uint8=>char')'; fclose(fid);
    b(b==char(0)) = ' ';
    d = regexp(b,'ACQUISITION_DATE (\d+)/(\d+)/(\d+)','tokens','once');
    tm = regexp(b,'ACQUISITION_TIME (\d+):(\d+):(\d+)','tokens','once');
    fr = regexp(b,'ACQUISITION_SECOND_FRACTION ([\d\.]+)','tokens','once');
    t0(i) = datenum(str2double(d{3}),str2double(d{2}),str2double(d{1}), ...
        str2double(tm{1}),str2double(tm{2}),str2double(tm{3})+str2double(fr{1}));
end
fprintf('%d files, %s to %s\n', nf, datestr(min(t0)), datestr(max(t0)));

% ---- catalogue ----
C = [];
for s = 1:3
    T = readtable(sprintf('%sFORGE16AApril22AllStage %d.xlsx', root, s));
    od = datenum(T.OriginDate,'dd/mm/yyyy');
    ot = cellfun(@(x) sscanf(x,'%d:%d:%f')', T.OriginTime, 'uni', 0);
    ot = cell2mat(ot);
    to = od + (ot(:,1) + ot(:,2)/60 + ot(:,3)/3600)/24;
    C = [C; to T.X T.Y T.Depth T.MomMag T.Error T.PS_N T.SS_N s*ones(height(T),1)]; %#ok<AGROW>
end
ok = all(isfinite(C(:,1:5)),2); C = C(ok,:);
fprintf('%d catalogue events\n', size(C,1));

% ---- geometry (ft -> m), 56-32 PSS level from GES 2024 receiver file ----
ft = 0.3048;
rec = [2860.8 -10.3 5012.4]*ft;                % E N depth
xyz = [C(:,2) C(:,3) C(:,4)]*ft;               % X=E, Y=N, Depth
R = sqrt(sum((xyz - rec).^2,2));

% ---- match & extract ----
pre = 0.10; win = 1.0; fs = 4000;
sel = find(C(:,5) >= -1.0);                    % magnitude floor for SNR
[~,o] = sort(C(sel,5),'descend'); sel = sel(o);
W = zeros(win*fs, 6, 0, 'single'); keep = [];
for j = sel'
    tst = C(j,1) - pre/86400;
    i = find(t0 <= tst & t0 + 60/86400 >= tst + win/86400, 1);
    if isempty(i), continue; end
    [d, dt] = readseg2([root f(i).name]);
    i1 = round((tst - t0(i))*86400/dt) + 1;
    w = d(i1:i1+win*fs-1, :);
    W(:,:,end+1) = single(w - mean(w)); %#ok<SAGROW>
    keep(end+1) = j; %#ok<SAGROW>
    if numel(keep) >= 800, break; end
end
ev.t = C(keep,1); ev.xyz = xyz(keep,:); ev.M = C(keep,5); ev.err = C(keep,6)*ft;
ev.snrP = C(keep,7); ev.snrS = C(keep,8); ev.stage = C(keep,9); ev.R = R(keep);
ev.rec = rec; ev.fs = fs; ev.pre = pre; ev.W = W;
if ~exist('../data','dir'), mkdir('../data'); end
save('../data/forge2022_events.mat','ev','-v7.3');
fprintf('extracted %d events, R = %.0f-%.0f m, M = %.2f..%.2f\n', numel(keep), min(ev.R), max(ev.R), min(ev.M), max(ev.M));
