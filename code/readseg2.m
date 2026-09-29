function [data, dt, hdr, fhdr] = readseg2(filename)
% READSEG2  Read a SEG-2 file (Pullan 1990). Returns data [NS x N] (double),
% sample interval dt (s), per-trace header strings hdr{N}, file header fhdr.
fid = fopen(filename,'r','ieee-le');
c = onCleanup(@() fclose(fid));
fread(fid,1,'int16');                 % block id
fread(fid,1,'int16');                 % revision
M = fread(fid,1,'uint16');            % size of trace pointer subblock
N = fread(fid,1,'uint16');            % number of traces
fseek(fid,32,'bof');
trp = fread(fid,M/4,'uint32');
fh = fread(fid, trp(1)-32-M, 'uint8=>char')';
fhdr = regexprep(fh, char(0), ' ');
data = []; hdr = cell(N,1); dt = NaN;
for i = 1:N
    fseek(fid, trp(i), 'bof');
    fread(fid,1,'uint16');                % trace id
    X  = fread(fid,1,'uint16');           % descriptor size
    fread(fid,1,'uint32');                % data block size
    NS = fread(fid,1,'uint32');
    DFC = fread(fid,1,'uint8');
    fseek(fid, trp(i)+32, 'bof');
    h = fread(fid, X-32, 'uint8=>char')';
    hdr{i} = regexprep(h, char(0), ' ');
    k = strfind(hdr{i}, 'SAMPLE_INTERVAL');
    if ~isempty(k), dt = sscanf(hdr{i}(k+15:end), '%f', 1); end
    if i == 1, data = zeros(NS, N); end
    fseek(fid, trp(i)+X, 'bof');
    switch DFC
        case 1, data(:,i) = fread(fid, NS, 'int16');
        case 2, data(:,i) = fread(fid, NS, 'int32');
        case 4, data(:,i) = fread(fid, NS, 'float32');
        case 5, data(:,i) = fread(fid, NS, 'float64');
        otherwise, error('unsupported DFC %d', DFC);
    end
end
end
