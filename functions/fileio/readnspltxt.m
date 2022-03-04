function [data, fs] = readnspltxt(filename)
%% READNSPLTXT Read nspl .txt data into matlab.
% Use as£º
%   [data, fs] = readnspltxt(filename);
% Input: 
%   - filename: data file name
% Output:
%   - data: data matrix
%   - fs: samplingrate

fid = fopen(filename);

fs = fscanf(fid,'%f', [1 1]);

nChn = fscanf(fid,'%d', [1 1]);

nPoint = fscanf(fid,'%d', [1 1]);

tmp = fscanf(fid,'%d', [1 1]);

%read data
format = '%';
for i = 1:nChn-1
    format = strcat(format,'f %');
end
format = strcat(format,'f\n');

data=fscanf(fid,format,[nChn inf]);
data= data';

fclose(fid);

end
