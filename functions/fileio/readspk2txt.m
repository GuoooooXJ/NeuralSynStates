function [data, hdr] = readspk2txt(filename)
%% READSPK2TXT Read spike2 .txt data into matlab.
% Use as:
%   [data, hdr] = readspk2txt(filename);
% Input:
%   - filename: data file name
% Output:
%   - data: data matrix
%   - hdr: header information

    temp = importdata(filename);
    data = temp.data;
    
    hdr.channels = temp.colheaders;
    
    % Use the first 11 time points to estimate simpling rate.
    hdr.samplingrate = 10/(data(11,1)-data(1,1));
    hdr.ndim = size(data);
end
