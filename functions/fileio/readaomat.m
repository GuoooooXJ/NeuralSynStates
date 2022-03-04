function [data,hdr] = readaomat(filename,channel,type,section)
%READAOMAT Read AlphaOmega .mat format.
%  Use as:
%    [data,hdr] = readaomat(filename,channel,type,section);
%  Input:
%    - filename, input file name.
%    - channel, select channels e.g. [1,3,7].
%    - type, 'LFP' or 'RAW'.
%    - section, optional, select a time range in e.g. [10,60] (s).
%  Output:
%    - data, data matrix, each column as a channel.
%    - hdr, structure, saving header information.
%
%  Author   :  NIE Yingnan
%  Created  : Sept 9, 2020
%  Modified : Sept 9, 2020

if nargin==3
    section = 'all';
end

% Load fs
type = strcat('C',type);

fsVar = strcat(type, '_', num2str(channel(1), '%03d'),'_KHz');
load(filename, fsVar);
fs = eval(fsVar)*1000;

nChan = length(channel);

chanName = [];

if isequal(section,'all')
    % Load channel lengths
    chanLength = zeros(nChan,1); % The length of channels maybe different.
    for i=1:nChan
        chanVar = strcat(type, '_', num2str(channel(i), '%03d'));
        load(filename, chanVar);
        temp = eval(chanVar);
        chanName{i} = chanVar;
        chanLength(i) = length(temp);
    end

    chanLength = min(chanLength);

    data = zeros(chanLength,nChan);
    for i=1:nChan
        chanVar = strcat(type, '_', num2str(channel(i), '%03d'));
        temp = eval(chanVar);
        temp = temp(1:chanLength);
        data(:,i) = double(temp');
    end
else
    chanLength = (section(2)-section(1))*fs;
    data = zeros(chanLength,nChan);
    for i=1:nChan
        chanVar = strcat(type, '_', num2str(channel(i), '%03d'));
        load(filename, chanVar);
        temp = eval(chanVar);
        temp = temp(section(1)*fs+1:section(2)*fs);
        chanName{i} = chanVar;
        data(:,i) = double(temp');
    end
end

hdr.fs = fs;
hdr.nchan = nChan;
hdr.chans = chanName;
hdr.event = [];
hdr.note = [];

end

