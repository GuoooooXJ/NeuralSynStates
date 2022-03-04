function mapfile2nspl(filenameIn, filenameOut, type, chan, time)
%MAPFILE2NSPL Convert Alpha Omega mapfile to nspl txt file
%  Use as:
%    mapfile2nspl(filenameIn, filenameOut, type, chan);
%  Input:
%    - filenameIn, input file name (mapfile)
%    - filenameOut, output file name (nspl)
%    - type, channel type, 'LFP','RAW'
%    - chan, channel list, e.g. [1,4,7]
%    - time, optional, select a time range [0,50], unit second
%
%  Author   : NIE Yingnan
%  Created  : Nov.19, 2019
%  Modified : Aug 18, 2020

type = strcat('C',type);

fsVar = strcat(type, '_', num2str(chan(1), '%03d'),'_KHz');
load(filenameIn, fsVar);
fs = eval(fsVar)*1000;

chanLength = zeros(length(chan),1); % The length of channels maybe different.
for i=1:length(chan)
    chanVar = strcat(type, '_', num2str(chan(i), '%03d'));
    load(filenameIn, chanVar);
    temp = eval(chanVar);
    chanLength(i) = length(temp);
end

chanLength = min(chanLength);

data = [];
for i=1:length(chan)
    chanVar = strcat(type, '_', num2str(chan(i), '%03d'));
    temp = eval(chanVar);
    temp = temp(1:chanLength);
    data = [data;temp];
end

data = double(data');

if nargin==5
    data = data(time(1)*fs+1:time(2)*fs,:);
end

savenspltxt(filenameOut, data, fs);

end

