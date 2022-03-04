function [pwr, freq] = getpower(signal, fs, resl, foi, normalize)
%% GETPSD Calculate PSD.
% Use as:
%   [pwr, freq] = getpower(signal, fs, resl, foi, normalize);
% Input:
%   - signal, signal vector or signal matrix (each column as one subject)
%   - fs, the sampling rate of the signal file
%   - resl, the frequency resolution of the output PSD
%   - foi, frequency band of interest, i.e. [2 90]
%	- normalize, set as 'yes' or 'no', normalization will be done by
%	deviding the intergrated power of foi
% Output:
%   - pwr, the estimate of power for each frequency
%   - freq, the frequecny index of PSD
%
% Author   : NIE Yingnan
% Created  : Aug.7, 2020
% Modified : Aug.7, 2020

win = ceil(fs/resl);
overlap = ceil(0.75*win);

[pwr,freq]=pwelch(signal,win,overlap,[],fs,'power');

pwr = pwr(freq>=foi(1)&freq<=foi(2),:);
freq = freq(freq>=foi(1)&freq<=foi(2),:);

if isequal(normalize,'yes')
    nSubj = size(pwr,2);
    for i=1:nSubj
        pwr(:,i) = pwr(:,i)./(sum(pwr(:,i))*resl);
    end
end

end
