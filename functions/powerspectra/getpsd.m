function [psd, freq] = getpsd(signal, fs, resl, foi, normalize)
%% GETPSD Calculate PSD.
% Use as:
%   [pwr, freq] = getpsd(signal, fs, resl, foi, normalize);
% Input:
%   - signal, signal vector or signal matrix (each column as one subject)
%   - fs, the sampling rate of the signal file
%   - resl, the frequency resolution of the output PSD
%   - foi, frequency band of interest, i.e. [2 90]
%	- normalize, set as 'yes' or 'no', normalization will be done by
%	deviding the intergrated power of foi
% Output:
%   - psd, calculated PSD
%   - freq, the frequecny index of PSD
%
% Author: NIE Yingnan
% Mar 20th, 2020

win = ceil(fs/resl);
overlap = ceil(0.75*win);

[psd,freq]=pwelch(signal,win,overlap,[],fs);

psd = psd(freq>=foi(1)&freq<=foi(2),:);
freq = freq(freq>=foi(1)&freq<=foi(2),:);

if isequal(normalize,'yes')
    nSubj = size(psd,2);
    for i=1:nSubj
        psd(:,i) = psd(:,i)./(sum(psd(:,i))*resl);
    end
end

end
