function thr = dbsthreshold(raw, fs, freq)
%DBSTHRESHOLD Get the threshold for DBS artifact detection.
%The amplitude of postive and negtive pulses must be the same:
%    __
%___|  |   ___
%      |__|
%
%   Use as:
%       thr = dbsthreshold(raw, fs, freq);
%   Input:
%       - raw, raw signal with DBS artifacts
%       - fs, sampling rate
%   Output:
%       - thr, threshold for DBS artifact detection

time = length(raw)/fs;
nPulse = time*freq; % The theoretical number of pulses.

raw(raw<0) = 0; % Only consider positive pulses

pks = findpeaks(raw);
pks = maxk(pks,ceil(nPulse*0.8));

ampl = mean(pks);

thr = 0.7*ampl;

end

