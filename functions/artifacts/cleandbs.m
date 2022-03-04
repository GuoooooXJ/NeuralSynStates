function signal = cleandbs(raw, fsraw, fs, fdbs, width, initial, varargin)
%CLEANDBS Remove DBS artifacts using resampling method.
%   Use as:
%       signal = cleandbs(raw, fsraw, fs, fdbs, width, initial);
%       signal = cleandbs(raw, fsraw, fs, fdbs, width, initial,thr);
%   Input:
%       - raw, raw signal with DBS artifacts
%       - fsraw, sampling rate of raw signal
%       - fs, sampling rate of output signal
%       - fdbs, the frequency of DBS stimulation
%       - width, the width of the artifact (us) (usually 2*pulsewidth)
%       - initial, DBS form 1 or -1,     __                   __
%                                    ___|  |   ___  or ___   |  |___
%                                          |__|           |__|
%       - thr, threshold (optional)
%   Output:
%       - signal, signal without DBS artifacts
%       - fs, sampling rate of output signal
%
%   Author   : NIE Yingnan
%   Created  : Aug 18, 2020
%   Modified : Aug 22, 2021

if fsraw<20*10^3 % 20k
    error("The sampling rate of raw signal must bigger than 20kHz.");
end

% if 1/fs*10^6<(width/2)
%     fprintf("The max output sampling rate for your signal is %d.\n", floor(2/(width/10^6)));
%     error("Invalid input.");
% end

raw = raw*sign(initial);

% Detrend
raw = detrend(raw);

% Threshloding
if nargin==6
    threshold = dbsthreshold(raw,fsraw,fdbs);
elseif nargin==7
    threshold = varargin{1};
else
    error("Too many input parameters.");
end



time = (1/fsraw:1/fsraw:length(raw)/fsraw)';

% Find artifacts
lenSignal = length(raw);
idx = 1;
nPeaks = 0;
gap = floor(0.75/fdbs*fsraw); % Skip 3/4 interval of the DBS pulses
while idx<lenSignal
    if raw(idx)>threshold
        nPeaks = nPeaks+1;
        peaks(nPeaks) = idx;
        idx = idx+gap;
    end
    idx = idx+1;
end


% Remove artifacts
usPreCut =500; %us
usInterval = 1/fdbs*10^6-usPreCut-width; % us

pPrecut = ceil(fsraw*usPreCut/10^6);
pInterval = ceil(fsraw*usInterval/10^6);

sampleVector = [];
timeVector = [];


for i=2:nPeaks % Ignore the first & one
    peakIdx = peaks(i);
    startIdx = peakIdx-pPrecut-pInterval+1;
    endIdx = peakIdx-pPrecut;
    
    sampleVector = [sampleVector;raw(startIdx:endIdx)];
    timeVector = [timeVector;time(startIdx:endIdx)];
end

% Resampling
[signal, ty] = resample(sampleVector,timeVector,fs,'pchip');

signal = signal*sign(initial);

signal = detrend(signal);

figure;
subplot(2,1,1);
hold on;
plot(time, raw);
plot(time(peaks), raw(peaks), 'o');
hold off;

subplot(2,1,2);
plot(ty, signal);

end

