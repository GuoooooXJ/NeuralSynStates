function count = countdurations(states, fs, ss, segmt)
%COUNTDURATIONS Count the neural state duration.
%   Use as:
%       durations = countdurations(states, fs, ss);
%   Input:
%       states, state vector of 1 or 0 (is or isnot the state)
%       fs, the frequency of sampling
%       ss, the state to be measured
%       segmt, segment lim [0,   0.1;
%                           0.1, 0.2;
%                           ...]
%   Output:
%       durations, vector of duraions
%
%   Author   : NIE Yingnan
%   Created  : Dec 14th, 2020
%   Modified : Dec 14th, 2020

% Calculate duration
durations = getdurations(states, fs, ss);

% Count
N = length(durations);
period = length(states)/fs;

count = zeros(1,size(segmt,1));
if N~=0
    for i=1:size(segmt,1)
        lLim = segmt(i,1);
        hLim = segmt(i,2);

        idx = durations>lLim&durations<=hLim;

        count(i) = sum(idx)/period*60; % Times per min.
    end
end

end

