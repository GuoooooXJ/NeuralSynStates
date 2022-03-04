function durations = getdurations(states, fs, ss)
%GETDURATIONS Calculate the neural state measures.
%   Use as:
%       durations = getdurations(states, fs, ss1);
%   Input:
%       states, state vector of 1 or 0 (is or isnot the state)
%       fs, the frequency of sampling
%       ss, the state to be measured
%   Output:
%       durations, vector of duraions
%
%   Author   : NIE Yingnan
%   Created  : Dec 14th, 2020
%   Modified : Dec 14th, 2020

state = (states==ss);

% Detect edge by minus shifted array, onset point will be 1 and
% offset will be -1.
temp1 = [state;0];
temp2 = [0;state];
edge = temp1-temp2;
edge(end) = [];

% Onset & offset time
onset = find(edge==1)/fs;
offset = find(edge==-1)/fs;

% Remove the last incomplete state
if ~isequal(size(onset),size(offset))
    onset(end) = [];
end

% % Remove short 1 state
% duration = offset-onset;
% idx = duration<0.1;
% onset(idx) = [];
% offset(idx) = [];

% Calculate duration
durations = offset-onset;

end

