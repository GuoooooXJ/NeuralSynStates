function out = nsmeasures(states, fs, ss, measure, varargin)
%NSMEASURES Calculate the neural state measures.
%   Use as:
%       out = nsmeasures(states, fs, ss, measure);
%   Input:
%       states, state vector of 1 or 0 (is or isnot the state)
%       fs, the frequency of sampling
%       ss, the state to be measured
%       measure, 'rate', 'duration', or 'occurrence'
%   Output:
%       out, the measure value
%
%   Author   : NIE Yingnan
%   Created  : Mar 6th, 2020
%   Modified : July 28th, 2020

state = (states==ss);

% Detect edge by minus shifted array, onset point will be 1 and
% offset will be -1.
temp1 = [state,0];
temp2 = [0,state];
edge = temp1-temp2;
edge(end) = [];

% Onset & offset time
onset = find(edge==1)/fs;
offset = find(edge==-1)/fs;

% Remove the last incomplete state
if ~isequal(size(onset),size(offset))
    onset(end) = [];
end

% Remove short 1 state
duration = offset-onset;
idx = duration<0.1;
onset(idx) = [];
offset(idx) = [];

switch measure
    case 'rate'
        % Calculate occurrence rate (1/min)
        out = length(onset)/(length(state)/fs)*60;
    case 'duration'
        % Calculate duration
        duration = offset-onset;
        out = mean(duration);
    case 'occurrence'
        % Calculate occurrence
        duration = offset-onset;
        out = sum(duration)/(length(state)/fs)*100;
    case 'longns'
        duration = offset-onset;
        
    otherwise
        help nsmeasures;
end



end

