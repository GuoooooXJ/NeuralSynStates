function rate = nsswitchrate(states, fs, direction, S1, S2)
%NSMEASURES Calculate the rate that NS switches between one state and
%another.
%   Use as:
%       rate = nsswitchrate(states, fs, direction, S1, S2)
%   Input:
%       states, state vector
%       fs, the frequency of sampling
%       direction, '>' - (S1 -> S2), '<' - (S1 <- S2), or '<>' - (S1 <-> S2)
%       S1, state 1
%       S2, state 2
%
%   Output:
%       rate, switch rate (1/min)
%
%   Author   : NIE Yingnan
%   Created  : July 28th, 2020
%   Modified : July 28th, 2020

time = length(states)/fs;

switch direction
    case '>'
        idxS1 = find(states==S1);
        if idxS1(end)>=length(states)
            idxS1(end) = [];
        end
        nextState = states(idxS1+1);
        nextStateS2 = (nextState==S2);
        nS2 = sum(nextStateS2);
        rate = nS2/time*60;
        
    case '<'
        rate = nsswitchrate(states, fs, '>', S2, S1);
    case '<>'
        rate1 = nsswitchrate(states, fs, '>', S1, S2);
        rate2 = nsswitchrate(states, fs, '<', S1, S2);
        rate = rate1+rate2;
end

end

