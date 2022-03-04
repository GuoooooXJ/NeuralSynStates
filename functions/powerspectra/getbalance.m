function [balance, fidx] = getbalance(pwr, freq)
%% GETBALANCE Calculate balance matrix.
% Use as:
%   [balance, fidx] = getbalance(pwr, freq);
% Input:
%   - pwr, power matrix, each column is a power vector
%   - freq, frequency index of the power matrix
% Output:
%   - balance, balance matrix (3D)
%   - fidx, frequency index of output
%
% Author: NIE Yingnan
% Mar 31th, 2020

[nFreq,nSubj] = size(pwr);

win = 3; % Unit: points
nWin = nFreq-win+1;

intpwr = zeros(nWin, nSubj); % Integrated power
fidx = zeros(nWin,1); % Frequency index

for i = 1:nWin
    idx = i:i+win-1;
    
    % Integrated power
    intpwr(i,:) = sum(pwr(idx,:)).*(freq(2)-freq(1));
    fidx(i) = mean(freq(idx));
end

balance = zeros(nWin,nWin,nSubj);

for i = 1:nWin
    for j = i+1:nWin
        % Calculate balance
        balance(i,j,:) = intpwr(j,:)./intpwr(i,:);
    end
end

end
