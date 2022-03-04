function [h, p]= pwrtestfoi(pwr1, pwr2, freq, foi, varargin)
%% PWRTESTROI Paired test if two group's power (from a given foi) 
% are different.
% Use as:
%   [h, p] = pwrtestroi(pwr1, pwr2, freq, foi, options);
%   Correlation analysis between each frequency band and score.
% Input:
%   - pwr1, power matrix, each column is a power vector
%   - pwr2, power matrix, each column is a power vector
%   - freq, frequency index of the power matrix
%	- foi, frequency of interest
%   options:
%       '-thr', threshold for the hypothesis test, default 0.05
%       '-method', 'ttest' for paired t-test, 'Wilconxon' for Wilconxon 
%       signed rank test, default ttest
% Output:
%   - h, testing result, 1 - is different, 0 - not different
%   - p, the p-value of the hypothesis test
%
% Author   : NIE Yingnan
% Created  : Aug 8, 2020
% Modified : Aug 8, 2020

thr = 0.05;
method = 'ttest';

if nargin>4
    for i=1:2:nargin-4
        option = varargin{i};
        switch option
            case '-thr'
                thr = varargin{i+1};
            case '-method'
                method = varargin{i+1};
        end
    end
end

if ~isequal(size(pwr1),size(pwr2))
    error("The size of pwr1&pwr2 not match!\n");
end

idx = (freq>foi(1))&(freq<foi(2));

% Integrated power
intpwr1 = sum(pwr1(idx,:)).*(freq(2)-freq(1));
intpwr2 = sum(pwr2(idx,:)).*(freq(2)-freq(1));

switch method
    case 'ttest'
        [h,p] = ttest(intpwr1,intpwr2,'tail','both','alpha',thr);
    case 'Wilconxon'
        [p,h] = signrank(intpwr1,intpwr2,'alpha',thr);
end

end
