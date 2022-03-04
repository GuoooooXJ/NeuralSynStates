function [r, p] = pwrcorr(varargin)
%% PWRCORR Power spectrum analysis.
% Use as:
%   [r, p] = pwrcorr(pwr, freq, score, file);
%   Correlation analysis between each frequency band and score.
% Input:
%   - pwr, power matrix, each column is a power vector
%   - freq, frequency index of the power matrix
%   - score, a column vector of assessment scores
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - r, the correlation coefficients of all frequency bands
%   - p, the p-value of all frequency bands
%
% Also use as:
%   [r, p] = pwrcorr(pwr1, pwr2, freq, score1, score2, file);
%   Correlation between the power difference (pwr 1-2) and symptom improvements
%   (score 1-2).
% Input:
%   - pwr1, power matrix of state 1, each column is a power vector
%   - pwr2, power matrix of state 2, each column is a power vector
%   - freq, frequency index of the power matrix
%   - score1, a column vector of assessment scores of state 1
%   - score2, a column vector of assessment scores of state 1
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - r, the correlation coefficients of all frequency bands
%   - p, the p-value of all frequency bands

% Author: Nie Yingnan
% Mar 21th, 2020

switch nargin
    case 4
        pwr = varargin{1};
        freq = varargin{2};
        score = varargin{3};
        file = varargin{4};
        
        [nFreq,nSubj] = size(pwr);

        % Integrated power
        win = 3; % Unit: points
        nWin = nFreq-win+1;

        intpwr = zeros(nWin, nSubj); % Integrated power
        intfreq = zeros(nWin,1); % Frequency index

        for i = 1:nWin
            idx = i:i+win-1;
            intpwr(i,:) = sum(pwr(idx,:)).*(freq(2)-freq(1));
            intfreq(i) = mean(freq(idx));
        end

        % Correlation between power and score
        [r,p] = corr(intpwr', score,'type','spearman');

        % Plot
        fig = figure;
        set(get(gca, 'Title'), 'String', 'Correlation between power and score');
        hold on;
        plot(intfreq, r);
        plot(intfreq, p);
        hold off;
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300',strcat(file, '-pwstat.jpeg'));
            save(strcat(file, '-pwrcorr.mat'),'r','p');
        end
    case 6
        pwr1 = varargin{1};
        pwr2 = varargin{2};
        freq = varargin{3};
        score1 = varargin{4};
        score2 = varargin{5};
        file = varargin{6};
        
        if ~isequal(size(pwr1),size(pwr2))
            disp("The dimension of pwr1 and pwr2 not match !");
            return
        end
        
        if ~isequal(size(score1),size(score2))
            disp("The dimension of score1 and score2 not match !");
            return
        end
        
        [nFreq,nSubj] = size(pwr1);

        % Integrated power
        win = 3; % Unit: points
        nWin = nFreq-win+1;

        intpwr1 = zeros(nWin, nSubj); % Integrated power
        intpwr2 = zeros(nWin, nSubj);
        intfreq = zeros(nWin,1); % Frequency index

        for i = 1:nWin
            idx = i:i+win-1;
            % Integrated power
            intpwr1(i,:) = sum(pwr1(idx,:)).*(freq(2)-freq(1));
            intpwr2(i,:) = sum(pwr2(idx,:)).*(freq(2)-freq(1));
            intfreq(i) = mean(freq(idx));
        end
        
        % Calculate difference
        diffpwr = intpwr1-intpwr2; % Power difference between two states
        diffscore = score1-score2; % Score difference between two states
        
        % Correlation between power difference and score difference
        [r,p] = corr(diffpwr', diffscore, 'type', 'spearman');
        
        % Plot
        fig = figure;
        set(get(gca, 'Title'), 'String', 'Correlation between difference of power and score');
        hold on;
        plot(intfreq, r);
        plot(intfreq, p);
        hold off;
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300',strcat(file, '-pwrdiffcorr.jpeg'));
            save(strcat(file, '-pwrdiffcorr.mat'),'r','p');
        end
    otherwise
        help pwrcorr;
end
            
end
