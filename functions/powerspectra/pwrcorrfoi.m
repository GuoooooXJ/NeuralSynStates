function [r, p] = pwrcorrfoi(varargin)
%% PWRCORRFOI Power spectrum analysis within FOI.
% Use as:
%   [r, p] = pwrcorr(pwr, freq, score, foi, file);
%   Correlation analysis between each frequency band and score.
% Input:
%   - pwr, power matrix, each column is a power vector
%   - freq, frequency index of the power matrix
%   - score, a column vector of assessment scores
%   - foi, frequency band of interest
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - r, the correlation coefficients of all frequency bands
%   - p, the p-value of all frequency bands
%
% Also use as:
%   [r, p] = pwrcorr(pwr1, pwr2, freq, score1, score2, foi, file);
%   Correlation between the power difference (pwr 1-2) and symptom improvements
%   (score 2-1).
% Input:
%   - pwr1, power matrix of state 1, each column is a power vector
%   - pwr2, power matrix of state 2, each column is a power vector
%   - freq, frequency index of the power matrix
%   - score1, a column vector of assessment scores of state 1
%   - score2, a column vector of assessment scores of state 1
%   - foi, frequency band of interest
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - r, the correlation coefficients of all frequency bands
%   - p, the p-value of all frequency bands

% Author: Nie Yingnan
% Mar 31th, 2020

switch nargin
    case 5
        pwr = varargin{1};
        freq = varargin{2};
        score = varargin{3};
        foi = varargin{4};
        file = varargin{5};
                
        idx = freq>foi(1)&freq<foi(2);
        intpwr(:) = sum(pwr(idx,:)).*(freq(2)-freq(1));
        intpwr = intpwr';

        % Correlation between power and score
        [hnorm1,~] = lillietest(intpwr,'alpha',0.05);
        [hnorm2,~] = lillietest(score,'alpha',0.05);
        
        if hnorm1||hnorm2
            % Spearman
            method = 'Spearman';
            [r, p] = corr(intpwr,score,'type','Spearman');
        else
            % Pearson
            method = 'Pearson';
            [r, p] = corr(intpwr,score,'type','Pearson');
        end
        
        fprintf("R: %f, p: %f\n",r,p);
        if p<0.05
            fig = figure;
            scatter(intpwr,score);
            if ~isequal(file,'no')
                print(fig,'-djpeg','-r300', strcat(file,'pwrcorr.jpeg'));
            end
        end
        
        if ~isequal(file,'no')
            save(strcat(file, 'pwrcorr.mat'),'intpwr', 'p', 'r','method');
        end
 
    case 7
        pwr1 = varargin{1};
        pwr2 = varargin{2};
        freq = varargin{3};
        score1 = varargin{4};
        score2 = varargin{5};
        foi = varargin{6};
        file = varargin{7};
        
        if ~isequal(size(pwr1),size(pwr2))
            disp("The dimension of pwr1 and pwr2 not match !");
            return
        end
        
        if ~isequal(size(score1),size(score2))
            disp("The dimension of score1 and score2 not match !");
            return
        end
        
        idx = freq>foi(1)&freq<foi(2);
        intpwr1(:) = sum(pwr1(idx,:)).*(freq(2)-freq(1));
        intpwr2(:) = sum(pwr2(idx,:)).*(freq(2)-freq(1));
        intpwr1 = intpwr1';
        intpwr2 = intpwr2';

        % Correlation between power and score
        pwrdiff = intpwr2-intpwr1; % Power change
        scorediff = score1-score2; % Imporvement
        
        [hnorm1,~] = lillietest(pwrdiff,'alpha',0.05);
        [hnorm2,~] = lillietest(scorediff,'alpha',0.05);
        
        if hnorm1||hnorm2
            % Spearman
            method = 'Spearman';
            [r, p] = corr(pwrdiff,scorediff,'type','Spearman');
        else
            % Pearson
            method = 'Pearson';
            [r, p] = corr(pwrdiff,scorediff,'type','Pearson');
        end
        
        fprintf("R: %f, p: %f\n",r,p);
        if p<0.05
            fig = figure;
            scatter(pwrdiff,scorediff);
            if ~isequal(file,'no')
                print(fig,'-djpeg','-r300', strcat(file,'pwrdiffcorr.jpeg'));
            end
        end
                
        if ~isequal(file,'no')
            save(strcat(file, 'pwrdiffcorr.mat'),'pwrdiff', 'p', 'r','method');
        end
    otherwise
        help pwrcorrfoi;
end
            
end
