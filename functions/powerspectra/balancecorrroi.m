function [r, p] = balancecorrroi(varargin)
%% BALANCECORRROI Balance (power ratio) analysis within a ROI.
% Use as:
%   [r, p] = balancecorrroi(balance, score, roi, file);
% Input:
%   - balance, balance matrix (3D)
%   - score, a column vector of assessment scores
%   - roi, ROI matrix
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - r, the correlation coefficients of all frequency bands
%   - p, the p-value of all frequency bands
%
% Also use as:
%   [r, p] = balancecorrroi(balance1, balance2, score1, score2, roi, file);
%   Correlation between the balance difference (balance 1-2) and symptom 
%   improvements(score 1-2).
% Input:
%   - balance1, power ratio matrix of state 1
%   - balance2, power ratio matrix of state 2
%   - score1, a column vector of assessment scores of state 1
%   - score2, a column vector of assessment scores of state 1
%   - roi, ROI matrix
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - r, the correlation coefficients of all frequency bands
%   - p, the p-value of all frequency bands

% Author: Nie Yingnan
% Mar 31th, 2020

switch nargin
    case 4
        balance = varargin{1};
        score = varargin{2};
        roi = varargin{3};
        file = varargin{4};
        
        nSubj = size(balance,3);
             
        % Correlation between power ratio and score
        [m,n] = find(roi==1);
        tmp = zeros(length(m),nSubj); % Balance matrix of ROI
        
        for i=1:length(m)    
        	tmp(i,:) = balance(m(i),n(i),:); 
        end
        
        % Mean power ratio
        meanbalance = mean(tmp,1)';
        
        % Normality test
        [hnorm1,~] = lillietest(meanbalance,'alpha',0.05);
        [hnorm2,~] = lillietest(score,'alpha',0.05);

        if hnorm1||hnorm2
            % Spearman
            method = 'Spearman';
            [r, p] = corr(meanbalance,score,'type','Spearman');
        else
            % Pearson
            method = 'Pearson';
            [r, p] = corr(meanbalance,score,'type','Pearson');
        end
        
        fprintf("R: %f, p: %f\n",r,p);
        if p<0.05
            fig = figure;
            scatter(meanbalance,score);
            if ~isequal(file,'no')
                print(fig,'-djpeg','-r300', strcat(file,'balancecorrroi.jpeg'));
            end
        end
        
        if ~isequal(file,'no')
            save(strcat(file, 'balancecorrroi.mat'),'meanbalance', 'p', 'r','method');
        end

    case 6
        balance1 = varargin{1};
        balance2 = varargin{2};
        score1 = varargin{3};
        score2 = varargin{4};
        roi = varargin{5};
        file = varargin{6};
        
        if ~isequal(size(balance1),size(balance2))
            disp("The dimension of pwr1 and pwr2 not match !");
            return
        end
        
        if ~isequal(size(score1),size(score2))
            disp("The dimension of score1 and score2 not match !");
            return
        end
        
        nSubj = size(balance1,3);

        [m,n] = find(roi==1);
        tmp1 = zeros(length(m),nSubj); % Balance matrix of ROI
        tmp2 = zeros(length(m),nSubj);
        for i=1:length(m)    
        	tmp1(i,:) = balance1(m(i),n(i),:); 
            tmp2(i,:) = balance2(m(i),n(i),:); 
        end
        
        meanbalance1 = mean(tmp1,1)';
        meanbalance2 = mean(tmp2,1)';
        
        % Correlation between difference of power ratio and score
        balancediff = meanbalance2-meanbalance1; % Change
        scorediff = score1-score2; % Improvement

        % Normality test
        [hnorm1,~] = lillietest(balancediff,'alpha',0.05);
        [hnorm2,~] = lillietest(scorediff,'alpha',0.05);

        if hnorm1||hnorm2
            % Spearman
            method = 'Spearman';
            [r, p] = corr(balancediff,scorediff,'type','Spearman');
        else
            % Pearson
            method = 'Pearson';
            [r, p] = corr(balancediff,scorediff,'type','Pearson');
        end
        
        fprintf("R: %f, p: %f\n",r,p);
        if p<0.05
            fig = figure;
            scatter(balancediff,scorediff);
            if ~isequal(file,'no')
                print(fig,'-djpeg','-r300', strcat(file,'balancediffcorrroi.jpeg'));
            end
        end
        
        if ~isequal(file,'no')
            save(strcat(file, 'balancediffcorrroi.mat'),'meanbalance1','meanbalance2', 'p', 'r','method');
        end

    otherwise
        help balancecorrroi;
end

end
