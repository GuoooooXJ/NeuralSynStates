function [r, p] = balancecorr(varargin)
%% BALANCECORR Balance (power ratio) analysis.
% Use as:
%   [r, p] = balancecorr(pwr, freq, score, file);
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
%   [r, p] = balancecorr(pwr1, pwr2, freq, score1, score2, file);
%   Correlation between the balance difference (balance 1-2) and symptom 
%   improvements(score 1-2).
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
        
        % Correlation between balance and score
        r = zeros(nWin,nWin);
        p = ones(nWin,nWin);
        balance = zeros(nWin,nWin,nSubj); % Balance
        for i=1:nWin
            for j=i:nWin
                tmp = intpwr(j,:)./intpwr(i,:);
                balance(i,j,:) = tmp; 
                [r(i,j), p(i,j)]=corr(tmp', score,'type','spearman');
                % 正态检验 spearman
                % TODO
                
            end
        end

        % Mean balance
        meanpr = mean(balance,3);

        if ~isequal(file,'no')
            save(strcat(file, '-balancecorr.mat'),'r','p','meanpr');
        end

        % Plot
        fig = figure;
        imagesc([intfreq(1), intfreq(end)],[intfreq(1), intfreq(end)],meanpr);
        colormap(jet);
        colorbar;
        axis('square');
        set(gca,'Box','on');
        title('Mean Power Ratio','Fontsize',20);
        xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        set(gca,'FontSize',20,'fontname','Arial');
        set(gca, 'xLim',[intfreq(1), intfreq(end)]); 
        set(gca, 'yLim',[intfreq(1), intfreq(end)]); 
        set(gca, 'xTick',intfreq(1):10:intfreq(end));
        set(gca, 'xTicklabel',intfreq(1):10:intfreq(end));
        set(gca, 'yTick',intfreq(1):10:intfreq(end));
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300', strcat(file,'-balancecorr-meanpr.jpeg'));
        end

        fig = figure;
        imagesc([intfreq(1), intfreq(end)],[intfreq(1), intfreq(end)],r);
        colormap(jet);
        caxis([-1,1]);
        colorbar;
        axis('square');
        set(gca,'Box','on');
        title('Corrlation Coefficients','Fontsize',20);
        xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        set(gca,'FontSize',20,'fontname','Arial');
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300', strcat(file,'-balancecorr-r.jpeg'));
        end

        p05 = p<0.05;
        p01 = p<0.01;
        hotmap = p05+p01;

        color = [0,0,0.515625;1,1,0;0.5,0,0];

        fig = figure;
        imagesc([intfreq(1),intfreq(end)],[intfreq(1),intfreq(end)],hotmap);
        colormap(color);
        hc = colorbar;
        set(hc,'YTick',[0,2/3,4/3,2]);
        set(hc,'YTickLabel',{'1','0.05','0.01','0'});
        axis('square');
        set(gca,'Box','on');
        title('Hot area','Fontsize',20);
        xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        set(gca,'FontSize',20,'fontname','Arial');
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300', strcat(file,'-balancecorr-p.jpeg'));
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
        
        % Correlation between difference of balance and score
        r = zeros(nWin,nWin);
        p = ones(nWin,nWin);
        balancediff = zeros(nWin,nWin,nSubj); % Balance difference
        scorediff = score1-score2;
        for i=1:nWin
            for j=i:nWin
                tmp = intpwr1(j,:)./intpwr1(i,:)-intpwr2(j,:)./intpwr2(i,:);
                balancediff(i,j,:) = tmp; 
                [r(i,j), p(i,j)]=corr(tmp', scorediff);
            end
        end

        % Mean balance
        meanpr = mean(balancediff,3);
        
        if ~isequal(file,'no')
            save(strcat(file, '-balancediffcorr.mat'),'r','p','meanpr');
        end

        % Plot
        fig = figure;
        imagesc([intfreq(1), intfreq(end)],[intfreq(1), intfreq(end)],meanpr);
        colormap(jet);
        colorbar;
        axis('square');
        set(gca,'Box','on');
        title('Mean Power Ratio','Fontsize',20);
        xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        set(gca,'FontSize',20,'fontname','Arial');
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300', strcat(file,'-balancediffcorr-meanpr.jpeg'));
        end

        fig = figure;
        imagesc([intfreq(1), intfreq(end)],[intfreq(1), intfreq(end)],r);
        colormap(jet);
        caxis([-1,1]);
        colorbar;
        axis('square');
        set(gca,'Box','on');
        title('Corrlation Coefficients','Fontsize',20);
        xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        set(gca,'FontSize',20,'fontname','Arial');
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300', strcat(file,'-balancediffcorr-r.jpeg'));
        end

        p05 = p<0.05;
        p01 = p<0.01;
        hotmap = p05+p01;

        color = [0,0,0.515625;1,1,0;0.5,0,0];

        fig = figure;
        imagesc([intfreq(1),intfreq(end)],[intfreq(1),intfreq(end)],hotmap);
        colormap(color);
        hc = colorbar;
        set(hc,'YTick',[0,2/3,4/3,2]);
        set(hc,'YTickLabel',{'1','0.05','0.01','0'});
        axis('square');
        set(gca,'Box','on');
        title('Hot area','Fontsize',20);
        xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
        set(gca,'FontSize',20,'fontname','Arial');
        if ~isequal(file,'no')
            print(fig,'-djpeg','-r300', strcat(file,'-balancediffcorr-p.jpeg'));
        end
        
    otherwise
        help balancecorr;
end

end
