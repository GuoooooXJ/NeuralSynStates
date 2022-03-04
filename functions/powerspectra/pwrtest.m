function [f, p] = pwrtest(pwr1, pwr2, freq, thr, file)
%% PWRTEST Paired test if two group's power are different.
% Use as:
%   [f, p] = pwrtest(pwr1, pwr2, freq, thr, file);
%   Correlation analysis between each frequency band and score.
% Input:
%   - pwr1, power matrix, each column is a power vector
%   - pwr2, power matrix, each column is a power vector
%   - freq, frequency index of the power matrix
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - f, frequency index of output
%   - p, the p-value of all frequency bands
%
% Author: NIE Yingnan
% Mar 23th, 2020

if ~isequal(size(pwr1),size(pwr2))
    printf("The size of pwr1&pwr2 not match!\n");
    return
end

[nFreq,nSubj] = size(pwr1);

win = 3; % Unit: points
nWin = nFreq-win+1;

intpwr1 = zeros(nWin, nSubj); % Integrated power
intpwr2 = zeros(nWin, nSubj);
f = zeros(nWin,1); % Frequency index
p = zeros(nWin,1); % p-value of each frequency

hnorm1 = zeros(nWin,1); % Normality test
hnorm2 = zeros(nWin,1);

for i = 1:nWin
    idx = i:i+win-1;
    
    % Integrated power
    intpwr1(i,:) = sum(pwr1(idx,:)).*(freq(2)-freq(1));
    intpwr2(i,:) = sum(pwr2(idx,:)).*(freq(2)-freq(1));
    f(i) = mean(freq(idx));
    
    % Normality test
    [hnorm1(i),~] = lillietest(intpwr1(i,:),'alpha',thr);
    [hnorm2(i),~] = lillietest(intpwr2(i,:),'alpha',thr);
end

% if any(hnorm1)||any(hnorm2)
    % Wilcoxon signed rank test
    testmethod = 'Wilconxon';
    for i = 1:nWin
        [p(i),~] = signrank(intpwr1(i,:),intpwr2(i,:),'alpha',thr);
    end
% else
%     % Paired t-test
%     testmethod = 'ttest';
%     for i = 1:nWin
%         [~,p(i)] = ttest(intpwr1(i,:),intpwr2(i,:),'tail','both','alpha',thr);
%     end
% end

meanpwr1 = mean(intpwr1,2);
meanpwr2 = mean(intpwr2,2);

% Plot
fig = figure;
set(get(gca, 'Title'), 'String', testmethod);

hold on;
% Plot lines to auto set the axes
plot(f, meanpwr1,'LineWidth',2,'Color','#FE2E64');
plot(f, meanpwr2,'LineWidth',2,'Color','#0080FF');
ylim = get(gca,'YLim'); 
fig1 = area(f,(p<thr)*ylim,'FaceColor','#F3F781','EdgeColor','#F3F781');
% Plot again to make them visible
fig2 = plot(f, meanpwr1,'LineWidth',2,'Color','#FE2E64');
fig3 = plot(f, meanpwr2,'LineWidth',2,'Color','#0080FF');
legend([fig1(1);fig2(1);fig3(1)],'Significant','State1','State2');
set(gca,'Layer','top');
hold off;

if ~isequal(file,'no')
    print(fig,'-djpeg','-r300',strcat(file, 'pwrtest.jpeg'));
    save(strcat(file, 'pwrtest.mat'),'meanpwr1', 'meanpwr2','f','p','testmethod');
end

end
