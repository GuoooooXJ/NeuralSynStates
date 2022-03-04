function [h,p] = balancetest(balance1, balance2, freq, thr, file)
%% BALANCETEST Paired test if two group's balance are different.
% Use as:
%   [h, p] = balancetest(pwr1, pwr2, freq, thr, file);
% Input:
%   - balance1, balance matrix
%   - balance2, balance matrix
%   - freq, frequency index of the power matrix
%	- file, set as 'no'(not save) or 'filename' to save the results
% Output:
%   - h, the hypothesis of all frequency bands
%   - p, the p-value of all frequency bands
%
% Author: NIE Yingnan
% Mar 31th, 2020

if ~isequal(size(balance1),size(balance2))
    printf("The size of balance1&balance2 not match!\n");
    return
end

nFreq = length(freq);

hnorm1 = zeros(nFreq,nFreq); % Normality test
hnorm2 = zeros(nFreq,nFreq);
for i = 1:nFreq
    for j = i+1:nFreq
        % Normality test
        tmp1(:) = balance1(i,j,:);
        tmp2(:) = balance2(i,j,:);
        [hnorm1(i),~] = lillietest(tmp1,'alpha',thr);
        [hnorm2(i),~] = lillietest(tmp2,'alpha',thr);
    end
end

p = ones(nFreq,nFreq); % p-value of each frequency
h = zeros(nFreq,nFreq); % Significant different
if any(hnorm1)|any(hnorm2)
    % Wilcoxon signed rank test
    testmethod = 'Wilconxon';
    for i = 1:nFreq
        for j = i+1:nFreq
            tmp1(:) = balance1(i,j,:);
            tmp2(:) = balance2(i,j,:);
            [p(i,j),h(i,j)] = signrank(tmp1,tmp2,'alpha',thr);
        end
    end
else
    % Paired t-test
    testmethod = 'ttest';
    for i = 1:nFreq
        for j = i+1:nFreq
            tmp1(:) = balance1(i,j,:);
            tmp2(:) = balance2(i,j,:);
            [h(i,j),p(i,j)] = ttest(tmp1,tmp2,'tail','both','alpha',thr);
        end
    end
end

meanbalance1 = mean(balance1,3);
meanbalance2 = mean(balance2,3);
diffbalance = meanbalance2-meanbalance1;

if ~isequal(file,'no')
    save(strcat(file, 'balancetest.mat'),'meanbalance1','meanbalance2','diffbalance', 'p', 'h','testmethod');
end


% Plot
fig = figure;
imagesc([freq(1), freq(end)],[freq(1), freq(end)],meanbalance1);
colormap(jet);
colorbar;
axis('square');
set(gca,'Box','on');
title('Mean Balance 1','Fontsize',20);
xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
set(gca,'FontSize',20,'fontname','Arial');
if ~isequal(file,'no')
    print(fig,'-djpeg','-r300', strcat(file,'balancetest_mean1.jpeg'));
end

fig = figure;
imagesc([freq(1), freq(end)],[freq(1), freq(end)],meanbalance2);
colormap(jet);
colorbar;
axis('square');
set(gca,'Box','on');
title('Mean Balance 2','Fontsize',20);
xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
set(gca,'FontSize',20,'fontname','Arial');
if ~isequal(file,'no')
    print(fig,'-djpeg','-r300', strcat(file,'balancetest_mean2.jpeg'));
end

fig = figure;
imagesc([freq(1), freq(end)],[freq(1), freq(end)],diffbalance);
colormap(jet);
colorbar;
axis('square');
set(gca,'Box','on');
title('Balance Difference 2-1','Fontsize',20);
xlabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
ylabel('Frequency (Hz)','fontsize',20,'fontname','Arial');
set(gca,'FontSize',20,'fontname','Arial');
if ~isequal(file,'no')
    print(fig,'-djpeg','-r300', strcat(file,'balancetest_diff2-1.jpeg'));
end

p05 = p<0.05;
p01 = p<0.01;
hotmap = p05+p01;

color = [0,0,0.515625;1,1,0;0.5,0,0];

fig = figure;
imagesc([freq(1), freq(end)],[freq(1), freq(end)],hotmap);
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
    print(fig,'-djpeg','-r300', strcat(file,'balancetest_p.jpeg'));
end


end
