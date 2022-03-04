function [varargout] = NeuralSynCWT(signal, fs,  varargin)
%Neural Syncronization States Caculation
%   Detect synchronization metric for LFP signal With SWT.
%
%   Use as:
%       [SynMetric, tVector] = dnsdetection(signal, fs, node, options, values);

%Input:
%       signal      - signal vector
%       fs          - sampling rate of input signal
%       node        - target node, should match the target frequency band
%
%   Output:
%       states      - synchronization states vector
%       tVector     - time vector that matches the states vector
%
%   Author   : Xuanjun Guo
%   Created  : Dec 16, 2021
%   Modified : Dec 31, 2022

basis = 'amor';
step  = round(0.02 * fs);

% Get customized settings
nArg = nargin - 2;
if mod(nArg, 2) ~= 0
    help dnsdetection
    return
end


% Create time vector
time = (1:length(signal))./fs; 

%{
% Init progress bar
N = floor((length(signal)-(winlen-step))/step);
cont = 0;
pro = progress('Neural State', N);
%}

% Compute CWT coefficients sequence
[wt,f] = cwt(signal,basis,fs);
figure;
pcolor(time,f,abs(wt));shading interp;
title('CWT on PAGdata');


% Remove the points before and after to exclude edge effects
WTwithoutEdge = wt(:,fs:end-fs);

[FrequencyNumber,n] = size(WTwithoutEdge);
WTwithoutEdgeNew = WTwithoutEdge(1:round(FrequencyNumber/2),:);
frequency = f(1:round(FrequencyNumber/2));
[FrequencyNumber,n] = size(WTwithoutEdgeNew);

centralFreq = round(frequency(end));
winlen = round(1/centralFreq * fs * 4);
thrswin = round(1/centralFreq * fs * 60);

% Compute wavelet pocket coefficients for priori window
nStepThrs = ceil((thrswin-(winlen-step))/step);
thrs  = zeros(FrequencyNumber, nStepThrs);
sig   = zeros(FrequencyNumber, nStepThrs);
thrsel = zeros(FrequencyNumber, nStepThrs);
%thrCal = zeros(FrequencyNumber, nStepThrs);
%thrManuls = zeros(FrequencyNumber, nStepThrs);


% Detect synchronization state
startPoint = nStepThrs*step + 1;
nStep = floor((length(WTwithoutEdgeNew(1,:))-(startPoint-1)-(winlen-step))/step);

cfs = zeros(FrequencyNumber, step*nStep);
states = zeros(FrequencyNumber, nStep);
K1 = zeros(FrequencyNumber, nStep);
Kr = zeros(FrequencyNumber, nStep);

for frequencyBand = 1:FrequencyNumber-1

centralFreq = round(f(frequencyBand));
winlenCurrent = round(1/centralFreq * fs * 4);
thrswinCurrent = round(1/centralFreq * fs * 60);
thrswinlen = thrswinCurrent/fs;
    
for iStep = 1:nStep
    
    idx = (1:step)+(iStep-1)*step;

    % Extract coefficients
    cfsTempCurrent  = WTwithoutEdgeNew(frequencyBand,(startPoint:startPoint + winlenCurrent - 1)+(iStep-1)*step);
    cfsThrsCurrent  = WTwithoutEdgeNew(frequencyBand,(startPoint-thrswinCurrent:startPoint)+(iStep-1)*step);
    cfs(frequencyBand,idx) = cfsTempCurrent(end-step+1:end);
    
    % Compute threshold
    % According to the paper, Normlization before median
        
    % Control chart
    flag = 0;
    ControlChartNumber = 0;
    while(flag == 0)
    ControlChartNumber = ControlChartNumber+1;
    meancfsThrs  = mean(cfsThrsCurrent);
    sigmacfsThrs = std(cfsThrsCurrent);
    cfsThrsNew = cfsThrsCurrent(cfsThrsCurrent<meancfsThrs+2.58*sigmacfsThrs & cfsThrsCurrent>meancfsThrs-2.58*sigmacfsThrs);
        % 1%-2.58 2%-2.33
        %{
        tVectorPlot = (1:length(cfsThrs))./fs;
        tVectorPlotNew = (1:length(cfsThrsNew))./fs;
        figure;
        plot(tVectorPlot,cfsThrs,'- -');hold on;
        plot(tVectorPlotNew,cfsThrsNew);
        hold off;
        xlabel('Time'); % x轴注解
        ylabel('Coefficients'); % y轴注解
        title(['Before And After Control Chart in Progress ',num2str(ControlChartNumber),'time']); % 图形标题
        legend('Before','After');
        %}
    if length(cfsThrsNew) == length(cfsThrsCurrent)
        cfsThrsCurrent = cfsThrsNew;
        flag = 1;
    else
        cfsThrsCurrent = cfsThrsNew;
    end
    
    end
    
    cfsThrsNorm = zscore(cfsThrsCurrent);
    sigma = median(abs(cfsThrsCurrent))/0.6745;
    thr   = sigma * (0.3936+0.1829*log(length(cfsThrsNorm) - 2));
    thrs(frequencyBand,idx)   = sigma * (0.3936+0.1829*log(length(cfsThrsNorm) - 2));
    sig(frequencyBand,idx)    = sigma;
    thrsel(frequencyBand,idx) = (0.3936+0.1829*log(length(cfsThrsNorm) - 2));
    %thrCal(idx) = (0.3936+0.1829*log(thrswin - 2));
        
    % State discriminiation
    cfsTempHilbert  = hilbert(cfsTempCurrent);
    cfsTempEnvelop  = sqrt(real(cfsTempHilbert).^2 + imag(cfsTempHilbert).^2);
    En        = sum(cfsTempEnvelop*(1/fs));
    RandomSig = rand(size(cfsThrsCurrent))*2*thr-thr;
    EnRandom  = sum(abs(RandomSig).^2*(1/fs))/thrswinlen;
    
    
    Synstate = En/EnRandom;
    states(frequencyBand,iStep) = 1-exp(-Synstate);
    K1(frequencyBand,iStep) = En;
    Kr(frequencyBand,iStep) = EnRandom;
    
end
        

    %{
    % Update progress bar
    cont = cont+1;
    pro = progress(pro, cont);
    %}
end


% Scale state arrary
startPoint = startPoint+winlen-step;
stateScaled = zeros(FrequencyNumber, step*nStep);
K1Scaled = zeros(FrequencyNumber, step*nStep);
KrScaled = zeros(FrequencyNumber, step*nStep);
for frequencyBand = 1:FrequencyNumber
	for iStep = 1:nStep
        stateScaled(frequencyBand,(1:step)+step*(iStep-1)) = ones(1, step)*states(frequencyBand,iStep);
        K1Scaled(frequencyBand,(1:step)+step*(iStep-1)) = ones(1, step)*K1(frequencyBand,iStep);
        KrScaled(frequencyBand,(1:step)+step*(iStep-1)) = ones(1, step)*Kr(frequencyBand,iStep);
    end
end
states = stateScaled;
K1 = K1Scaled;
Kr = KrScaled;


tVector = time(startPoint:startPoint+nStep*step-1);

%cfs = [cfsPriori, cfs];
%thrs = [zeros(1,length(cfsPriori)), thrs];



% plot Frequency Band
figure;
subplot(2,1,1);
%tVector = (1:length(DwithoutEdge(FrequencyBand,:)))./fs; 
plot(tVector,cfs(30,:));
%xlim([PlotStart PlotStart+10]);
hold on
plot(tVector,thrsel(30,:));
hold on
plot(tVector,thrs(30,:));
hold off;
%xlim([PlotStart PlotStart+10]);
xlabel('Time'); % x轴注解
ylabel('Taget Frequency Band'); % y轴注解
title([ ' CWT coeffcients frequency = ',num2str(round(frequency(30))),'hz']); % 图形标题
legend('coefficients','minimaxi','threshold');
subplot(2,1,2);
plot(tVector,states(30,:));
%xlim([PlotStart PlotStart+10]);
xlabel('Time'); % x轴注解
ylabel('Sates'); % y轴注解
title('Syncronization states'); % 图形标题


figure;
plot(tVector,K1(30,:));
hold on;
plot(tVector,Kr(30,:));
hold off;
%xlim([PlotStart PlotStart+10]);
xlabel('Time'); % x轴注解
ylabel('Values'); % y轴注解
title(['K1 & Kr frequency = ',num2str(round(frequency(30))),'hz']); % 图形标题
legend('K1','Kr');


figure;
pcolor(tVector,frequency,states);shading interp;
title('Syncgram Minimaxi '); % 图形标题
xlabel('Time'); % x轴注解
ylabel('Frequency'); % y轴注解


figure;
pcolor(tVector,frequency,K1);shading interp;
title('K1 Minimaxi '); % 图形标题
xlabel('Time'); % x轴注解
ylabel('Frequency'); % y轴注解


figure;
pcolor(tVector,frequency,Kr);shading interp;
title('K0 Minimaxi '); % 图形标题
xlabel('Time'); % x轴注解
ylabel('Frequency'); % y轴注解


figure;
pcolor(tVector,frequency,thrs);shading interp;
title('Threshold Minimaxi '); % 图形标题
xlabel('Time'); % x轴注解
ylabel('Frequency'); % y轴注解

figure;
pcolor(tVector,frequency,sig);shading interp;
title('Sigma Minimaxi '); % 图形标题
xlabel('Time'); % x轴注解
ylabel('Frequency'); % y轴注解


% Return results
switch nargout
    case 2
        varargout{1} = stateScaled;
        varargout{2} = tVector;
    case 4
        varargout{1} = stateScaled;
        varargout{2} = tVector;
        varargout{3} = cfs;
        varargout{4} = thrs;
    otherwise
        help NeuralSynCWT
        return

end

