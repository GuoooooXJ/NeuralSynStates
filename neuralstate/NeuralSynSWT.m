function [varargout] = NeuralSynSWT(signal, fs, node, centralFreq, frequencyBandPlot, varargin)
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

basis = 'sym8';
step  = round(1/centralFreq * fs);
winlen = round(1/centralFreq * fs * 4);
windowTemp = winlen/fs;
windowThrs = 1/centralFreq * 60;
thrswin = round(1/centralFreq * fs * 60);

% Get customized settings
nArg = nargin - 5;
if mod(nArg, 2) ~= 0
    help dnsdetection
    return
end

    
% Cut signal to fit SWT
layerNum = node(1);
FrequencyBand = node(2);
n  = fix(length(signal)/(2^layerNum));
signal = signal(1:n*(2^layerNum));

% Create time vector
time = (1:length(signal))./fs; 

% Init progress bar
N = floor((length(signal)-(winlen-step))/step);
cont = 0;
pro = progress('Neural State', N);


% Compute SWT coefficients sequence

[swa,swd] = swt(signal,layerNum,basis);

% Remove the points before and after to exclude edge effects
DwithoutEdge = swd(:,fs:end-fs);
AwithoutEdge = swa(:,fs:end-fs);

% Compute wavelet pocket coefficients for priori window
nStepThrs = ceil((thrswin-(winlen-step))/step);
thrs  = zeros(1, nStepThrs);
thrManuls = zeros(1, nStepThrs);
sig   = zeros(1, nStepThrs);
thrsel = zeros(1, nStepThrs);
thrCal = zeros(1, nStepThrs);
%cfsPriori = zeros(1, step*nStep);

%{
for iStep = 1:nStep
    
    idx = (1:step)+(iStep-1)*step;
    
    % Extract coefficients
    cfsThrs = DwithoutEdge(FrequencyBand,(1:thrswin)+(iStep-1)*step);
    %cfsTemp = cfsTemp(end-step+1:end);
    
    % Control chart
    flag = 0;
    ControlChartNumber = 0;
    while(flag == 0)
    ControlChartNumber = ControlChartNumber+1;
    meancfsThrs  = mean(cfsThrs);
    sigmacfsThrs = std(cfsThrs);
    cfsThrsNew = cfsThrs(cfsThrs<meancfsThrs+2.58*sigmacfsThrs & cfsThrs>meancfsThrs-2.58*sigmacfsThrs);
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
    if length(cfsThrsNew) == length(cfsThrs)
        cfsThrs = cfsThrsNew;
        flag = 1;
    else
        cfsThrs = cfsThrsNew;
    end
    
    end

    
    cfsThrsNorm = zscore(cfsThrs);
    sigma = median(abs(cfsThrsNorm))/0.6745;
    thrManul   = sigma * (0.3936+0.1829*log(thrswin - 2));
    thr        = sigma * thselect(cfsThrsNorm,'rigrsure') ;
    thrs(idx)      = thr;
    thrManuls(idx) = thrManul;
    sig(idx)    = sigma;
    thrsel(idx) = thselect(cfsThrsNorm,'rigrsure');
    thrCal(idx) = (0.3936+0.1829*log(thrswin - 2));
    
    RandomSig = rand(size(cfsThrs))*2*thr-thr;
    EnRandom  = sum(abs(RandomSig)*(1/fs));
    Kr(iStep) = EnRandom/windowThrs;
    
    EnNorm      = 2*thr*winlen*(1/fs);
    
    
    
    
    % Update progress bar
    cont = cont+1;
    pro = progress(pro, cont);
end
%}


% Detect synchronization state
startPoint = nStepThrs*step + 1;
startPointMean = nStepThrs*step + 1;
endPoint = startPoint + winlen - 1;
nStep = floor((length(DwithoutEdge(FrequencyBand,:))-(startPoint-1)-(winlen-step))/step);

cfs = zeros(1, step*nStep);
states = zeros(1, nStep);
K1 = zeros(1, nStep);
K0 = zeros(1, nStep);
Kr = zeros(1, nStep);
Ks = zeros(1, nStep);


for iStep = 1:nStep
    
    idx = (1:step)+(iStep-1)*step;

    % Extract coefficients
    cfsTemp = DwithoutEdge(FrequencyBand,(startPoint:endPoint)+(iStep-1)*step);
    cfsThrs = DwithoutEdge(FrequencyBand,(1:thrswin)+(iStep-1)*step);
    cfs(idx) = cfsTemp(end-step+1:end);
    
    % Compute threshold
    % According to the paper, Normlization before median


    % Control chart
    flag = 0;
    ControlChartNumber = 0;
    while(flag == 0)
    ControlChartNumber = ControlChartNumber+1;
    meancfsThrs  = mean(cfsThrs);
    sigmacfsThrs = std(cfsThrs);
    cfsThrsNew = cfsThrs(cfsThrs<meancfsThrs+2.58*sigmacfsThrs & cfsThrs>meancfsThrs-2.58*sigmacfsThrs);
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
    if length(cfsThrsNew) == length(cfsThrs)
        cfsThrs = cfsThrsNew;
        flag = 1;
    else
        cfsThrs = cfsThrsNew;
    end
    
    end
        
    cfsThrsNorm = zscore(cfsThrs);
    sigma = median(abs(cfsThrs))/0.6745;
    %thrManul   = sigma * (0.3936+0.1829*log(thrswin - 2));
    thr        = mean(thrs(1:step*MeanStep)+(iStep-1)*step);
    thrs(idx)      = sigma * thselect(cfsThrsNorm,'minimaxi');
    %thrManuls(idx) = thrManul;
    sig(idx)    = sigma;
    thrsel(idx) = thselect(cfsThrsNorm,'minimaxi');
    %thrCal(idx) = (0.3936+0.1829*log(thrswin - 2));
        
    % State discriminiation
    %{
    cfsTempAbs  = abs(cfsTemp);
    En          = sum(cfsTempAbs*(1/fs));
    EnNorm      = 2*thr*winlen*(1/fs);
    RandomSig = rand(size(cfsThrs))*2*thr-thr;
    EnRandom  = sum(abs(RandomSig)*(1/fs));
    
    
    Synstate = En/(0.288*EnNorm);
    states(iStep) = 1-exp(-Synstate);
    K1(iStep) = En/windowTemp;
    K0(iStep) = 0.288*EnNorm/windowTemp;
    Kr(iStep) = EnRandom/windowThrs;
    %}
    
    % State discriminiation
    cfsTempHilbert  = hilbert(cfsTempCurrent);
    cfsTempEnvelop  = sqrt(real(cfsTempHilbert).^2 + imag(cfsTempHilbert).^2);
    En        = sum(cfsTempEnvelop*(1/fs));
    EnRandom  = 1/3*(thr).^2;
    
    
    Synstate = En/EnRandom;
    states(iStep) = 1-exp(-Synstate);
    K1(iStep) = En;
    Kr(iStep) = EnRandom;
   
    % Update progress bar
    cont = cont+1;
    pro = progress(pro, cont);
end


% Scale state arrary
startPoint = startPointMean+winlen-step;
stateScaled = zeros(1, step*nStep);
K1Scaled = zeros(1, step*nStep);
K0Scaled = zeros(1, step*nStep);
for iStep = 1:nStep
    stateScaled((1:step)+step*(iStep-1)) = ones(1, step)*states(iStep);
    K1Scaled((1:step)+step*(iStep-1)) = ones(1, step)*K1(iStep);
    K0Scaled((1:step)+step*(iStep-1)) = ones(1, step)*K0(iStep);
    KrScaled((1:step)+step*(iStep-1)) = ones(1, step)*Kr(iStep);
end
states = stateScaled;
K1 = K1Scaled;
K0 = K0Scaled;
Kr = KrScaled;


tVector = time(startPointMean:startPointMean+nStep*step-1);
endPoint = startPointMean+length(states)-1;

%cfs = [cfsPriori, cfs];
%thrs = [zeros(1,length(cfsPriori)), thrs];


%{
% plot Frequency Band
figure;
subplot(3,1,1);
%tVector = (1:length(DwithoutEdge(FrequencyBand,:)))./fs; 
plot(tVector,cfs);
xlabel('Time'); % x轴注解
ylabel('Taget Frequency Band'); % y轴注解
title([ frequencyBandPlot,' window = ',num2str(windowThrs),'s']); % 图形标题
%xlim([PlotStart PlotStart+10]);
hold on
plot(tVector,thrManuls);
hold off;
%xlim([PlotStart PlotStart+10]);
legend('coefficients','minimaxi');
%legend('coefficients','threshold','thselect','thrManul');

subplot(3,1,2);
plot(tVector,cfs);
xlabel('Time'); % x轴注解
ylabel('Taget Frequency Band'); % y轴注解
title([ frequencyBandPlot,' window = ',num2str(windowThrs),'s']); % 图形标题
%xlim([PlotStart PlotStart+10]);
hold on
plot(tVector,thrs);
hold off;
%xlim([PlotStart PlotStart+10]);
legend('coefficients','rigrsure');


subplot(3,1,3);
plot(tVector,states);
%xlim([PlotStart PlotStart+10]);
xlabel('Time'); % x轴注解
ylabel('Sates'); % y轴注解
title('Syncronization states'); % 图形标题
%}

figure;
plot(tVector,K1);hold on;
plot(tVector,K0);hold on;
plot(tVector,Kr);
hold off;
xlabel('Time'); % x轴注解
ylabel('Values'); % y轴注解
title(['Rigrsure K1 & K0 & Kr in ', frequencyBandPlot,' window = ',num2str(windowThrs),'s']); % 图形标题
legend('K1','K0','Kr');


figure;
plot(tVector,sig);hold on;
%plot(tVector,thsel);hold on;
plot(tVector,thrMean);
plot(tVector,thrsel);
hold off;
xlabel('Time'); % x轴注解
ylabel('Values'); % y轴注解
title(['Sigma & Thselect in ', frequencyBandPlot,' window = ',num2str(windowThrs),'s']); % 图形标题
legend('Sigma','minimaxi','rigrsure');
%legend('Sigma','Thselect','ThrMaunl');


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
    case 6
        varargout{1} = stateScaled;
        varargout{2} = tVector;
        varargout{3} = cfs;
        varargout{4} = thrs;
        varargout{5} = startPoint;
        varargout{6} = endPoint;
    otherwise
        help NeuralSynSWT
        return

end

