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
step  = 1;

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
title('CWT on data');

% 2hz以上
% Remove the points before and after to exclude edge effects
index = find(f>2);
WTwithoutEdge = wt(index,fs:end-fs);
frequency = f(index);

% 90hz以下
index = find(frequency<90);
WTwithoutEdge = WTwithoutEdge(index,:);
frequency = frequency(index);


centralFreq = round(frequency(end));
thrswin = round(1/centralFreq * fs * 60);

[FrequencyNumber,n] = size(WTwithoutEdge);

% Compute wavelet pocket coefficients for priori window
nStepThrs = ceil(thrswin/step);
nStepTemp = ceil((thrswin/2)/step);
thrs  = zeros(FrequencyNumber, nStepThrs);
sig   = zeros(FrequencyNumber, nStepThrs);
thrsel = zeros(FrequencyNumber, nStepThrs);
%thrCal = zeros(FrequencyNumber, nStepThrs);
%thrManuls = zeros(FrequencyNumber, nStepThrs);


% Detect synchronization state
%nStep = floor((length(WTwithoutEdgeNew(1,:))-(startPoint-1)-(winlen-step))/step);

startPointThrs = nStepThrs*step + 1;
startPointTemp = nStepTemp*step + 1;
%startPoint = nStepThrs*step + 1;
endPoint = startPointTemp + step - 1;
%endPoint = startPoint + winlen - 1;
nStep = floor((length(WTwithoutEdge(1,:))-(startPointThrs-1))/step);
%nStep = floor((length(DwithoutEdge(FrequencyBand,:))-(startPoint-1)-(winlen-step))/step);


cfs = zeros(FrequencyNumber, step*nStep);
states = zeros(FrequencyNumber, nStep);
Esig = zeros(FrequencyNumber, nStep);
Eback = zeros(FrequencyNumber, nStep);

for frequencyBand = 1:FrequencyNumber-1

centralFreq = round(f(frequencyBand));
thrswinCurrent = round(1/centralFreq * fs * 60);

Signal = WTwithoutEdge(frequencyBand,:);
SignalHibert  = hilbert(Signal);
SignalEnvelop = sqrt(real(SignalHibert).^2 + imag(SignalHibert).^2);
    
for iStep = 1:nStep
    
    idx = (1:step)+(iStep-1)*step;

    % Extract coefficients
    cfsThrs  = Signal((1:thrswinCurrent)+(iStep-1)*step);
    cfs(frequencyBand,idx) = Signal((startPointTemp:endPoint)+(iStep-1)*step);
    
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

    sigma = median(abs(cfsThrs))/0.6745;
    thr   = sigma * (0.3936+0.1829*log(length(cfsThrs) - 2));
    thrs(frequencyBand,idx)   = sigma * (0.3936+0.1829*log(length(cfsThrs) - 2));
    sig(frequencyBand,idx)    = sigma;
    thrsel(frequencyBand,idx) = (0.3936+0.1829*log(length(cfsThrs) - 2));
    %thrCal(idx) = (0.3936+0.1829*log(thrswin - 2));
        
    % State discriminiation
    En        = SignalEnvelop((startPointTemp:endPoint)+(iStep-1)*step);
    EnNorm    = 0.288*2*thr;
    
    
    Synstate = En/EnNorm;
    states(frequencyBand,idx) = 1-exp(-Synstate);
    Esig(frequencyBand,idx) = En;
    Erandom(frequencyBand,idx) = EnNorm;
    
end
        

    %{
    % Update progress bar
    cont = cont+1;
    pro = progress(pro, cont);
    %}
end

%{
% Scale state arrary
startPoint = startPointTemp;
stateScaled = zeros(FrequencyNumber, step*nStep);
EsigScaled = zeros(FrequencyNumber, step*nStep);
ErandomScaled = zeros(FrequencyNumber, step*nStep);
for frequencyBand = 1:FrequencyNumber
	for iStep = 1:nStep
        stateScaled(frequencyBand,(1:step)+step*(iStep-1)) = ones(1, step)*states(frequencyBand,iStep);
        EsigScaled(frequencyBand,(1:step)+step*(iStep-1)) = ones(1, step)*Esig(frequencyBand,iStep);
        ErandomScaled(frequencyBand,(1:step)+step*(iStep-1)) = ones(1, step)*Erandom(frequencyBand,iStep);
    end
end
states = stateScaled;
Esig = EsigScaled;
Erandom = ErandomScaled;
%}

tVector = time(startPointTemp:startPointTemp+nStep*step-1);

%cfs = [cfsPriori, cfs];
%thrs = [zeros(1,length(cfsPriori)), thrs];


%{
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
plot(tVector,Esig(30,:));
hold on;
plot(tVector,Erandom(30,:));
hold off;
%xlim([PlotStart PlotStart+10]);
xlabel('Time'); % x轴注解
ylabel('Values'); % y轴注解
title(['Etotal & Erandom frequency = ',num2str(round(frequency(30))),'hz']); % 图形标题
legend('Etotal','Erandom');
%}

figure;
pcolor(tVector,frequency,states);shading interp;
title('Syncgram Minimaxi '); % 图形标题
xlabel('Time'); % x轴注解
ylabel('Frequency'); % y轴注解


figure;
pcolor(tVector,frequency,Esig);shading interp;
title('E signal Minimaxi '); % 图形标题
xlabel('Time'); % x轴注解
ylabel('Frequency'); % y轴注解


figure;
pcolor(tVector,frequency,Eback);shading interp;
title('E random activity Minimaxi '); % 图形标题
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

