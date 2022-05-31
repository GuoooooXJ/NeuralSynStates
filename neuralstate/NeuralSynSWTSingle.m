function [varargout] = NeuralSynSWTSingle(signal, fs, node, centralFreq ,startFreq, endFreq,varargin)
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
step  = 1;
windowThrs = 1/centralFreq * 60;
thrswin = round(1/centralFreq * fs * 60);

% Get customized settings
nArg = nargin - 6;
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
N = floor(length(signal)/step);
cont = 0;
pro = progress('Neural State', N);


% Compute SWT coefficients sequence
[swa,swd] = swt(signal,layerNum,basis);

%{
mzero=zeros(size(swd));
A=mzero;
    %Reduction of the Approximate coefficients of the last layer
A(layerNum,:)=iswt(swa,mzero,basis);

    %Reduction of the Detail coefficients of every layer
D=mzero;
for i=1:layerNum
   swcfs=mzero;
   swcfs(i,:)=swd(i,:);
   D(i,:)=iswt(mzero,swcfs,basis);
end
    %Reduction of the Approximate coefficients of every layer
for i=1:layerNum-1
    A(layerNum-i,:) = A(layerNum+1-i,:) + D(layerNum+1-i,:);
end
%}

%{
figure;
kp = 0;
for i = 1:layerNum
subplot(layerNum,2,kp+1), plot(time,swa(i,:));
title(['Approx. cfs level ',num2str(i)])
xlim([1 5]);
subplot(layerNum,2,kp+2), plot(time,swd(i,:));
title(['Detail cfs level ',num2str(i)])
xlim([1 5]);
kp = kp + 2;
end
%}

% Remove the points before and after to exclude edge effects
DwithoutEdge = swd(:,fs:end-fs);
AwithoutEdge = swa(:,fs:end-fs);

Signal = DwithoutEdge(FrequencyBand,:);
SignalHibert  = hilbert(Signal);
SignalEnvelop = sqrt(real(SignalHibert).^2 + imag(SignalHibert).^2);

% Compute wavelet pocket coefficients for priori window
nStepThrs = ceil(thrswin/step);
nStepTemp = ceil((thrswin/2-step/2)/step);
thrs  = zeros(1, nStepThrs);
sig   = zeros(1, nStepThrs);
thrsel = zeros(1, nStepThrs);
%cfsPriori = zeros(1, step*nStep);


% Detect synchronization state
startPointThrs = nStepThrs*step + 1;
startPointTemp = nStepTemp*step + 1;
%startPoint = nStepThrs*step + 1;
endPoint = startPointTemp + step - 1;
%endPoint = startPoint + winlen - 1;
%nStep = floor((length(DwithoutEdge(FrequencyBand,:))-(startPointThrs-1))/thrswin);
nStep = floor((length(Signal)-(startPointThrs-1))/step);

cfs = zeros(1, step*nStep);
states = zeros(1, nStep);
K1 = zeros(1, nStep);
Kr = zeros(1, nStep);

%{
states = zeros(1, step*nStep);
K1 = zeros(1, step*nStep);
Kr = zeros(1, step*nStep);
%}

for iStep = 1:nStep
    
    idx = (1:step)+(iStep-1)*step;

    % Extract coefficients
    cfsThrs  = Signal((1:thrswin)+(iStep-1)*step);
    cfs(idx) = Signal((startPointTemp:endPoint)+(iStep-1)*step);
    
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
    plot(tVectorPlot,cfsThrs,'-');hold on;
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
    %thrManul   = sigma * (0.3936+0.1829*log(thrswin - 2));
    thr        = sigma * (0.3936+0.1829*log(length(cfsThrs) - 2));
    thrs(idx)  = thr;
    %thrManuls(idx) = thrManul;
    sig(idx)    = sigma;
    %thrsel(idx) = thselect(cfsThrsNorm,'minimaxi');
    thrsel(idx) = (0.3936+0.1829*log(length(cfsThrs) - 2));
        
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
    En        = SignalEnvelop((startPointTemp:endPoint)+(iStep-1)*step);
    EnNorm    = 0.288*2*thr;
    
    Synstate = En/EnNorm;
    
    states(idx) = 1-exp(-Synstate);
    K1(idx) = En;
    Kr(idx) = EnNorm;
   
    % Update progress bar
    cont = cont+1;
    pro = progress(pro, cont);
end

startPoint = startPointTemp;
%startPoint = startPoint+winlen-step;


% Scale state arrary

stateScaled = zeros(1, step*nStep);
K1Scaled = zeros(1, step*nStep);
KrScaled = zeros(1, step*nStep);
for iStep = 1:nStep
    stateScaled((1:step)+step*(iStep-1)) = ones(1, step)*states(iStep);
    K1Scaled((1:step)+step*(iStep-1)) = ones(1, step)*K1(iStep);
    KrScaled((1:step)+step*(iStep-1)) = ones(1, step)*Kr(iStep);
end
states = stateScaled;
K1 = K1Scaled;
Kr = KrScaled;


tVector = time(startPoint:startPoint+nStep*step-1);
endPoint = startPoint+length(states)-1;
%tVector = time(startPointTemp:startPointTemp+nStep*step-1);
%endPoint = startPointTemp+length(states)-1;

%cfs = [cfsPriori, cfs];
%thrs = [zeros(1,length(cfsPriori)), thrs];

%{
statesRe = resample(states, 1/0.02, fs);
statesRe = resample(statesRe, fs, 1/0.02);
statesRe = statesRe(1:length(states));
%}



% plot Frequency Band
figure;
PlotStart = windowThrs/2;
%PlotStart = startPointTemp/fs;
subplot(2,1,1);
plot(tVector,cfs);
xlabel('Time'); % x轴注解
ylabel('Taget Frequency Band'); % y轴注解
title([ num2str(startFreq),'hz -' ,num2str(endFreq),'hz window = ',num2str(windowThrs),'s' ]); % 图形标题
xlim([PlotStart PlotStart+10]);
hold on
plot(tVector,thrs);
hold off;
xlim([PlotStart PlotStart+10]);
legend('频段系数','minimax阈值');
subplot(2,1,2);
%plot(tVector,statesRe);
plot(tVector,states);
xlim([PlotStart PlotStart+10]);
%xlim([PlotStart PlotStart+10]);
xlabel('Time'); % x轴注解
ylabel('Sates'); % y轴注解
title('同步化状态'); % 图形标题


figure;
plot(tVector,K1);
xlim([PlotStart PlotStart+10]);
hold on;
plot(tVector,Kr);xlim([PlotStart PlotStart+10]);
hold off;
xlabel('Time'); % x轴注解
ylabel('Values'); % y轴注解
title(['信号能量和随机活动能量 ', num2str(startFreq) , 'hz -' , num2str(endFreq),'hz window = ',num2str(windowThrs),'s']); % 图形标题
legend('信号能量','随机活动能量');


figure;
plot(tVector,sig);xlim([PlotStart PlotStart+10]);hold on;
plot(tVector,thrsel);xlim([PlotStart PlotStart+10]);hold on;
plot(tVector,thrs);xlim([PlotStart PlotStart+10]);
hold off;
xlabel('Time'); % x轴注解
ylabel('Values'); % y轴注解
title(['噪声方差 & 阈值 in ', num2str(startFreq) , 'hz -' , num2str(endFreq),'hz window = ',num2str(windowThrs),'s']); % 图形标题
legend('噪声方差','minimaxi阈值','总阈值');
%legend('Sigma','Thselect','ThrMaunl');


%{
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
%}

% Return results
switch nargout
    case 2
        varargout{1} = states;
        varargout{2} = tVector;
    case 4
        varargout{1} = states;
        varargout{2} = tVector;
        varargout{3} = cfs;
        varargout{4} = thrs;
    case 6
        varargout{1} = states;
        varargout{2} = tVector;
        varargout{3} = cfs;
        varargout{4} = thrs;
        varargout{5} = startPointTemp;
        varargout{6} = endPoint;
    otherwise
        help NeuralSynSWT
        return
        

end

