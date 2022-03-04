function [varargout] = dnsdetectionSWTnew(signal, fs, node, varargin)

%DNSDETECTION Dynamic Neural State Detection .
% States are divided into 3 types Strong syn, Weak syn and de-syn
%   Detect synchronization state for LFP signal With SWT.
%
%   Use as:
%       [states, tVector] = dnsdetection(signal, fs, node, options, values);
%   or:
%       [states, tVector, cfs, thrs] = ...
%           dnsdetection(signal, fs, node, options, values);
%   or:
%       [states, tVector, cfs, thrs, startPoint, endPoint] = ...
%           dnsdetection(signal, fs, node, options, values);
%
%   Input:
%       signal      - signal vector
%       fs          - sampling rate of input signal
%       node        - target node, should match the target frequency band
%
%   Available options:    %options('n1',1) like this
%       'n1'           - set N1 value, default N1=2
%       'n2'           - set N2 value, default N2=3
%       'basis'        - set wavelet basis, default 'sym8'
%       'step'         - set step length, default 8 points
%       'winlen'       - set window length, default 128 points
%       'thrswinlen'   - set priori window length, default 2 seconds
%       'weight'       - set threshold weight, default 1
%   
%   Output:
%       states      - synchronization states vector, 0 or 1
%       tVector     - time vector that matches the states vector,coeffients
%       vector
%       cfs         - wavelect pocket coeffients at target node
%       thrs        - threshods
%       startPoint  - the start point of the detection (exclude priori
%                     window)
%       endPoint    - the end point of the detection (exclude the points
%                     at the end of the signal, if they are not enough for
%                     a new window)
%
%   Reference: Yingnan Nie
%
%   Author   : Xuanjun Guo
%   Created  : Nov 24, 2021
%   Modified : Nov 24, 2021

% Set default values

% Set default values
N1 = 2;
N2 = 3;
PlotStart = 3;
basis = 'sym8';
step  = 32;
winlen = 128;
thrswin = round(PlotStart*fs);
weight = 1;



% Get customized settings 
nArg = nargin - 3;
if mod(nArg, 2) ~= 0
    help dnsdetection
    return
end

for iVar = 1:2:nArg
    var = varargin{iVar};
    value = varargin{iVar+1};
    
    switch var
        case 'n1'
            N1 = value;
        case 'n2'
            N2 = value;
        case 'basis'
            basis = value;
        case 'step'
            step = value;
        case 'winlen'
            winlen = value;
        case 'thrswin'
            thrswin = round(value*fs);
        case 'weight'
            weight = value;
        otherwise
            help dnsdetection
            return
    end
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
% Remove the points before and after to exclude edge effects
DwithoutEdge = swd(:,fs:end-fs);
AwithoutEdge = swa(:,fs:end-fs);

% Compute wavelet pocket coefficients for priori window and syncron threshold

%nStep = ceil((thrswin-(winlen-step))/step);
startPoint = thrswin + 1;
endPoint = startPoint + winlen - 1;
nStep = floor((length(DwithoutEdge(FrequencyBand,:))-(startPoint-1)-(winlen-step))/step);

cfsThrs = zeros(1, step*nStep);
cfsThrsOnSyn = zeros(1, step*nStep);

for iStep = 1:nStep
    
    idx = (1:step)+(iStep-1)*step;
    
    % Extract coefficients
    cfsTemp = DwithoutEdge(FrequencyBand,(1:thrswin)+(iStep-1)*step);
    
    % Compute threshold
    %Amplitude
    m = median(abs(cfsTemp))/0.6745;
    thr = thselect(cfsTemp./m, 'minimaxi')*m*weight;
    cfsThrs(idx) = thr;
    %Syncro
    cfsTemp2 = cfsTemp.*cfsTemp;
    cfsTempEn = sum(cfsTemp2);
    cfsThrsOnSyn(idx) = cfsTempEn * (winlen/thrswin) * 2 ;%超过了本来应该能量占比的2倍
    
    % Update progress bar
    %cont = cont+1;
    %pro = progress(pro, cont);
end



% Detect synchronization state
states = zeros(1, nStep);
cfs = zeros(1, step*nStep);

for iStep = 1:nStep
    
    idx = (1:step)+(iStep-1)*step;
    currentIdx = ceil(median(idx));

    % Extract coefficients
    cfsTemp = DwithoutEdge(FrequencyBand,(startPoint:endPoint)+(iStep-1)*step);
    cfsTemp2 = cfsTemp.*cfsTemp;
    cfsTempEn = sum(cfsTemp2);
    cfs(idx) = cfsTemp(end-step+1:end);
    
    % State discriminiation
    if ~isempty(find(abs(cfsTemp)>cfsThrs(currentIdx), 1)) && cfsTempEn > cfsThrsOnSyn(currentIdx)
        state = 2;
        %{
        idxTemp=find(cfsTemp>thr);
        cfsTemp(idxTemp)=cfsTemp(idxTemp)-thr;
        idxTemp=find(cfsTemp<-thr);
        cfsTemp(idxTemp)=cfsTemp(idxTemp)+thr;
        %}
        
    elseif ~isempty(find(abs(cfsTemp)>cfsThrs(currentIdx), 1)) && cfsTempEn <= cfsThrsOnSyn(currentIdx)
        state = 1;
        %{
        idxTemp=find(cfsTemp>thr);
        cfsTemp(idxTemp)=cfsTemp(idxTemp)-thr;
        idxTemp=find(cfsTemp<-thr);
        cfsTemp(idxTemp)=cfsTemp(idxTemp)+thr;
        %}
        
    else
        state = 0;
    end
    
    states(iStep) = state;
    
    % Update progress bar
    cont = cont+1;
    pro = progress(pro, cont);
end



% Correct states with N1/N2
statesCorrected = zeros(size(states));
prioriState = 0;

for iStep = 1:nStep
    % Boundary
    if (iStep+N1)>nStep
        statesN1 = states(iStep+1:end);
    else
        statesN1 = states(iStep+1:iStep+N1);
    end
    
    if (iStep+N2)>nStep
        statesN2 = states(iStep+1:end);
    else
        statesN2 = states(iStep+1:iStep+N2);
    end
    
    % Correct current state with posterior states
    if states(iStep)~=0
        if all(statesN1)
            state = states(iStep);
        else
            state = prioriState;
        end
    else
        if any(statesN2)
            state = prioriState;
        else
            state = 0;
        end
    end
    
    statesCorrected(iStep) = state;
    prioriState = state;
end


% Scale state arrary
startPoint = startPoint+winlen-step;
stateScaled = zeros(1, step*nStep);
for iStep = 1:nStep
    stateScaled((1:step)+step*(iStep-1)) = ones(1, step)*statesCorrected(iStep);
end
tVector = time(startPoint:startPoint+nStep*step-1);

%states = stateScaled;

endPoint = startPoint+length(stateScaled)-1;

%cfs = [cfsPriori, cfs];
%thrs = [zeros(1,length(cfsPriori)), thrs];

% plot Frequency Band
figure;
subplot(2,1,1);
%tVector = (1:length(DwithoutEdge(FrequencyBand,:)))./fs; 
plot(tVector,cfs);
xlabel('Time'); % x轴注解
ylabel('Taget Frequency Band'); % y轴注解
title('Stationary Wavelet Transform in 4-8hz Frequency band'); % 图形标题
xlim([PlotStart 10000/fs]);
hold on
plot(tVector,cfsThrs);
xlim([PlotStart 10000/fs]);
legend('coefficients','Amthreshold');

subplot(2,1,2);
plot(tVector,stateScaled);
ylim([-3 3]);
xlim([PlotStart 10000/fs]);
xlabel('Time'); % x轴注解
ylabel('Sates'); % y轴注解
title('States based on SWT in 4-8hz Frequency band'); % 图形标题


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
        help dnsdetection
        return


end