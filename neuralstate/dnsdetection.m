function [varargout] = dnsdetection(signal, fs, node, varargin)
%DNSDETECTION Dynamic Neural State Detection.
%   Detect synchronization state for LFP signal.
%
%   Use as:
%       [states, tStates] = dnsdetection(signal, fs, node, options, values);
%   or:
%       [states, tStates, cfs, tCfs, thrs] = ...
%           dnsdetection(signal, fs, node, options, values);
%   or:
%       [states, tStates, cfs, tCfs, thrs, startPoint, endPoint] = ...
%           dnsdetection(signal, fs, node, options, values);
%
%   Input:
%       signal      - signal vector
%       fs          - sampling rate of input signal
%       node        - target node, should match the target frequency band
%
%   Available options:    %options('n1',1) like this
%       'n1'        - set N1 value, default N1=2
%       'n2'        - set N2 value, default N2=3
%       'basis'     - set wavelet basis, default 'sym8'
%       'step'      - set step length, default 8 points
%       'winlen'    - set window length, default 128 points
%       'thrswin'   - set priori window length, default 2 seconds
%       'cutnode'   - cut nodes from the wavelet pocket tree {[x1,y1],..}
%       'weight'    - set threshold weight, default 1
%   
%   Output:
%       states      - synchronization states vector, 0 or 1
%       tStates     - time vector that matches the states vector
%       cfs         - wavelect pocket coeffients at target node
%       tCfs        - time vector that matches the cfs vector
%       thrs        - threshods
%       startPoint  - the start point of the detection (exclude priori
%                     window)
%       endPoint    - the end point of the detection (exclude the points
%                     at the end of the signal, if they are not enough for
%                     a new window)
%
%   Reference: Huichun Luo, 2018.
%
%   Author   : NIE Yingnan
%   Created  : June 13, 2020
%   Modified : June 15, 2020

% Set default values
N1 = 2;
N2 = 3;
basis = 'sym8';
step  = 8;
winlen = 128;
thrswin = round(2*fs);
entropy = 'shannon';
cutnodes = [];
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
        case 'cutnode'
            cutnodes = value;
        case 'weight'
            weight = value;
        otherwise
            help dnsdetection
            return
    end
end


time = (1:length(signal))./fs; % Create time vector
targetLevel = node(1);

nCfsWin = ceil(winlen/2^targetLevel);
nCfsStep = ceil(step/winlen*nCfsWin);

% Compute decompose level
if ~isempty(cutnodes)
    maxLev = 0;
    for iNode = 1:length(cutnodes)
        temp = cutnodes{iNode}(1);
        if temp>maxLev
            maxLev = temp;
        end
    end
    
    if node(1)>maxLev
        maxLev = node(1);
    end
    
    wpdecLevel = maxLev;
else
    wpdecLevel = node(1);
end

% Init progress bar
N = floor((length(signal)-(winlen-step))/step);
cont = 0;
pro = progress('Neural State', N);

% Set DWT mode
ST = dwtmode('status','nodisp');

if ~isequal(ST, 'sym')
    dwtmode('sym');
end

% Compute wavelet pocket coefficients for priori window
nStep = ceil((thrswin-(winlen-step))/step);
cfsThrs = zeros(1, nCfsStep*nStep);
cfsPriori = zeros(1, nCfsStep*nStep);
tCfsPriori = zeros(1, nCfsStep*nStep);

for iStep = 1:nStep
    signalWin = signal((1:winlen)+(iStep-1)*step);
    timeWin = time((1:winlen)+(iStep-1)*step);
    timeStep = selmednums(timeWin, step);
    
    idx = (1:nCfsStep)+(iStep-1)*nCfsStep;
    tCfsPriori(idx) = resample(timeStep, nCfsStep, step);
    
    % Compute wavelet pocket tree
    wpTree = wpdec(signalWin, wpdecLevel, basis, entropy);
    
    % Cut nodes
    if ~isempty(cutnodes)
        wpTree = wptcutnodes(wpTree, cutnodes);
        wpTree = wpjoin(wpTree, node);
    end
    
    % Extract coefficients
    cfsTemp = wpcoef(wpTree, node);
    cfsTemp = selmednums(cfsTemp, nCfsWin);
    cfsTemp = cfsTemp(end-nCfsStep+1:end);
    cfsPriori(idx) = cfsTemp;
    cfsThrs(idx) = cfsTemp;
    
    % Update progress bar
    cont = cont+1;
    pro = progress(pro, cont);
end


% Detect synchronization state
startPoint = nStep*step + 1;
endPoint = startPoint + winlen - 1;
nStep = floor((length(signal)-(startPoint-1)-(winlen-step))/step);

cfs = zeros(1, nCfsStep*nStep);
states = zeros(1, nStep);
thrs = zeros(1, nCfsStep*nStep);
tCfs = zeros(1, nCfsStep*nStep);

for iStep = 1:nStep
    % Segment signal
    signalWin = signal((startPoint:endPoint)+(iStep-1)*step);
    timeWin = time((startPoint:endPoint)+(iStep-1)*step);
    timeStep = selmednums(timeWin, step);
    
    idx = (1:nCfsStep)+(iStep-1)*nCfsStep;
    tCfs(idx) = resample(timeStep, nCfsStep, step);
    
    % Compute threshold
    m = median(abs(cfsThrs))/0.6745;
    thr = thselect(cfsThrs./m, 'minimaxi')*m*weight;
    thrs(idx) = thr;
    
    % Compute wavelet pocket tree
    wpTree = wpdec(signalWin, wpdecLevel, basis, entropy);
    
    % Cut nodes
    if ~isempty(cutnodes)
        wpTree = wptcutnodes(wpTree, cutnodes);
    end
    
    % Extract coefficients
    cfsTemp = wpcoef(wpTree, node);
    cfsTemp = selmednums(cfsTemp, nCfsWin);
    cfsTemp = cfsTemp(end-nCfsStep+1:end);
    cfs(idx) = cfsTemp;
    
    % State discriminiation
    if ~isempty(find(abs(cfsTemp)>thr, 1))
        state = 1;
        idxTemp=find(cfsTemp>thr);
        cfsTemp(idxTemp)=cfsTemp(idxTemp)-thr;
        idxTemp=find(cfsTemp<-thr);
        cfsTemp(idxTemp)=cfsTemp(idxTemp)+thr;
    else
        state = 0;
    end
    
    states(iStep) = state;
    
    % Shift thresholding window
    cfsThrs(1:end-nCfsStep) = cfsThrs(nCfsStep+1:end);
    cfsThrs(end-nCfsStep+1:end) = cfsTemp;
    
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
            state = 1;
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
tStates = time(startPoint:startPoint+nStep*step-1);
%states = stateScaled;

endPoint = startPoint+length(stateScaled)-1;

cfs = [cfsPriori, cfs];
tCfs = [tCfsPriori, tCfs];
thrs = [zeros(1,length(cfsPriori)), thrs];

% Return results
switch nargout
    case 2
        varargout{1} = stateScaled;
        varargout{2} = tStates;
    case 5
        varargout{1} = stateScaled;
        varargout{2} = tStates;
        varargout{3} = cfs;
        varargout{4} = tCfs;
        varargout{5} = thrs;
    case 7
        varargout{1} = stateScaled;
        varargout{2} = tStates;
        varargout{3} = cfs;
        varargout{4} = tCfs;
        varargout{5} = thrs;
        varargout{6} = startPoint;
        varargout{7} = endPoint;
    otherwise
        help dnsdetection
        return
end
