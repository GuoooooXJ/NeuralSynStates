function [count,timeMax,timeMin,timeMean,timeStd,PeakMax,PeakMin,PeakMean,PeakStd] = zero_crossings(x,thr)
% x是时域信号，必须是1位的行向量或者列向量;
% thr是设定的阈值
% count是该函数计算信号x过零率的值
 
% initial value
count = 0;
idx = [];
Peak = [];
 
% error checks
if(length(x) == 1)
    error('ERROR: input signal must have more than one element');
end
 
if((size(x, 2) ~= 1) && (size(x, 1) ~= 1))
    error('ERROR: Input must be one-dimensional');
end
    
% force signal to be a vector oriented in the same direction
x = x(:);
 
num_samples = length(x);
for i=2:num_samples
    % Any time you multiply to adjacent values that have a sign difference
    % the result will always be negative.  When the signs are identical,
    % the product will always be positive.
    if(((x(i)*x(i-1)) < 0)&&(abs(x(i)-x(i-1))>thr))
        count = count + 1;
        idx(count) = i;
    end
end

idxNext = circshift(idx,-1);

for k=1:length(idx)-1
    segment = x(idx(k):idxNext(k)); 
    a = max(segment);
    Peak(k) = max(segment);
end


TimeInter = diff(idx);
timeMax = max(TimeInter);
timeMin = min(TimeInter);
timeMean = mean(TimeInter);
timeStd  = std(TimeInter);

PeakInter = diff(Peak);
PeakMax = max(PeakInter);
PeakMin = min(PeakInter);
PeakMean = mean(PeakInter);
PeakStd  = std(PeakInter);


    
