function  plotsignal(data, var)
%PLOTSIGNAL Plot signal
%  Use as:
%    plotsignal(data, time);
%  Or:
%    plotsignal(data, fs);
%  Input:
%    -data, nPoint x nChan matrix
%    -time, time vector
%    -fs, frequency of sampling
%
%  Author: Yingnan Nie
%  Date: Dec.08, 2019

nPoint = size(data,1);
nChan = size(data,2);

if isscalar(var)
    fs = var;
    time = 1/fs:1/fs:nPoint/fs;
elseif isrow(var)
    time = var(:);
else
    time = var;
end
    
fig = figure;
plot(time, data);

for i=1:nChan
    subplot(nChan,1,i);
    plot(time, data(:,i));   
end

set(fig,'windowkeypressfcn',@keypressfcn);

end

function keypressfcn(h,evt)
    axes = get(gcf,'children');
    x = get(axes(1),'XLim');
    
    key = evt.Key;
    switch key
        case 'uparrow'
            lim = [x(1), x(1)+(x(2)-x(1))*0.9];
            for i=1:length(axes)
                set(axes(i), 'XLim', lim);
            end
        case 'downarrow'
            lim = [x(1), x(1)+(x(2)-x(1))*1.1];
            for i=1:length(axes)
                set(axes(i), 'XLim', lim);
            end
        case 'leftarrow'
        	lim = x-(x(2)-x(1))/10;
            for i=1:length(axes)
                set(axes(i), 'XLim', lim);
            end
        case 'rightarrow'
            lim = x+(x(2)-x(1))/10;
            for i=1:length(axes)
                set(axes(i), 'XLim', lim);
            end
        case 'z'
            for i=1:length(axes)
            	y = get(axes(i), 'YLim');
                lim = y*0.9;
                set(axes(i), 'YLim', lim);
            end
        case 'x'
            for i=1:length(axes)
            	y = get(axes(i), 'YLim');
                lim = y*1.1;
                set(axes(i), 'YLim', lim);
            end
    end
end