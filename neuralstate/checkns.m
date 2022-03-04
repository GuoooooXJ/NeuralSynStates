function checkns(raw, time, state, tstate, cfs, thrs, tcfs)
%CHECKNS Plot the intermediate results for inspection
%   Use as:
%       checkns(raw, traw, state, tstate, cfs, thrs, tcfs);
%   Input:
%       raw    - raw signal
%       time   - time vector of raw signal
%       state  - synchronization states vector
%       tstate - time vector that matches the states vector
%       cfs    - wavelect pocket coeffients at target node
%       thrs   - threshods
%       tcfs   - time vector that matches the cfs vector       
%   See dnsdetection.m
%
%   Interaction:
%       'uparrow'    - zoom in
%       'downarrow'  - zoom out
%       'leftarrow'  - move left
%       'rightarrow' - move right
%
%   Author   : NIE Yingnan
%   Created  : June 15, 2020
%   Modified : June 15, 2020

section = [0,20];

fig = figure;

% Raw signal 
subplot(311);
plot(time, raw);
set(gca, 'XLim', section);

% Coefficients and thresholds 
subplot(312);
hold on;
plot(tcfs, cfs);
plot(tcfs, thrs, 'r');
plot(tcfs, -thrs, 'r');
hold off;
set(gca, 'XLim', section);

% Neural states 
subplot(313);
plot(tstate, state);
set(gca, 'XLim', section);
set(gca, 'YLim', [-0.2,1.2]);

set(fig,'windowkeypressfcn',@keypressfcn);

end

function keypressfcn(~,evt)
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
    end

end
