function [roi] = selectroi(var, shape)
%% SELECTROI Select ROI (Region Of Interest) from a martix.
% Use as:
%   [roi] = selectroi(var, shape);
% Input:
%   - var, the input variable, must be a matrix
%   - shape, the shape of the ROI, can be 'rectangle','circle','ellipse',or
%   'polygon'
% Output:
%   - roi, selected ROI (a matrix with the same size of var)
%
% Author: NIE Yingnan
% Mar 31th, 2020

switch shape
    case 'rectangle'
    case 'circle'
    case 'ellipse'
    case 'polygon'
    otherwise
        help selectroi;
        return
end

[ylim,xlim,~] = size(var);

% Plot the figure
figure;
imagesc(var);
colormap(jet);
colorbar;
axis('square');
set(gca,'Box','on');
title('Select the ROI then press Enter');

% Begin interactive placement of the ROI
roi = zeros(size(var));
switch shape
    case 'rectangle'
        h = drawrectangle();
        pause;
        pos = h.Position;
        for y=1:ylim
            for x=1:xlim
                roi(y,x) = x>pos(1)&&x<(pos(1)+pos(3))&&y>pos(2)&&y<(pos(2)+pos(4));
            end
        end
    case 'circle'
        h = drawcircle();
        pause;
        x0 = h.Center(1);
        y0 = h.Center(2);
        r = h.Radius;
        for y=1:ylim
            for x=1:xlim
                d2 = (x-x0)^2+(y-y0)^2;
                roi(y,x) = d2<r^2;
            end
        end
    case 'ellipse'
        h = drawellipse();
        pause;
        x0 = h.Center(1);
        y0 = h.Center(2);
        a = h.SemiAxes(1);
        b = h.SemiAxes(2);
        for y=1:ylim
            for x=1:xlim
                d = (x-x0)^2/a^2+(y-y0)^2/b^2;
                roi(y,x) = d<1;
            end
        end
    case 'polygon'
        h = drawpolygon();
        pause;
        XV = h.Position(:,1);
        YV = h.Position(:,2);
        X = repmat(1:xlim,ylim,1);
        Y = repmat((1:ylim)',1,xlim);
        roi = inpolygon(X,Y,XV,YV);
end

end
