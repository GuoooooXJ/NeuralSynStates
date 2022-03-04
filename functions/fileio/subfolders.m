function [subfolder] = subfolders(parent)
%SUBFOLDERS Get subfolder names
%   Use as:
%       [subfolder] = subfolders(parent);
%   Input:
%       parent, parent folder
%   Output:
%       subfolder, name list of subfolders
%   
%   Author: NIE Yingnan
%   Date: Mar.16th, 2020
  
dirnames = dir(parent);
subfolder = {};
for i=1:length(dirnames)
        if( isequal(dirnames(i).name, '.')||...
        isequal(dirnames(i).name, '..')||...
        ~dirnames(i).isdir)
            continue;
        end
        subfolder = [subfolder; dirnames(i).name];
end
 
end