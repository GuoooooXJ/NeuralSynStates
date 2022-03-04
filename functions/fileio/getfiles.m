function list = getfiles(folder,extension)
%GETFILES Get a list of file names of the input folder.
%   Use as:
%       list = getfiles(folder,extension);
%   Input:
%       folder, the path of the folder, i.e. '.txt'
%       extension, specify the extension of the files
%   Output:
%       list, the list of the file names
%
%   Author   : NIE Yingnan
%   Created  : June 23th, 2020
%   Modified : June 23th, 2020
if(nargin==1)
    extension = [];
end

type = '*';
type = strcat(type,extension);

dirnames = dir(fullfile(folder, type));
list = [];
for i=1:length(dirnames)
        if ~dirnames(i).isdir
            list{i,1} = dirnames(i).name;
        end     
end

end

