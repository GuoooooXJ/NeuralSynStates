function obj = progress(var1, var2)
%PROGRESS Show a text progess bar in console.
%   Init: 
%       obj = progress(title, N); 
%   	- title, the title of the processing function
%       - N, the total number of steps of the processing function
%       - obj, the object used for saving the progress status
%   Update one step:
%       obj = progress(obj); 
%       - obj, the object used for saving the progress status
%   Or update with current step i:
%       obj = progress(obj, i); 
%       - i, current progress step number
%       - obj, the object used for saving the progress status
%
%   Author   : NIE Yingnan
%   Created  : June 15, 2020
%   Modified : June 24, 2020

% Init
lmax = 50;

if ischar(var1)&&isnumeric(var2)&&(var2>0)
    obj.status = 'inprogress';
    obj.N = var2;
    obj.cont = 0;
    obj.i = 0;
    
    fprintf(1,'\n   %s\n',var1);
    str = repmat(' ',1,lmax-4);
    fprintf(1,'  |-%s-|\n',str);
    fprintf(1,'%s','  ');

% Update
elseif nargin==1&&isequal(var1.status,'inprogress')
    obj = var1;
    obj.i = obj.i+1;
    obj = progress(obj, obj.i);
    
% Update
elseif isequal(var1.status,'inprogress')&&isnumeric(var2)&&(var2>=0)
    obj = var1;
    obj.i = var2;
    rate = var2/obj.N;
    
    if rate>=1
        diff = lmax - obj.cont;
        str = repmat('*',1,diff);
        fprintf(1,'%s',str);
        fprintf(1,'%s\n\n',' done');
        
        obj.status = 'done';
    else
        cont = floor(rate*lmax);
        diff = cont - obj.cont;
        obj.cont = cont;
        
        str = repmat('*',1,diff);
        fprintf(1,'%s',str);
    end

% Invalid input
else
    error("Error (progress.m): invalid input.");
end

end

