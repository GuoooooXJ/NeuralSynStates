function [] = savenspltxt(filename, data, fs)
%% SAVENSPLTXT Save data matrix into a nspl .txt format file
% Use as:
%   savenspltxt(filename, data, fs);
% Input: 
% 	- filename: filename that will be created
%   - data: data matrix, each cloumn is a channel
%   - fs: sampling rate   

    fid = fopen(filename, 'w');
    fprintf(fid, '%d\t%d\t%d\t%d\t\n', fs, size(data,2), size(data,1), 0);

    for i = 1:size(data, 1)
        for j = 1:size(data, 2)-1
            fprintf(fid, '%f\t', data(i, j));
        end
        fprintf(fid, '%f\n', data(i, end));
    end
    
    fclose(fid);
    
end

