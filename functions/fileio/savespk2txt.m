function [] = savespk2txt(filename, data, channels)
%% SAVESPK2TXT Save data matrix into a spike2 .txt format file
% Use as:
%   savespk2txt(filename, data, channels);
% Input: 
%	- filename; filename that will be created
%	- data; data matrix, each cloumn is a channel
%	- channels; cell contains channel names    

    fid = fopen(filename, 'w');
    
    for i = 1:length(channels)-1
        fprintf(fid, '%s\t', channels{i});
    end
    fprintf(fid, '%s\n', channels{end});

    for i = 1:size(data, 1)
        for j = 1:size(data, 2)-1
            fprintf(fid, '%f\t', data(i, j));
        end
        fprintf(fid, '%f\n', data(i, end));
    end
    
    fclose(fid);
    
end

