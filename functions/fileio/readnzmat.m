function [data, hdr] = readnzmat(filename)
%READNZMAT Save data to file as NZ .mat format.
%  Use as:
%    [data, hdr] = readnzmat(filename)
%  Input:
%    - filename, input file name.
%  Output:
%    - data, data matrix, each column as a channel.
%    - hdr, structure, saving header information
%         .fs - sampling rate
%         .nchan - number of channels
%         .chans - cell, channel names, optional
%         .event - cell, {idx1, 'msg1'; idx2, 'msg2'}, optional
%         .note - string for notes, optional
%
%  Author   :  NIE Yingnan
%  Created  : Sept 9, 2020
%  Modified : Sept 9, 2020

% TODO: check input format

data = load(filename,'data');
hdr = load(filename,'hdr');

data = data.data;
hdr = hdr.hdr;

end

