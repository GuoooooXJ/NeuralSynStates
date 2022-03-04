function mapfile2nzmat(fileIn,fileOut,channel,type,section)
%MAPFILE2NZMAT Convert AlphaOmega .mat format to NZ .mat format.
%  Use as:
%    mapfile2nzmat(fileIn,fileOut,channel,type,section);
%  Input:
%    - fileIn, input filename
%    - fileOut, output filename
%    - channel, select channels e.g. [1,3,7].
%    - type, 'LFP' or 'RAW'.
%    - section, optional, select a time range in e.g. [10,60] (s).
%
%  Author   :  NIE Yingnan
%  Created  : Sept 9, 2020
%  Modified : Sept 9, 2020

if nargin==4
    section = 'all';
end

[data,hdr] = readaomat(fileIn,channel,type,section);
savenzmat(fileOut,data,hdr);

end

