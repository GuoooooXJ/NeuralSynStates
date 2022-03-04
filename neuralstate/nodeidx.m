function [idxOut] = nodeidx(idxIn)
%NODEIDX Convert wavelet pocket tree node index between 'N' and '[X, Y]'
%   Matlab use a single number to index the wavelet poecket tree node,
%   but a index with two numbers like [X, Y] that indicates the level of
%   wavelet pocket tree node is more friendly for human. NODEIDX is
%   a function that can convert the index between these two formats.

    if length(idxIn)==1
        maxLevelNum = 0;
        level = 0;
        while(maxLevelNum<idxIn)
            level = level+1;
            maxLevelNum = maxLevelNum + 2^level;
        end
        num = 2^level - (maxLevelNum - idxIn) - 1;
        idxOut = [level, num];
    else
        temp = 0;
        nLevel = idxIn(1);
        for iLevel = 0:nLevel-1
            temp = temp + 2^iLevel;
        end
        idxOut = temp + idxIn(2);
    end

end

