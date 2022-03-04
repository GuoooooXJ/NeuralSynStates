function arrayout = selmednums(arrayin, n)
%SELMEDNUMS Select n numbers from input array.
%
%   Use as:
%       arrayout = selmednums(arrayin, n);


nArrayin = length(arrayin);

if xor(iseven(nArrayin), iseven(n))
    startPoint = ceil((nArrayin-n)/2);
    arrayout = arrayin(startPoint:startPoint+n-1);
else
    startPoint = (nArrayin-n)/2+1;
    arrayout = arrayin(startPoint:startPoint+n-1);
end

