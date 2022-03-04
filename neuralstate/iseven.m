function ret = iseven(n)
%ISEVEN If n is even return true, otherwise return false.
%
%   Use as:
%       ret = iseven(n);

if mod(n, 2) == 0
    ret = true;
else
    ret = false;
end