function [a, sgm] = ARLSS(x,p)

N = length(x);
A = zeros(N-1, p);
for m = 1:N-1
    for n = 1:p
        if n > m
            A(m, n) = 0;
        else
            A(m, n) = x(m-n+1);
        end
    end
end
a = inv(A'* A)*(A')*(- x(2:end)');

sgm = norm(-x(2:end)'-A*a,2);
end