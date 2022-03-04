function Rx = Rxx(x,P)
% x      - 自变量
% P      - 自相关数目
% Rxx[m] = sum(x[n]*x[n+m])
%         (n = 0->inf,m = 0->P)

N  = length(x);

sizex = size(x);
if sizex(1)==1
    % x 为行向量那么Rx也为行向量
    Rx = zeros(1,P+1); 
    
    for m = 0:P
        Rx(m+1) = x(m+1:N)*x(1:N-m)';
    end
    
elseif sizex(2)==1
    % x 为列向量那么Rx也为列向量
    Rx = zeros(P+1,1); 
    for m = 0:P
        Rx(m+1) = x(m+1:N)'*x(1:N-m);
    end
    
end
