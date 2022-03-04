function [A,Sgm2] = AR(x,P)

%  x   - 输入数据，列向量
%  P   - AR 阶数,标量
%  A   - AR 系数[a1 a2 ...aP]，行向量;
% Sgm2 - 方差
A  = zeros(1,P);  % 初始化参数矩阵，行向量
Rx = Rxx(x,P);    % 自相关估计，size = P+1

% p = 1 时候的Yule-Walker方程解,注意，MATLAB下标从1开始
A(1) = - Rx(1+1)/Rx(0+1);
Sgm2  = Rx(0+1)*(1 - A(1)^2);

% p = 1->P 时候的递推,注意，MATLAB下标从1开始
for p = 1:P-1
    k = 1:p;
    K = -(Rx((p+1)+1) + A(k)*Rx((p+1-k)+1))/Sgm2;
    Sgm2 = Sgm2 * (1 - K*K);
    
    A(k) = A(k) + K * A(p+1-k);
    A(p+1) = K;
end
