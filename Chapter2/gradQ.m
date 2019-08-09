function g = gradQ(Q,H)
% function g = gradQ(Q,H)
%
% Calculation of the gradient of the MIMO MAC rates with respect to the transmit covariance matrices Q given the channel H
%
% Inputs
% Q: n x n x K set of covariance matrices
% H: m x n x K set of channel matrices
% Outputs
% g: n x n x K set of derivatives V_1, ..., V_K 

K = size(Q,3);
m = size(H,1);
n = size(H,2);
g = zeros(n,n,K);
I_m=eye(m);
sum_2=0;
for l = 1:K
    sum = H(:,:,l)*Q(:,:,l)*H(:,:,l)';
    sum_2 = sum_2 + sum;
end
for i = 1: K
    V=(1/log(2))*H(:,:,i)'*(I_m+sum_2)^(-1)*H(:,:,i);
    g(:,:,i)=V;
end
end
% TODO
