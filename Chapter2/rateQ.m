function R = rateQ(Q,H)
% function R = rateQ(Q,H)
%
% Calculation of the MIMO MAC rates given the channel H and the transmit covariance matrices Q
%
% Inputs
% Q: n x n x K set of covariance matrices
% H: m x n x K set of channel matrices
% Outputs
% R: achievable sum rate

K = size(Q,3);
m = size(H,1);
I_m=eye(m);
sum_2=0;
for i = 1:K
    sum = H(:,:,i)*Q(:,:,i)*H(:,:,i)';
    sum_2 = sum_2 + sum;
end
R=log2(det(I_m+sum_2));
R=real(R);
end
% TODO
