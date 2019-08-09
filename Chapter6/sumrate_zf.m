function [R q] = sumrate_zf(H,P,a)
% function [R q] = sumrate_zf(H,P,a)
%
% Sumrate with zero-forcing filter
%
% Input
%  H: m x 1 x K set of channel vectors
%  P: maximum sum power
%  a: K x 1 column vector which indicates the active users, i.e., a(k)==1 if
%     user k is active and a(k)==0 otherwise
%
% Output
%  R: sum rate of the waterfilling solution
%  q: power allocation of the waterfilling solution
%

% TO DO
[m,~,K] = size(H);
H = reshape(H,m,K);
for i=1:K
    if a(i)==0
        H(:,i)=0;
    end
end
H(:,all(H==0,1)) = [];
K=size(H,2);
Vk = inv(H'*H)*H';
%vk = Vk'*a;
sigma_square=zeros(K,1);
for i=1:K
    sigma_square(i)=Vk(i,:)*Vk(i,:)';
end
% sigma = vk'*vk;
one_vector = ones(size(sigma_square));
[R,q] = waterfilling(one_vector,sigma_square,P);

q_tem = zeros(size(a));
q_tem(a==1) = q;
q = q_tem;
    





end