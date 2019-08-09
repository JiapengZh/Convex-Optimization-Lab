function [R,s,mu] = worstCaseRate(p,N)
% Function: [R,s,mu] = worstcaseRate(p,N)
%
% Calculation of the worst case rate to the sum rate 
% minimization problem under a limited sum noise power.
%     R =    min sum(log2(1+psi_i/sigma_i) 
%            s.t.:   sigma_i >= 0   for all i
%                    sum(sigma_i) <= N
% The sulution tried to be calculated via CVX.
% Note that the psi_i has to be strictly positive.
%
% Inputs
% p: vector of useful received power psi_1,...,psi_m
% N: available sum noise power N
%
% Outputs
% R: value of the worst-case sum rate R
% s: vector of worst-case powers sigma_1,...,sigma_m
% mu: value of the optimal Lagrangian multiplier mu


% TODO
p=p(:);
mu_p_high=((2*N+sum(p))^2-sum(p.^2))/(4*sum(p));

disp(mu_p_high);
fun=@(x) -2*N-sum(p)+sum((p.^(2)+4*p*x).^(1/2));

mu_p=fzero(fun,[0 mu_p_high]);
disp(mu_p);
mu=1/(log(2)*mu_p);
s=zeros(size(p));
for i=1:size(p)
   if p(i)==0
       s(i)=0;
   elseif p(i)>0
       s(i)=-p(i)/2+1/2*(p(i)^2+4*p(i)*mu_p)^(1/2);
   end
end
R=sum(log(1+p./s))/log(2);
end