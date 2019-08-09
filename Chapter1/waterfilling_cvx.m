function [R,p,mu] = waterfilling_cvx(h,s,P)
% function [R,p,mu] = waterfilling_cvx(h,s,P)
%
% Calculation of the optimal solution to the following sum rate
% maximization problem under a limited sum power budget.
%     R =    max sum(log2(1+|hi|^2/si*pi) 
%            s.t.:   pi >= 0   for all i
%                    sum(pi) <= P
% The sulution is calculated via CVX.
%
% Inputs
% h: vector of channel coefficients h1,...,hn
% s: vector of noise variances sigma1^2,...,sigma^2
% P: available sum transmit power P
% Outputs
% R: value of the maximal sum rate R
% p: vector of optimal transmit powers p1,...,pn
% mu: value of th optimal Lagrangian multiplier mu

h = h(:);
s = s(:);
gama = abs(h).^2./s;
R=0;
% determine n
n = length(h);
cvx_begin
    variable p(n);
    dual variable mu;
    R = sum(log(1+gama.*p)./log(2));
    maximize(R);
    subject to
        mu: sum(p) <= P;
        p >= 0;
cvx_end
mu.*(P-sum(p));
% TODO
end
