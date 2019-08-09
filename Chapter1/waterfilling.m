function [R,p,mu] = waterfilling(h,s,P)
% function [R,p,mu] = waterfilling(h,s,P)
%
% Calculation of the optimal solution to the following sum rate
% maximization problem under a limited sum power budget.
%     R =    max sum(log2(1+|hi|^2/si*pi) 
%            s.t.:   pi >= 0   for all i
%                    sum(pi) <= P
% The sulution is calculated via the waterfilling procedure.
%
% Inputs
% h: vector of channel coefficients h1,...,hn
% s: vector of noise variances sigma1^2,...,sigma^2
% P: available sum transmit power P
% Outputs
% R: value of the maximal sum rate R
% p: vector of optimal transmit powers p1,...,pn
% mu: value of th optimal Lagrangian multiplier mu
tic
h = h(:);
s = s(:);

% determine n
n = length(h);
gamma=abs(h).^2./s;
[gamma_sorted,index] = sort(gamma,'descend');
c=1./gamma_sorted;
for k_star = 1:n
    w_star=(P+sum(c(1:k_star)))/k_star;
    if k_star < n
        if c(k_star + 1) >= w_star
            if c(k_star) > w_star
                k_star = k_star - 1;
                w_star=(P+sum(c(1:k_star)))/k_star;
                break;
            else
                break;
            end
        end
    else
        if c(k_star) > w_star
            k_star = k_star - 1;
            w_star=(P+sum(c(1:k_star)))/k_star;
            break;
        else
            break;
        end
    end
end
mu = 1/(log(2)*w_star);
p_sorted = zeros(n,1);
for i=[1:k_star]
    p_sorted(i)=w_star-c(i);
end
p = zeros(size(p_sorted));
p(index)=p_sorted;
R=sum((log2(1+gamma.*p)));
toc
end
% TODO