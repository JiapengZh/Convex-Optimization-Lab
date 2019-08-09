function [alpha,x,u,s,util_best] = primal_recovery(X,U,S)
% [alpha,x,u,s] = primal_recovery(X,U,S)
%
% Compute general primal recovery
%
% Input:
%  X: primal flow variables over iterations
%  U: primal capacity variables over iterations
%  S: primal 's' values over iterations
%
% Output:
%  alpha: coefficients of the optimal convex combination
%  x: feasible flow variables
%  u: feasible capacity variables
%  s: feasible 's' values
%

% TODO

iter = size(S,2);
n_com = size(S,1);
arc = size(X,1);


X_ = zeros(arc,n_com,iter+1);
X_(:,:,1) = 0;
X_(:,:,2:iter+1) = X;

S_ = zeros(n_com,iter+1);
S_(:,1) = 0;
S_(:,2:iter+1) = S;

U_ = zeros(arc,iter+1);
U_(:,1) = 0;
U_(:,2:iter+1) = U;

util = zeros(iter+1+arc,1);
util(1) = 0;
g = zeros(arc+1,iter+1+arc);
g(:,1) = 0;
for i = 2:iter+1
    
    ratesum = 0;
    for k = 1:n_com
        ratesum = ratesum + log(1+S_(k,i));
        
    end
    util(i) = ratesum;
    
    g(1:arc,i) = U_(:,i) - sum(X_(:,:,i),2);
    g(arc+1,i) = 1;
    
    
end 
g(1:arc,iter+2:end) = -eye(arc);
g(arc+1,iter+2:end) = 0;
util(iter+2:end) = 0;
b = zeros(arc+1,1);
b(arc+1) = 1;
c = util;
A = g;

[alpha,~,~] = linprog_cvx(-c,A,b);
util_best = c' * alpha;

x_matrix = zeros(arc,iter+1,n_com);
x = zeros(arc,n_com);
s = zeros(n_com,1);
for i = 1:iter+1
    for k = 1:n_com
    x_ = X_(:,k,i);
    x_matrix(:,i,k) = x_;
    %s_ = S_(k,i);
    
    end
end
for k = 1:n_com
    x(:,k) = x_matrix(:,:,k) * alpha(1:iter+1);
    s(k) = S_(k,:) * alpha(1:iter+1);
end

u = U_ * alpha(1:iter+1);


end
    
    
        