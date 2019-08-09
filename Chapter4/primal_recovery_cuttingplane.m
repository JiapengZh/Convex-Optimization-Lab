function [x,u,s] = primal_recovery_cuttingplane(X,U,S,alpha)
% [x,u,s] = primal_recovery_cuttingplane(X,U,S,alpha)
%
% Compute general primal recovery
%
% Input:
%  X: primal flow variables over iterations
%  U: primal capacity variables over iterations
%  S: primal 's' values over iterations
%  alpha: dual variables of cutting plane master problem
%
% Output:
%  x: feasible flow variables
%  u: feasible capacity variables
%  s: feasible 's' values
%

% TODO
iter = size(S,2);
n_com = size(S,1);
arc = size(X,1);


x_matrix = zeros(arc,iter,n_com);
x = zeros(arc,n_com);
s = zeros(n_com,1);

alpha_ = [alpha{:}];
alpha_ = alpha_';
for i = 1:iter
    for k = 1:n_com
    x_ = X(:,k,i);
    x_matrix(:,i,k) = x_;
    %s_ = S(k,i);
    
    end
    %alpha_(i) = alpha{i};
    
end
for k = 1:n_com
    x(:,k) = x_matrix(:,:,k) * alpha_(1:iter);
    s(k) = S(k,:) * alpha_(1:iter);
end

u = U * alpha_(1:iter);

