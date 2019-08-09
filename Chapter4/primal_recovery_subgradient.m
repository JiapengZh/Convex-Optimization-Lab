function [x,u,s] = primal_recovery_subgradient(X,U,S)
% [x,u,s] = primal_recovery_subgradient(X,U,S)
%
% Compute general primal recovery
%
% Input:
%  X: primal flow variables over iterations
%  U: primal capacity variables over iterations
%  S: primal 's' values over iterations
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

x = zeros(arc,n_com);
s = zeros(n_com,1);
u = zeros(arc,1);
for i  = 1:iter
    for k = 1:n_com
        x(:,k) = x(:,k) + X(:,k,i);
        s(k) = s(k) + S(k,i);
    end
    u = u + U(:,i);
end

x = x/iter;
s = s/iter;
u = u/iter;
end


