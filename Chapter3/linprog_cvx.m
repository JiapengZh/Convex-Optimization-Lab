function [x,y,s] = linprog_cvx(c,A,b)
% function x = linprog_cvx(c,A,b)
%
% Solve a linear program in standard form with CVX
%
% Inputs
% c,A,b: problem data which defines the linear program in standard form
%
% Outputs:
% x,y,s: resulting estimates of the optimal primal and dual variables

% TODO
m = size(A,1);
n = size(A,2);
cvx_begin
    variable x(n);
    dual variable y;
    dual variable s;
    fx = c' * x;
    minimize fx;
    subject to
        -x <= 0 : s;
        A * x == b : y;
        
cvx_end
end