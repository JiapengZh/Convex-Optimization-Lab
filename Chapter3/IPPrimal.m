function [x y] = IPPrimal(c,A,b, beta, tau, zeta)
% function [x y] = IPPrimal(c,A,b,beta,tau,zeta)
%
% Linar Program solver using a primal method.
%
% Inputs
% c,A,b: problem data which defines the linear program in standard form
% beta: optional parameter that controls the progression of the algorithm (default: beta = 0.5)
% tau: optional parameter defining a starting value for tau (default: tau = 1)
% zeta: optional parameter controlling the termination of the algorithm (default: zeta = 0.0001)
%
% Outputs
% x,y: resulting estimates of the optimal primal and dual variables

% TODO
%[xopt,yopt,sopt] = linprog(c,A,b);
m = size(A,1);
n = size(A,2);
x = ones(n,1);
y = zeros(m,1);


while tau >= zeta

    X = diag(x);
    onev = ones(n,1);
    inv = X^(-2);
    i = 1;

    F_g = zeros(n+m,n+m);
    F_g(1:n,1:n) = (-tau) * inv;
    F_g(1:n,n+1:n+m) = A';
    F_g(n+1:n+m,1:n) = A;
    F_g(n+1:n+m,n+1:n+m) = zeros(m,m);

    F = zeros(n+m,1);
    F(1:n) = -(A' * y + tau * (inv^(0.5)) * onev - c);
    F(n+1:n+m) = -(A * x - b);


    direc = F_g^(-1) * F;
    deltax = direc(1:n);
    minimum = x(1)/deltax(1);
    for t = 2:n
        if deltax(t) <= 0
            minimum = min(minimum,x(t)/deltax(t));
        end
    end
    
    %stack = max(x ./ (-deltax),0);
    %[stack,~] = sort(stack,'ascend');
    %while 1
    %    if stack(i) == 0
    %        i = i + 1;
    %    else
    %        break;
    %    end
    %end
    
    %alpha = 0.99995 * min(stack(i),1);
    
    alpha = 0.99995 * min(minimum,1);
    x = x + alpha * direc(1:n);
    y = y + alpha * direc(n+1:n+m);

    tau = beta * tau;

end
end



