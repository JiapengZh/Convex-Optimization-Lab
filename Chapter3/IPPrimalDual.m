function [x y s] = IPPrimalDual(c, A, b, sigma, zeta)
% function [x y s] = IPPrimalDual(c,A,b,sigma,zeta)
%
% Linar Program solver using a primal-dual method.
%
% Inputs
% c,A,b: problem data which defines the linear program in standard form
% sigma: optional parameter defining the centering parameter (default: sigma = 0.5)
% zeta: optional parameter controlling the termination of the algorithm (default: zeta = 0.0001)
%
% Outputs
% x,y,s: resulting estimates of the optimal primal and dual variables

% TODO
m = size(A,1);
n = size(A,2);
x = ones(n,1);
y = zeros(m,1);
s = zeros(n,1);
div = c' * x - b' * (A'^(-1) * (c - s));
mu = div / it;
while mu >= zeta

    X = diag(x);
    S = diag(s);
    onev = ones(n,1);
   
    F_g = zeros(3*n,3*n);
    F_g(1:n,1:n) = zeros(n,n);
    F_g(1:n,n+1:n+m) = A';
    F_g(1:n,n+m+1:2*n+m) = eye(n);
    F_g(n+1:n+m,1:n) = A;
    F_g(n+1:n+m,n+1:n+m) = zeros(m,m);
    F_g(n+1:n+m,n+m+1:2*n+m) = zeros(m,n);
    F_g(n+m+1:2*n+m,1:n) = S;
    F_g(n+m+1:2*n+m,n+1:n+m) = zeros(n,m);
    F_g(n+m+1:2*n+m,n+m+1:2*n+m) = X;

    F_aff = zeros(2*n+m,1);
    rc = A' * y + s - c;
    rb = A * x - b;
    rs = X * S * onev;
    F_aff(1:n) = -rc;
    F_aff(n+1:n+m) = -rb;
    F_aff(n+m+1:2*n+m) = -rs;


    direc_aff = F_g^(-1) * F_aff;
    deltax_aff = direc_aff(1:n);
    deltas_aff = direc_aff(n+m+1:2*n+m);
    
    F_c = zeros(2*n+m,1);
    F_c(1:n) = zeros(n,1);
    F_c(n+1:n+m) = zeros(m,1);
    F_c(n+m+1:2*n+m) = sigma * mu * onev;
  
    direc_c = F_g^(-1) * F_c;
    deltax_c = direc_c(1:n);
    deltas_c = direc_c(n+m+1:2*n+m);
    
    deltax = deltax_aff + deltax_c;
    deltas = deltas_aff + deltas_c;
    direc = direc_aff + direc_c;
    
    minimumx = x(1)/deltax(1);
    for t = 2:n
        if deltax(t) <= 0
            minimumx = min(minimumx,x(t)/deltax(t));
        end
    end
    
    minimums = s(1)/deltas(1);
    for t = 2:n
        if deltas(t) <= 0
            minimums = min(minimums,s(t)/deltas(t));
        end
    end

    alpha = 0.99995 * min(min(minimumx,minimums),1);
    x = x + alpha * direc(1:n);
    y = y + alpha * direc(n+1:n+m);
    s = s + alpha * direc(n+m:2*n+m);
    
    
    div = c' * x - b' * y;
    mu = div / n;

end
end