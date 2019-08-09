function [x y s] = IPPredictorCorrector(c,A,b,zeta)
% function [x y s] = IPPredictorCorrector(c,A,b,zeta)
%
% Linar Program solver using a primal-dual predictor-corrector interior point method.
%
% Inputs
% c,A,b: problem data which defines the linear program in standard form
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
mu = div / n;

while mu >= zeta

    X = diag(x);
    S = diag(s);
    onev = ones(n,1);

    F_g = zeros(2*n+m,2*n+m);
    F_g(1:n,1:n) = zeros(n,n);
    F_g(1:n,n+1:n+m) = A';
    F_g(1:n,n+m+1:2*n+m) = eye(n);
    F_g(n+1:n+m,1:n) = A;
    F_g(n+1:n+m,n+1:n+m) = zeros(m,m);
    F_g(n+1:n+m,n+m+1:2*n+m) = zeros(n,n);
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
    deltay_aff = direc_aff(n+1:n+m);
    deltas_aff = direc_aff(n+m+1:2*n+m);
    X_aff = diag(deltax_aff);
    S_aff = diag(deltas_aff);
    
    minimumx_aff = x(1)/deltax_aff(1);
    for t = 2:n
        if deltax_aff(t) <= 0
            minimumx_aff = min(minimumx_aff,x(t)/deltax_aff(t));
        end
    end
    
    alpha_aff = 0.99995 * min(min(minimumx_aff,minimums_aff),1);
    div_aff = c' * (x + alpha_aff * deltax_aff) - b' * ((A'^(-1)) * (c - s - alpha_aff * deltas_aff));
    mu_aff = div_aff / n;
    sigma = (mu_aff/mu)^3;
    F_cc = zeros(2*n+m,1);
    F_cc(1:n) = zeros(n,1);
    F_cc(n+1:n+m) = zeros(m,1);
    F_cc(n+m+1:2*n+m) = sigma * mu * onev - X_aff * S_aff * onev;
  
    direc_cc = F_g^(-1) * F_cc;
    deltax_cc = direc_cc(1:n);
    deltay_cc = direc_cc(n+1:n+m);
    deltas_cc = direc_cc(n+m+1:2*n+m);
    
    deltax = deltax_aff + deltax_cc;
    deltas = deltas_aff + deltas_cc;
    deltay = deltay_aff + deltay_cc;
    direc = direc_aff + direc_cc;
    
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

    x = x + alpha * deltax;
    y = y + alpha * deltay;
    s = s + alpha * deltas;
    
    div = c' * x - b' * (A'^(-1) * (c - s));
    mu = div / n;

end
end