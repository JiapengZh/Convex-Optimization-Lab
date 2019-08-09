function  [P,T,status] = PMinVectorBC_cvx(H,rho,sigma)
% Power Minimization for the vector broadcast channel using second
% order cone programming with CVX and SDPT3.
%
% [P,T,status] = PMinVectorBC_cvx(H,rho,sigma)
%
% Solves the following optimization problem:
%               minimize      P
%               subject to    norm(T(:))^2 <= P
%                             rate(T,sigma) >= rho
%               variables     P, T
%
% The reformulation of the optimization problem into a second order
% cone probelm according to [WieselTSP2006] is used.
%
% Inputs: 
%       H     : K times N matrix of channel vectors
%       rho   : minimum rate target vector
%       sigma : vector of noise variances
%
% Outputs:
%       P     : minimum transmit power
%       T     : optimal transmit precoder
%       status: status information of the optimization
%               (e.g.,'Problem','Infeasible','Solved', ...)
%
%******************************************************************
    
    % Input verification
    K = size(H,1);
    N = size(H,2);
    rho = rho(:);
    sigma = sigma(:);

    % Output initialization
    P = Inf;
    T = NaN(N,K);
    status = 'Problem';

    % QoS feasibility check (under infinite power budget)
    % ***TO DO***
    if sum(2.^(-rho))>K-rank(H) & 0<2.^(-rho)<=1

    
    % Parameter reformulation for the SOCP with SDPT3.
    % ***TO DO***
    
    % Problem formulation with CVX
    cvx_clear;
    
    cvx_begin
    cvx_solver sdpt3;
    % ***TO DO***
    variable  p;
    variable T(N,K) complex;
    minimize( p );
    subject to
       %a = [P;T(:)];
       %{a,0} <In> complex_lorentz(N*K+1);
       %norm(T(:),2)<=P;
         vec_T = T(:);
        {vec_T,p} <In> complex_lorentz(N*K);
        for i=1:K
            e=zeros(K,1);
            e(i)=1;
            XX=H*T;
            gamma = 2^rho(i)-1;
            a=sqrt(1+1/gamma)*real(XX(i,i));
            b=[T'*H'*e;sigma(i)];
            %norm(b,2)<=a;
            %b = [sqrt(1+1/(2^(rho(i))-1))*real(XX(i,i));T'*H'*e;sigma(i)];
            %{ b, 0 } <In> complex_lorentz(K+2);
            { b, a } <In> complex_lorentz(K+1);
            imag(XX(i,i))==0;
          
        end
        
    cvx_end
    P = value(p)^2;
    % Determine P and T with CVX output
    % ***TO DO***
    status = 0;
           
    else
         status = 1;
    end
end