function  [P,T,status] = PMinVectorBC_sdpt4(H,rho,sigma)
% Power Minimization for the vector broadcast channel using second
% order cone programming directly with SDPT3.
%
% [P,T,status] = PMinVectorBC_sdpt3(H,rho,sigma)
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

    % QoS feasibility check (under infinite power budget)
    % ***TO DO***
    if sum(2.^(-rho))>K-rank(H) & 0<2.^(-rho)<=1
    % Problem formulation with YALMIP
    % ***TO DO***
    sdpvar p;

    T = sdpvar(N,K,'full','complex');
    %T = sdpvar(N,K);

    vec_T = T(:);
    cons = [cone(vec_T,p)];
    %cons= [norm(vec_T,2) <= p];
        for i=1:K
            e=zeros(K,1);
            e(i)=1;
            XX=H*T;
            gamma = 2^rho(i)-1;
            a=sqrt(1+1/gamma)*real(XX(i,i));
            b=[T'*H'*e;sigma(i)];
            %a = [sqrt(1+1/(2^(rho(i))-1))*real(XX(i,i));T'*H'*e;sigma(i)];
            cons=[cons,cone(b,a)];
            %cons=[cons,norm(b,2) <= a]

        end

    obj=p;

    sol = optimize(cons,obj,sdpsettings('solver','sdpt3'));
    status = sol.problem;
    optobj = value(obj)
    P = optobj^2;
    T = value(T);
    % Determine P and T with YALMIP output
    % ***TO DO***

    else
        status=1;
    end
end