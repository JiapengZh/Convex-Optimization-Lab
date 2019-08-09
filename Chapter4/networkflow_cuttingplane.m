% [x,u,lambda,util,D,LB,UB] = networkflow_cuttingplane(M,terminal,source)
%
% network flow maximization, dual decomposition, subgradient
%
% Input
%  M: node-arc incidence matrix
%  terminal: index of the terminal
%  source: index of the source
% 
% Output
%
% x: flow solution  
% u: capacity solution  
% lambda: dual variables  
% util: utility function value  
% D: dual function value over iterations
% LB: lower bounds over iterations for utility function   
% UB: upper bounds over iterations for utility function

function [x,u,lambda,util,D,LB,UB,Util_pr,Util_cu,iter] = networkflow_cuttingplane(M,terminal,source)

  [n_nodes,n_flows] = size(M);  
  n_com = length(terminal);

  lambda = ones(n_flows,1);
  
  epsilon = 0.01;
  
  tic
  
  for iter = 1:1000

    % network layer problem
    
    for k = 1:n_com
        [s,x] = nl(M,terminal(k),source(k),lambda);
        S(k,iter) = s;
        X(:,k,iter) = x;
    end
    
    
    % physical layer problem
    u = pl(lambda);    
    
    % collect results
    %S(iter) = s;  
    %X(:,iter) = x;
    U(:,iter) = u;

    % compute dual value 
    
    %D(iter) =  log(1+s) - lambda'*x + lambda'*u;  % TODO
    %[best_dual,best_dual_idx] = min(D);
    %lambda_v(:,iter) = lambda;
    ratesum = 0;
    util = zeros(iter,1);
    for k = 1:n_com
        network(k,iter) = log(1+S(k,iter)) - lambda'*X(:,k,iter);
        ratesum = ratesum + network(k,iter);
        util(iter) = util(iter) + log(1+S(k,iter));
    end
    D(iter) = ratesum + lambda'*u;
    best_dual = min(D);
    lambda_v(:,iter) = lambda;

    % cutting plane master problem (hint: use S,X and U)
    
    cvx_begin
        variable r
        variable lmd(n_flows)
        dual variables alpha{iter}
        minimize r
        subject to
            lmd >= 0
            for i = 1:iter
                
                alpha{i} : r >= D(i) + (U(:,i) - sum(X(:,:,i),2))' * (lmd - lambda_v(:,i))
                
            end
            
    cvx_end
    
    
    lambda = lmd; % TODO
    % update upper and lower bound
    LB(iter) = r; % TODO
    UB(iter) = 
    ; % TODO
    
     
    % convergence criterion 
    criterion = (UB(iter) - LB(iter)) / abs(LB(iter)); % TODO
    
    fprintf(1, ['Iteration: ' num2str(iter,'%u') '\t \t Lower Bound: '  num2str(LB(iter),'%1.4f') '\t \t Upper Bound: '  num2str(UB(iter),'%1.4f') '\t \t convergence criterion: '  num2str(criterion,'%1.4f') '\n']);
    
    [~,~,~,~,util_best] = primal_recovery(X,U,S);
    Util_pr(iter) = util_best;
    
    [x,u,s] = primal_recovery_cuttingplane(X,U,S,alpha);
    Util_cu(iter) = 0;
    for k = 1:n_com
        Util_cu(iter) = Util_cu(iter) + log(1+s(k));
    end
    
    if (criterion < epsilon )
      break
    end

  end
  toc
  % Recover optimal primal values
  for com = 1:n_com
    c{com} = -M(terminal(com),:)';
    Mbar{com} = M(setdiff(1:n_nodes,[terminal(com),source(com)]),:);
  end
  d = zeros(n_nodes-2,1);
  cvx_begin quiet
    variable util(n_com)
    variable x(n_flows,n_com)
    maximize (  sum(util)  )
    subject to 
      for com = 1:n_com
        util(com) <= log(1 + c{com}'*x(:,com))
        Mbar{com}*x(:,com) == d 
        x(:,com) >= 0
      end
      sum(x,2) <= u
  cvx_end 

 
  fprintf(1, '\n ===========================================================\n');
  fprintf(1, ['Best Dual ' num2str(best_dual,'%1.4f') '\t Primal ' num2str(sum(util),'%1.4f') '\t Reconstructed ' num2str(Util_pr(end),'%1.4f') '\t Reconstructed_cu ' num2str(Util_cu(end),'%1.4f') '\n']);

  
plot(1:iter, D(1:iter), 'm-',1:iter,LB(1:iter),'c-',1:iter,UB(1:iter),'g-',1:iter, Util_pr(1:iter),'b-',1:iter, Util_cu(1:iter),'y-');
legend('Dual-Iter','Lower Bound','Upper Bound','Reconstructed-Iter','Reconstructed_cu-Iter');
  
  
end



