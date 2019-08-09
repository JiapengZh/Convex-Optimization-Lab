% [x,u,lambda,util,D] = networkflow_subgradient(M,terminal,source)
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
 


function [x,u,lambda,util,D,Util_pr,Util_pr_infeasible,iter] = networkflow_subgradient(M,terminal,source)

  [n_nodes,n_flows] = size(M);
  n_com = length(terminal);

  lambda = ones(n_flows,1);
  tic
  for iter = 1:300

    % network layer problem
    %[s,x] = nl(M,terminal,source,lambda);
    for k = 1:n_com
        [s,x] = nl(M,terminal(k),source(k),lambda);
        S(k,iter) = s;
        X(:,k,iter) = x;
    end
    % physical layer problem
    u = pl(lambda);
    U(:,iter) = u;
    % compute dual value 
    %D(iter) =  log(1+s) - lambda'*x + lambda'*u; % TODO
    %best_dual = min(D);
    ratesum = 0;
    util = zeros(iter,1);
    for k = 1:n_com
        network(k,iter) = log(1+S(k,iter)) - lambda'*X(:,k,iter);
        ratesum = ratesum + network(k,iter);
        util(iter) = util(iter) + log(1+S(k,iter));
    end
    D(iter) = ratesum + lambda'*u;
    best_dual = min(D);
    

    % subgradient update
    lambda = max(lambda - ((u - sum(X(:,:,iter),2))/iter),0);   % TODO 
    
    fprintf(1, ['Iteration: ' num2str(iter,'%u') '\t \t Best Dual: '  num2str(best_dual,'%1.4f')  '\n']);
    
    
    [~,~,~,~,util_best] = primal_recovery(X,U,S);
    Util_pr(iter) = util_best;
    
    [x,u,s] = primal_recovery_subgradient(X,U,S);
    Util_pr_infeasible(iter) = 0;
    for k = 1:n_com
        Util_pr_infeasible(iter) = Util_pr_infeasible(iter) + log(1+s(k));
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
  fprintf(1, ['Best Dual ' num2str(best_dual,'%1.4f') '\t Primal ' num2str(sum(util),'%1.4f') '\t Reconstructed ' num2str(Util_pr(end),'%1.4f') '\t Reconstructed_inf ' num2str(Util_pr_infeasible(end),'%1.4f') '\n']);

plot(1:iter, D(1:iter), 'm-',1:iter, Util_pr(1:iter),'c-',1:iter,Util_pr_infeasible(1:iter),'g-');
legend('Dual-Iter','Reconstructed-Iter','Reconstructedinf-Iter');

end


