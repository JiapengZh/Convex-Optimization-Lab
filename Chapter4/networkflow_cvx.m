
function [x,u,lambda,util] = networkflow_cvx(M,terminal,source)

[n_nodes,n_flows] = size(M);
n_com = length(terminal);

for com = 1:n_com
  c{com} = -M(terminal(com),:)';
  Mbar{com} = M(setdiff(1:n_nodes,[terminal(com),source(com)]),:);
end

d = zeros(n_nodes-2,1);

cvx_begin quiet
  variable util(n_com)
  variable x(n_flows,n_com)
  variable u(n_flows)  
  dual variable lambda
  maximize (  sum(util)  )
  subject to 
    for com = 1:n_com
      util(com) <= log(1 + c{com}'*x(:,com))
      Mbar{com}*x(:,com) == d 
      x(:,com) >= 0
    end
    lambda: sum(x,2) <= u
    norm(u) <= 1
    u >= 0
cvx_end 

fprintf(1, ['The optimal value found with CVX is ' num2str(sum(util),'%1.4f') '\n']);
end
