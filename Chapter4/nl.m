% function [s x] = nl(M,terminal,source,lambda)
%
% Verification for network layer subproblem 
%
% Input
% M: node-arc incidence matrix
% terminal: index of the terminal
% source: index of the source
% lambda: dual variables
%
% Output
%  s,x: resulting throughput and networkflow

function [s x] = nl(M,terminal,source,lambda)

  lambda = max(lambda(:),0);
  [n_nodes,n_flows] = size(M);
  c = M(terminal,:);
  c = -c(:);
  Mbar = M(setdiff(1:n_nodes,[terminal,source]),:);
  d = zeros(n_nodes-2,1);

  cvx_clear
  cvx_begin quiet
    variable x_lin(n_flows)
    maximize (  (c - lambda)'*x_lin  )
    subject to 
      Mbar*x_lin == d
      x_lin >= 0
      c'*x_lin == 10;
  cvx_end
  cvx_clear
  % we now know the flow, scale it appropriately 
  x_lin = max(x_lin,0);
  scale = (c'*x_lin - lambda'*x_lin) / (c'*x_lin * lambda'*x_lin);

  scale = max(scale,0);

  x_lin = min(scale, 10/(c'*x_lin)).*x_lin;

  if( c'*x_lin > 10.0001)
    display('This should not happen')
  end

  % dirty hack
  if sum(x_lin) == 0
    x_lin = ones(size(x_lin))* 10^-11;
  end

  x = double(x_lin);
  %s = double(log(1+ c'*x_lin));
  s = c'*x;
end
