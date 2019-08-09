function rate = pl_mimo(M,P,weights)

  [n_nodes,n_flows] = size(M,1);
  RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  H = randn(2,2,n_flows) + 1i*randn(2,2,n_flows); 
  

  rate = zeros(n_flowss,1);
  % for every node call the MIMO MAC function
  for node = 1:n_nodes
    % all incoming links to the node
    links = find(M(node,:) == -1);
    r = MIMO_MAC_cvx(H(:,:,links),P,weights(links));
    rate(links) = rate(links) + r(:)./n_nodes;
  end
  rate(weights < 0.0001) = 0;
end
