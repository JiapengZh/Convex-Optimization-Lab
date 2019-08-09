function rate = MIMO_MAC_cvx(H,P,weights)

    
  [n_tx, n_rx, n_users] = size(H);
  
  if sum(weights) == 0;
     rate = zeros(n_users,1);
     return
  end
 
  weights = weights/sum(weights);
  weights(weights < 0) = 0;
  if n_users ~= length(weights)
    error('Size of the channels does not match size of weights')
  end

  weights = weights(:);
  % find appropriate sorting of the users
  [muh,order] = sort(weights,'descend');

  % generate the difference between subsequent weights
  alphas = [-diff(muh);muh(n_users)];


  ind = zeros(n_users);
  for user = order(:)'
    ind(user,user:end) = 1;
  end
  
  cvx_begin sdp quiet
    cvx_solver sdpt3
    cvx_precision high
    variable X(n_tx,n_tx,n_users) hermitian
    variable Y(n_tx,n_tx,n_users,n_users) hermitian
    variable c(n_users)
    variable S(n_rx,n_rx,n_users) hermitian
    maximize ( sum(c) )
    subject to
      for user = 1:n_users
        c(user) <= alphas(user)*log_det(X(:,:,user))
        for i_user = 1:n_users
          Y(:,:,i_user,user) == ind(user,i_user) * H(:,:,i_user) * S(:,:,i_user) * H(:,:,i_user)'
        end
        X(:,:,user) == eye(n_tx) + sum(Y(:,:,:,user),3)
        X(:,:,user) >= 0;
        S(:,:,user) >= 0;
        trace(S(:,:,user)) <= P
      end
    
    c >= 0
  cvx_end

  rate = zeros(n_users,1);
  for user = 1:n_users
    rate(user) = log2( real( det(X(:,:,user)) ) / real( det(X(:,:,user) - H(:,:,user) * S(:,:,user) * H(:,:,user)') ));
  end
end
