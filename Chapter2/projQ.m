function QP = projQ(Q,P)
% function QP = projQ(Q, P)
%
% Projection of the covariance matrices in Q onto the set of covariance 
% matrices with transmit power lower or equal to P
% 
% Inputs
% Q: n x n x K set of covariance matrices
% P: available sum transmit power
% Outputs
% QP: n x n x K set of projected covariance matrices

K = size(Q,3);
n = size(Q,1);
U = zeros(n,n,K);
QP = zeros(n,n,K);
lambda = zeros(n*K,1);
for i = 1:K
    [U(:,:,i),mm] = eig(Q(:,:,i));
    lambda(n*(i-1)+1:n*i) = real(diag(mm));
end

lambda = max(lambda,0);
if sum(lambda)<P
    for i=1:K
        QP(:,:,i) = U(:,:,i)*diag((lambda((i-1)*n+1:i*n)))*U(:,:,i)';
    end
    return
end
[lambda_sorted,index] = sort(lambda);
for i = 1:K*n
    omega_star = lambda_sorted(i);
    if sum(max(lambda_sorted-omega_star,0))<P
        omega_star = (sum(lambda_sorted(i:end))-P)/(K*n-i+1);
        lambda_sorted = max(lambda_sorted-omega_star,0);
        break;
    end
end

lambda(index) = lambda_sorted;
for i = 1:K
    QP(:,:,i) = U(:,:,i)*diag((lambda((i-1)*n+1:i*n)))*U(:,:,i)';
end