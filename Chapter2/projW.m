function WP = projW(W, P)
% function WP = projW(W, P)
%
% Projection of the precoders in W onto the set of precoders with transmit 
% power lower or equal to P
%
% Inputs
% W: n x n x K set of precoders
% P: available sum transmit power
% Outputs
%WP: n x n x K set of projected precoders
n=size(W,1);
K=size(W,3);
WP=zeros(n,n,K);
summ=0;
for l = 1:K
    W_l=W(:,:,l);
    f_norm = norm(W_l,'fro');
    %f_norm_sq= sum(diag(W_l'*W_l));
    summ = summ+f_norm^2;
end
%if summ <= P
    %for i = 1:K
        %WP(:,:,i)=W(:,:,i);
    %end
%else
    for i = 1:K
        WP(:,:,i) = (P/summ)^0.5*W(:,:,i);
    end        
%end
    
% TODO
