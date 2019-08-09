function R=sumrate_SIMO(H,x,y)
% function R=sumrate_SIMO(H,x,y)
%
% Compute sum-rate
%
% Inputs
%  H: m x 1 x K set of channel vectors h_1, ..., h_K
%  x: K x 1 column vector
%  y: optional K x 1 column vector
%
% Outputs
%  R: SIMO sumrate with parameters x and y [cf. (6.4)] if y is not given, y
%  is set equal to x
%

% TO DO
if exist('y','var')==0
    y = x;
end

[m,n,K] = size(H);
R = 0;
for i = 1:K
    sum_term = 0;
    if i == 1 
        for o = 2:K
            sum_term = sum_term + H(:,:,o)*y(o,:)*H(:,:,o)';
        end
    else
        for j = 1:i-1
            sum_term = sum_term + H(:,:,j)*y(j,:)*H(:,:,j)';
        end
        for w = i+1:K;
            sum_term = sum_term + H(:,:,w)*y(w,:)*H(:,:,w)';
        end
    end
    
    R = R +log2(1+x(i,:)*H(:,:,i)'*inv(eye(m)+sum_term)*H(:,:,i));
end
end