function grad = gradient_SIMO(H,x)
% function grad = gradient_SIMO(H,x)
%
% Gradient of the sum-rate w.r.t. the power allocation
%
% Inputs
%  H: m x 1 x K set of channel vectors
%  x: column vector containing transmit powers x_1, ..., x_K
%
% Outputs
%  grad: gradient dR(x)/dx
%
K=size(H,3);
m=size(H,1);
grad=zeros(K,1);
for j=1:K
    sum_Cj=0;
    for l=1:K
        if l~=j
            sum_Cj=sum_Cj+H(:,:,l)*x(l)*H(:,:,l)';
        end
    end
    Cj=eye(m)+sum_Cj;
    part_2=0;
    for k=1:K
        if k~=j
            sum_Ck=0;
            for l=1:K
               if l~=k
                    sum_Ck=sum_Ck+H(:,:,l)*x(l)*H(:,:,l)';
               end
            end
        
            C_k=eye(m)+sum_Ck;
            
            part_2=part_2+(x(k)*(abs(H(:,:,k)'*inv(C_k)*H(:,:,j)))^2)/(log(2)*(1+x(k)*H(:,:,k)'*inv(C_k)*H(:,:,k)));    
        end
    end
    part_1=H(:,:,j)'*inv(Cj)*H(:,:,j)/(log(2)*(1+x(j)*H(:,:,j)'*inv(Cj)*H(:,:,j)));
    grad(j)=part_1-part_2;
    
end
end
% TO DO