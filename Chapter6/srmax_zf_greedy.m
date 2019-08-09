function [R q] = srmax_zf_greedy(H,P)
% function [R q] = srmax_zf_greedy(H,P)
%
% Greedy optimization of the zero-forcing sum-rate
%
% Input
%  H: m x 1 x K set of channel vectors
%  P: maximum sum power
%
% Output
%  R: sum rate achieved at the suboptimal zero-forcing solution
%  q: power allocation at the suboptimal zero-forcing solution
%

% TO DO
[m,~,K] = size(H);
s = zeros(K,1);
K_da=[];
R_max = 0;
q_max = 0;
s_tem = zeros(K,1);

for i = 1:m 
    if i > K
        break;
    end
    for k = 1:K
        
        if ismember(k,K_da)
            if k == K
                s_tem = zeros(K,1);
                for l = 1:length(K_da)
                    s_tem(K_da(l)) = 1;
                end
            end
            continue;
        end
        
        s = s_tem;
        s(k,1) = 1;
        [R q] = sumrate_zf(H,P,s);
        if R > R_max
            R_max = R;
            q_max = q;
            K_da(i,1) = k;
        end
        if k == K
            s_tem = zeros(K,1);
            for l = 1:length(K_da)
                s_tem(K_da(l)) = 1;
            end
        end
        
    end
    s = s_tem;
end
    R = R_max;
    q = q_max;
end
    