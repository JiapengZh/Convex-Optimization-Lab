function [R q] = srmax_zf_exhaustive(H,P)
% function [R q] = srmax_zf_exhaustive(H,P)
%
% Exhaustive optimization of the zero-forcing sum-rate
%
% Input
%  H: m x 1 x K set of channel vectors
%  P: maximum sum power
%
% Output
%  R: sum rate achieved at the globally optimal zero-forcing solution
%  q: power allocation at the globally optimal zero-forcing solution
%

% TO DO
[m,~,K]=size(H);
R_max = 0;
for num = 1:2^K-1
    
    value=num2str(dec2bin(num));
    value=[repmat('0',[1,K-length(value)]),value];

    s=zeros(K,1);
    %for j=1:value
    for i=1:K
        s(i)=str2num(value(i));
    end
    
    if length(s(s==1))>m
        continue;
    end
    
    [R_tem,q_tem]=sumrate_zf(H,P,s);
    if R_tem > R_max
        R_max = R_tem;
        q_max = q_tem;
    end
end
R = R_max;
q = q_max;
end
