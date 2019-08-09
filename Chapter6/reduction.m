function [a,b] = reduction(P,a,b)
% function [a,b] = reduction(P,a,b)
%
% Inputs
%  P: maximum sum power
%  a: column vector specifiying the lower boundary of the box B
%  b: column vector specifiying the upper boundary of the box B
%
% Outputs
%  
%  a: column vector specifiying the lower boundary of the reduced box B'
%  b: column vector specifiying the upper boundary of the reduced box B'
%

% TO DO
if isempty(a) == 0
    [m,n] = size(a);
    sum_term = a+(P-ones(m,1)'*a)*ones(m,1);
    matrix = [sum_term b];
    b = min(matrix,[],2);
end


end