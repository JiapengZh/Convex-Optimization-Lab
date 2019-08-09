function [R q] = srmax_gradient(H,P,q)
% function [R q] = srmax_gradient(H,P,q)
%
% Gradient optimization of the sum-rate
%
% Input
%  H: m x 1 x K set of channel vectors
%  P: maximum sum power
%  q: K x 1 column vector with the initial power allocation
%
% Output
%  R: sum rate achieved at the local optimum
%  q: power allocation at the local optimum
%

% TO DO