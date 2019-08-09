function R = c_mimo_mac(H1,H2,P1,P2,W)
% FUNCTION           c_mimo_mac.m
% SYNTAX   [R] = c_mimo_mac(H1,H2,P1,P2,W)
%
% EXPLANATION
% Weighted sum-rate maximization in the MIMO MAC using the Matlab
% function WSR_MIMO_MAC.m and various wheigts collected in W.
% 
% Solves the following optimization problem vor various w1 and w2:
%               maximize      w1*R1 + w2*R2
%               subject to    trace(Q1) <= P1
%                             trace(Q1) <= P1
%               variables     Q1, Q2
%
% INPUTS 
% H1    : channel gain matrix H1
% H2    : channel gain matrix H2
% P1    : available transmit power at user 1
% P2    : available transmit power at user 2
% W     : 2 x L matrix of weight vectors W = [w1,...,wL]
%
% OUTPUTS 
% R     : 2 x L matrix of optimal rate vectors R = [r1,...,rL]
%******************************************************************

%% PLAUSIBILITY CHECK

% check if channel matrix dimensions match
[m1,~] = size(H1);
[m2,~] = size(H2);
if ( m1~=m2 )
    error('Dimensions of channel matrices do not match!');
else
    m = m1;
end

% check if power constraints are valid
if ( (P1<=0) || (P2<=0) )
    error('P1 and P2 must be positive!')
end

% check if matrix of weight vectors is valid
[Z,L] = size(W);
if ( Z ~= 2 )
    error('Weight vector must be a matrix of dimension 2xL!')
elseif ( (W>=0) ~= ones(2,L) )
    error('Weight vector must have nonnegative entries!')
end

%% SOLVE WSR MAXIMIZATION PROBLEM FOR ALL L WEIGHT VECTORS

% define matrix R of dimension 2xL to store rate vectors obtained as the
% solutions of the L WSR maximization problems
R = zeros(2,L);

% ***TO DO***
for i=1:L
    [Q1,Q2] = wsr_mimo_mac_yalmip(H1,H2,P1,P2,W(:,i));
    if W(1,i)>=W(2,i)
        R(:,i)=[0.5*logdet(eye(m1)+H1*Q1*H1')./log(2);0.5*logdet(eye(m2)+H1*Q1*H1'+H2*Q2*H2')./log(2)-0.5*logdet(eye(m1)+H1*Q1*H1')./log(2)];
    else
        R(:,i)=[0.5*logdet(eye(m2)+H1*Q1*H1'+H2*Q2*H2')./log(2)-0.5*logdet(eye(m1)+H2*Q2*H2')./log(2);0.5*logdet(eye(m1)+H2*Q2*H2')./log(2)];
    end
    
end
r1_max = max(R(1,:));
r2_max = max(R(2,:));
zero_matrix = zeros(2,1);
r1_max_matrix = [r1_max;0];
r2_max_matrix = [0;r2_max];
matrix = [R zero_matrix r1_max_matrix r2_max_matrix]';
x=matrix(:,1);
y=matrix(:,2);
%K = convhull(matrix);
K = convhull(x,y);
plot(matrix(K,1),matrix(K,2),'r-',matrix(:,1),matrix(:,2),'bx')
end