function [Q1,Q2] = wsr_mimo_mac_cvx(H1,H2,P1,P2,w)
% FUNCTION           wsr_mimo_mac.m
% SYNTAX   [Q1,Q2] = wsr_mimo_mac(H1,H2,P1,P2,w)
%
% EXPLANATION
% Weighted sum-rate maximization in the MIMO MAC using the conic
% optimization toolbox CVX with solver SDPT3.
% 
% Solves the following optimization problem:
%               maximize      w1*R1 + w2*R2
%               subject to    trace(Q1) = P1
%                             trace(Q2) = P2
%               variables     Q1, Q2
%
% INPUTS 
% H1    : channel gain matrix H1
% H2    : channel gain matrix H2
% P1    : available transmit power at user 1
% P2    : available transmit power at user 2
% w     : weight vector w = [w1;w2]
%
% OUTPUTS
% Q1    : optimal transmit covariance matrix Q1 of user 1
% Q2    : optimal transmit covariance matrix Q2 of user 2
%******************************************************************

%% PLAUSIBILITY CHECK

% check if channel matrix dimensions match
[m1,n1] = size(H1);
[m2,n2] = size(H2);
if ( m1~=m2 )
    error('Dimensions of channel matrices do not match!');
else
    m = m1;
end

% check if power constraints are valid
if ( (P1<=0) || (P2<=0) )
    error('P1 and P2 must be positive!')
end

% check if weight vector is valid
[x,y] = size(w);
if ( (x==1) && (y==2) )
    w = w';
elseif ( (x~=2) || (y~=1) )
    error('Weight vector must be a vector of dimension 2!')
end
if ( (w(1)<0) || (w(2)<0) )
    error('Weight vector must have nonnegative entries!')
end

%% SET UP WSR MAXIMIZATION PROBLEM

% case dinstinction
% 1: w1 > w2
% 2: w1 < w2
% 3: w1 = w2

% case 1: w1 > w2
if ( w(1) > w(2) )
    
    % define effective weight vector alpha
    %  ***TO DO***
    % use CVX to solve problem
    cvx_begin sdp
    variable Q1(n1,n1) semidefinite
    variable Q2(n2,n2) semidefinite
    cvx_solver('sdpt3')
    cvx_expert('true')
    cvx_quiet('true')
    obj= (w(1)-w(2))/2*(log_det(eye(m1)+H1*Q1*H1')./log(2))+w(2)/2*log_det(eye(m1)+H1*Q1*H1'+H2*Q2*H2')./log(2);
    % ***TO DO***
    maximize(obj)
    subject to
        trace(Q1)<=P1;
        trace(Q2)<=P2;
    cvx_end
    
% case 2: w1 < w2
elseif ( w(1) < w(2) )
    
    % define effective weight vector alpha
    % ***TO DO***
    
    cvx_begin sdp
    variable Q1(n1,n1) semidefinite
    variable Q2(n2,n2) semidefinite
    cvx_solver('sdpt3')
    cvx_expert('true')
    cvx_quiet('true')
    obj= (w(2)-w(1))/2*(log_det(eye(m2)+H2*Q2*H2')./log(2))+w(1)/2*log_det(eye(m2)+H2*Q2*H2'+H1*Q1*H1')./log(2);
    % ***TO DO***
    maximize(obj)
    subject to
        trace(Q1)<=P1;
        trace(Q2)<=P2;
    cvx_end
    
    
% case 3: w1 = w2
else
    
    % use CVX to solve problem
    cvx_begin sdp
    variable Q1(n1,n1) semidefinite
    variable Q2(n2,n2) semidefinite
    cvx_solver('sdpt3')
    cvx_expert('true')
    cvx_quiet('true')
    obj= 0.5*log_det(eye(m1)+H1*Q1*H1'+H2*Q2*H2')./log(2);
    % ***TO DO***
    maximize(obj)
    subject to
        trace(Q1)<=P1;
        trace(Q2)<=P2;
    cvx_end
    
    
end
WSR = cvx_optval;

% done