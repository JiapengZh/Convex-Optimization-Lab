function [Q1,Q2] = wsr_mimo_mac_yalmip(H1,H2,P1,P2,w)
% FUNCTION           wsr_mimo_mac_yalmip.m
% SYNTAX   [Q1,Q2] = wsr_mimo_mac_yalmip(H1,H2,P1,P2,w)
%
% EXPLANATION
% Weighted sum-rate maximization in the MIMO MAC using the conic
% optimization solver SDPT3.
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
    % ***TO DO***
    alpha1 = (w(1)-w(2))/2;
    alpha2 = w(2)/2;
    Q1 = sdpvar(n1,n1);
    Q2 = sdpvar(n2,n2);
%     term1 = eye(m1)+H1*Q1*H1';
%     term2 = eye(m1)+H1*Q1*H1'+H2*Q2*H2';
%     X1 = sdpvar(m1,m1);
%     X2 = sdpvar(m2,m2);
    %cons=[trace(Q1)<=P1,trace(Q2)<=P2,X1-H1*Q1*H1'==eye(m1),X2-H1*Q1*H1'-H2*Q2*H2'==eye(m1)];
    cons=[Q1>=0,Q2>=0,trace(Q1)<=P1,trace(Q2)<=P2];
    %obj= -alpha1*logdet(eye(m1)+H1*Q1*H1')./log(2)-alpha2*logdet(eye(m1)+H1*Q1*H1'+H2*Q2*H2')./log(2);
    %obj = -alpha1/log(2)*logdet(term1)-alpha2/log(2)*logdet(term2);
    obj = -alpha1/log(2)*logdet(eye(m1)+H1*Q1*H1')-alpha2/log(2)*logdet(eye(m1)+H1*Q1*H1'+H2*Q2*H2');
    %obj= -alpha1/log(2)*logdet(X1)-alpha2/log(2)*logdet(X2);
    sol = optimize(cons,obj,sdpsettings('solver','sdpt3'));
    % use YALMIP to solve problem
    % ***TO DO***
    
    
% case 2: w1 < w2
elseif ( w(1) < w(2) )
    
    % define effective weight vector alpha
    % ***TO DO***
    alpha1 = (w(2)-w(1))/2;
    alpha2 = w(1)/2;
    Q1 = sdpvar(n1,n1);
    Q2 = sdpvar(n2,n2);
%     term1 = eye(m1)+H2*Q2*H2';
%     term2 = eye(m1)+H1*Q1*H1'+H2*Q2*H2';
    cons=[Q1>=0,Q2>=0,trace(Q1)<=P1,trace(Q2)<=P2];
%     X1 = sdpvar(m1,m1);
%     X2 = sdpvar(m2,m2);
%     cons=[trace(Q1)<=P1,trace(Q2)<=P2,X1-H2*Q2*H2'==eye(m1),X2-H1*Q1*H1'-H2*Q2*H2'==eye(m1)];
    %obj= -alpha1*logdet(eye(m1)+H2*Q2*H2')./log(2)-alpha2*logdet(eye(m1)+H2*Q2*H2'+H1*Q1*H1')./log(2);
    %obj = -alpha1*logdet(Q2)-alpha2*logdet(Q2);
    %obj= -alpha1/log(2)*logdet(term1)-alpha2/log(2)*logdet(term2);
    obj= -alpha1/log(2)*logdet(eye(m1)+H2*Q2*H2')-alpha2/log(2)*logdet(eye(m1)+H1*Q1*H1'+H2*Q2*H2');
    sol = optimize(cons,obj,sdpsettings('solver','sdpt3'));
    % use YALMIP to solve problem
    % ***TO DO***
    
% case 3: w1 = w2
else
    
    % use YALMIP to solve problem
    % ***TO DO***
    alpha2 = w(1)/2;
    Q1 = sdpvar(n1,n1);
    Q2 = sdpvar(n2,n2);
%     term = eye(m1)+H1*Q1*H1'+H2*Q2*H2';
%     X2 = sdpvar(m1,m1);
%     cons=[trace(Q1)<=P1,trace(Q2)<=P2,X2-H1*Q1*H1'-H2*Q2*H2'==eye(m1)];
    %obj= -alpha2*logdet(eye(m1)+H1*Q1*H1'+H2*Q2*H2')./log(2);
    %obj = -alpha2*logdet(Q1);
%     obj= -alpha2/log(2)*logdet(term);
    cons=[Q1>=0,Q2>=0,trace(Q1)<=P1,trace(Q2)<=P2];
    obj= -alpha2/log(2)*logdet(eye(m1)+H1*Q1*H1'+H2*Q2*H2');
    sol = optimize(cons,obj,sdpsettings('solver','sdpt3'));
end
Q1 = value(Q1);
Q2 = value(Q2);
% done
