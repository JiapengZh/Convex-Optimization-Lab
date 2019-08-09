function plotObjValues(H,P,nIter)

if nargin < 1
    m = 8;
    n = 4;
    K = 4;
    
    H = 1/sqrt(2)*(randn(m,n,K) + 1i*randn(m,n,K));
    nIter = 20;
    P = 1;
end

K = size(H,3);
n = size(H,2);

W0 = zeros(n,n,K);
Q0 = zeros(n,n,K);

for k = 1:K
    W0(:,:,k) = eye(n);
    Q0(:,:,k) = eye(n);
end

%TODO
%default = 'fixed';
gradw = @(W)gradW(W,H);
projw = @(W)projW(W,P);
ratew = @(W)rateW(W,H);
gradq = @(Q)gradQ(Q,H);
projq = @(Q)projQ(Q,P);
rateq = @(Q)rateQ(Q,H);
[x_op,val_op,value_op,count_op] = projGrad(ratew,gradw,projw,W0,nIter,3);
[x,val,value,count] = projGrad(ratew,gradw,projw,W0,nIter);
[x_optimal,val_optimal,value_optimal,count_optimal] = projGrad(ratew,gradw,projw,W0,nIter,1);
[x_armijo,val_armijo,value_armijo,count_armijo] = projGrad(ratew,gradw,projw,W0,nIter,2);

[x_op_,val_op_,value_op_,count_op_] = projGrad(rateq,gradq,projq,Q0,nIter,3);
[x_,val_,value_,count_] = projGrad(rateq,gradq,projq,Q0,nIter);
[x_optimal_,val_optimal_,value_optimal_,count_optimal_] = projGrad(rateq,gradq,projq,Q0,nIter,1);
[x_armijo_,val_armijo_,value_armijo_,count_armijo_] = projGrad(rateq,gradq,projq,Q0,nIter,2);
%plot(0:count_optimal-1,value_optimal(1:count_optimal),'bo')

%plot(0:count_op_-1,value_op_(1:count_op_),'k-',0:count_-1,value_(1:count_),'c-',0:count_optimal_-1,value_optimal_(1:count_optimal_),'g-',0:count_armijo_-1,value_armijo_(1:count_armijo_),'y-',0:count_op-1,value_op(1:count_op),'b-',0:count-1,value(1:count),'r-',0:count_optimal-1,value_optimal(1:count_optimal),'g-',0:count_armijo-1,value_armijo(1:count_armijo),'m-');
%legend('covariance openloop','covariance fixed','covariance optimal','covariance armijo','precoder openloop','precoder fixed','precoder optimal','precoder armijo');


plot(0:count_op_-1,value_op_(1:count_op_),'k-',0:count_-1,value_(1:count_),'c-',0:count_armijo_-1,value_armijo_(1:count_armijo_),'y-',0:count_op-1,value_op(1:count_op),'b-',0:count-1,value(1:count),'r-',0:count_armijo-1,value_armijo(1:count_armijo),'m-');
legend('covariance openloop','covariance fixed','covariance armijo','precoder openloop','precoder fixed','precoder armijo');



%plot(0:count_op-1,value_op(1:count_op),'bo',0:count-1,value(1:count),'r+',0:count_optimal-1,value_optimal(1:count_optimal),'g*',0:count_armijo-1,value_armijo(1:count_armijo),'m--');
%0:count_optimal-1,value_optimal(1:count_optimal),'g',0:count_armijo-1,value_armijo(1:count_armijo),'b'
%plot(0:count_armijo-1,value_armijo(1:count_armijo),'bx')
end