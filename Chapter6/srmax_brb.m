function [R q] = srmax_brb(H,P,epsilon)
K=size(H,3);
f = @(x,y) sumrate_SIMO(H,x,y);
red = @(a,b) reduction(P,a,b);
[q R] = brb(f,red,zeros(K,1),P*ones(K,1),epsilon);
end