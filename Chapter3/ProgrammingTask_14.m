%Programming Task 3.14

close all;

%Initialization
load 'IPNetworks.mat';

P0 = 1000;
N0 = 1;

dmax = 5;
res = zeros(9,10);
for net = 1:10
    
count = 1;
for alpha = 2:.5:6

%TODO

M = []; %network incident matrix
u = []; %vector that contains the capacity constraints

dis = IPWN(net).distance;
for i = 1:10
    for j = 1:10
        if dis(i,j) > dmax
            dis(i,j) = 0;
        elseif dis(i,j) == 0
        else
            gamma(i,j) = (P0/N0) * (dis(i,j)^(-alpha));
            U(i,j) = log2(1 + gamma(i,j));
        end
    end
end

[x_sort,xindex] = sort(IPWN(net).x,'ascend');
xmin = x_sort(1);
xmin_idx = xindex(1);
xmax = x_sort(10);
xmax_idx = xindex(10);
[y_sort,yindex] = sort(IPWN(net).y,'ascend');
ymin = y_sort(1);
ymin_idx = yindex(1);
ymax = y_sort(10);
ymax_idx = yindex(10);

t = 0;
for i = 1:10
    for j = 1:10
        if dis(i,j) ~= 0
            t = t + 1;
            M(i,t) = 1;
            M(j,t) = -1;
            u(t) = U(i,j);
        end
    end
end

a = size(M,2);
n = size(M,1);

M11 = M(1:min(xmin_idx,xmax_idx)-1,:);
M12 = M(min(xmin_idx,xmax_idx)+1:max(xmin_idx,xmax_idx)-1,:);
M13 = M(max(xmin_idx,xmax_idx)+1:end,:);
M1 = zeros(n-2,a);
M1(1:min(xmin_idx,xmax_idx)-1,:) = M11;
M1(min(xmin_idx,xmax_idx):max(xmin_idx,xmax_idx)-2,:) = M12;
M1(max(xmin_idx,xmax_idx)-1:end,:) = M13;
M2 = M1;


M31 = M(1:min(ymin_idx,ymax_idx)-1,:);
M32 = M(min(ymin_idx,ymax_idx)+1:max(ymin_idx,ymax_idx)-1,:);
M33 = M(max(ymin_idx,ymax_idx)+1:end,:);
M3 = zeros(n-2,a);
M3(1:min(ymin_idx,ymax_idx)-1,:) = M31;
M3(min(ymin_idx,ymax_idx):max(ymin_idx,ymax_idx)-2,:) = M32;
M3(max(ymin_idx,ymax_idx)-1:end,:) = M33;
M4 = M3;

A = [];
A(1:n-2,1:a) = M1;
A(n-1:2*n-4,a+1:2*a) = M2;
A(2*n-3:3*n-6,2*a+1:3*a) = M3;
A(3*n-5:4*n-8,3*a+1:4*a) = M4;
A(4*n-7:4*n-8+a,1:a) = eye(a);
A(4*n-7:4*n-8+a,a+1:2*a) = eye(a);
A(4*n-7:4*n-8+a,2*a+1:3*a) = eye(a);
A(4*n-7:4*n-8+a,3*a+1:4*a) = eye(a);
A(4*n-7:4*n-8+a,4*a+1:5*a) = eye(a);

b = zeros(4*n-8+a,1);
b(4*n-7:end) = u';

c = [];
e1 = zeros(10,1);
e2 = zeros(10,1);
e3 = zeros(10,1);
e4 = zeros(10,1);

e1(xmin_idx) = 1;
%e1(xmax_idx) = 1;
%e2(xmin_idx) = 1;
e2(xmax_idx) = 1;

%e3(ymin_idx) = 1;
e3(ymax_idx) = 1;
e4(ymin_idx) = 1;
%e4(ymax_idx) = 1;


c(1:a) = IPWN(net).weight(1) * e1' * M;
c(a+1:2*a) = IPWN(net).weight(2) * e2' * M;
c(2*a+1:3*a) = IPWN(net).weight(3) * e3' * M;
c(3*a+1:4*a) = IPWN(net).weight(4) * e4' * M;
c(4*a+1:5*a) = zeros(a,1);
c = c';


[x_cvx,y_cvx,s_cvx] = linprog_cvx(-c,A,b);

obj_cvx = c' * x_cvx;


res(count,net) = obj_cvx;
count = count+1;
end
end
plot(2:.5:6,res(1:9,1),'m-',2:.5:6,res(1:9,2),'b-',2:.5:6,res(1:9,3),'g-',2:.5:6,res(1:9,4),'c-',2:.5:6,res(1:9,5),'k-',2:.5:6,res(1:9,6),'r-',2:.5:6,res(1:9,7),'c-',2:.5:6,res(1:9,8),'y-',2:.5:6,res(1:9,9),'m--',2:.5:6,res(1:9,10),'c--');
legend('Network1','Network2','Network3','Network4','Network5','Network6','Network7','Network8','Network9','Network10');