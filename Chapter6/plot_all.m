Parr = [0.1 1 10 100 1000];
K = 5;
M = 2;

nIter = 100;

corr=eye(M);
R=[];

for iter = 1:nIter
    
    % channel and noise variance vectors
    H = permute(corr*sqrt(1/(2*M))*(randn(M,K)+1i*randn(M,K)),[1 3 2]);
    s = ones(K,1);
    
    
    for p = 1:length(Parr)
    q0 = ones(size(H,3),1)*Parr(p)/size(H,3);
        R(p,1,iter) = srmax_gradient(H,Parr(p),q0);
        R(p,2,iter) = srmax_zf_exhaustive(H,Parr(p));
        R(p,3,iter) = srmax_zf_greedy(H,Parr(p) );
        %R(p,4,iter) = srmax_brb(H,Parr(p) );
    end
    
    display(['iteration ' num2str(iter)])
    
end

figure;
plot(10*log10(Parr),mean(R,3));
legend('Your Gradient Solution','Your ZF Solution','Your Greedy ZF Solution','Global Optimum');
xlabel('P [dB]');
ylabel('R');
title(sprintf('%01d users, %01d transmit antennas, %01d channels realizations',K,M,nIter));