function plot_gradient()

load testData.mat

nIter = size(testSets,2);

gradSol = zeros(length(Parr),1);
optSol = zeros(length(Parr),1);

for sysSize = 1:size(testSets,1)-1
    in = 0;
    for p=1:length(Parr)
        gradSol(p)=0;
        optSol(p)=0;
        for iter=1:nIter
            in = in+1;
            H = testSets(sysSize,iter).H;
            P = Parr(p);
            
            q0 = ones(size(H,3),1)*P/size(H,3)+1;
            gradSol(p) = gradSol(p) + srmax_gradient(H,P,q0)/nIter;
            
            optSol(p) = optSol(p) + testSets(sysSize,iter).srmax_brb_R(p)/nIter;
        end
    end
figure
plot(10*log10(Parr),optSol,10*log10(Parr),gradSol);
legend('Global Optimum','Your Gradient Solution');
xlabel('P [dB]');
ylabel('R');

end