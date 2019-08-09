zeta = 1e-1;

load incidence.mat
terminal{1} = 4;
source{1} = 1;


terminal{2} = [4,3];
source{2} = [1,5];
ok = true;

for t = 1:length(terminal)

    
    if t == 1
        fprintf('\n\n*******Testing single flow*******\n\n')
    elseif t == 2
        fprintf('\n\n*******Testing multi flow*******\n\n')
    end
    
    [x,u,lambda,util,D,LB,UB,Util_pr,Util_cu,iter] = networkflow_cuttingplane(M,terminal{t},source{t});
    [xOpt,uOpt,lambdaOpt,utilOpt] = networkflow_cvx(M,terminal{t},source{t});
    
    figure(t)
    plot(1:length(D), D(1:iter), 'm-',1:length(LB),LB(1:iter),'c-',1:length(UB),UB(1:iter),'g-',1:length(Util_pr), Util_pr(1:iter),'b-',1:length(Util_cu), Util_cu(1:iter),'y-');
    legend('Dual-Iter','Lower Bound','Upper Bound','Reconstructed-Iter','Reconstructed_cu-Iter');
    
    if max(abs(xOpt-x))>zeta
        display(['FAIL, wrong flows, test No. ' num2str(t)])
        keyboard
        ok=false;
    end
    if max(abs(uOpt-u))>zeta
        display(['FAIL, wrong link capacities, test No. ' num2str(t)])
        keyboard
        ok=false;
    end
    if max(abs(lambdaOpt-lambda))>zeta
        display(['WARNING, wrong lambda, test No. ' num2str(t)])
        display(['this is usually not a problem'])
    end
    if max(abs(utilOpt-util))>zeta
        display(['FAIL, wrong utility, test No. ' num2str(t)])
        keyboard
        ok=false;
    end
end

if ok
    disp('OK')
end