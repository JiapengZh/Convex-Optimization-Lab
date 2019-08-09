zeta = 1e-1;

load incidence.mat
terminal{1} = 4;
source{1} = 1;
ok = true;

terminal{2} = [4,3];
source{2} = [1,5];

for t = 1:length(terminal)
 

    if t == 1
        fprintf('\n\n*******Testing single flow*******\n\n')
    elseif t == 2
        fprintf('\n\n*******Testing multi flow*******\n\n')
    end
    
    [x,u,lambda,util,D,Util_pr,Util_pr_infeasible,iter] = networkflow_subgradient(M,terminal{t},source{t});
    [xOpt,uOpt,lambdaOpt,utilOpt] = networkflow_cvx(M,terminal{t},source{t});
    
    figure(t)
    plot(1:length(D), D(1:iter), 'm-',1:length(Util_pr), Util_pr(1:iter),'b-',1:length(Util_pr_infeasible), Util_pr_infeasible(1:iter),'c-');
    legend('Dual-Iter','Reconstructed-Iter','Reconstructedinf-Iter');
    
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