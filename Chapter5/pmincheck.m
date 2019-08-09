function ok = pmincheck(H,rho,sigma, P,T,status,P_opt,T_opt,status_opt)

% allowed errors
epsilon_P = 1e-4;
epsilon_r = 1e-5;

fprintf('\n');
disp(['Your optimization status is: ',yalmiperror(status)]);

if status ~= status_opt
    disp(['The status should be:',yalmiperror(status_opt)]);
    ok = false;
    return;
end
if status_opt > 0
    disp('Additional Question: is the problem infeasible or');
    disp('is the solver just unable to solve the problem?');
    ok = true;
    return;r_verify
end

gains = abs(H*T).^2;
r = log2(1+diag(gains)./(sigma(:)+sum(gains,2)-diag(gains)));
P_verify = (abs(P-P_opt)/P_opt <= epsilon_P);
r_verify = (max(abs(r-rho(:))) <= epsilon_r);
if P_verify&&r_verify,
    ok = true;
else
    disp('FAIL, the accuracy of your power minimization is insufficient!');
    disp(['Minimum Transmit-Power:  P_opt = ',num2str(P_opt),'; P = ',num2str(P),';']);
    disp(['Maximum distance to rate targets: max(r-rho) = ',num2str(max(abs(r-rho(:)))),';']);
    ok = false;
end
end