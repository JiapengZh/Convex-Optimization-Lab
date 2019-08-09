function test_wsr(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-2;

ok = true;
try
    for t=1:length(testSetsWSR)
        [Q1,Q2] = wsr_mimo_mac_yalmip(testSetsWSR(t).H1,testSetsWSR(t).H2,testSetsWSR(t).P1,testSetsWSR(t).P2,testSetsWSR(t).w);
        %[Q1,Q2] = wsr_mimo_mac_cvx(testSetsWSR(t).H1,testSetsWSR(t).H2,testSetsWSR(t).P1,testSetsWSR(t).P2,testSetsWSR(t).w);
        if max(abs(testSetsWSR(t).Q1(:)-Q1(:)))>zeta
            display(['FAIL, wrong Q1, test No. ' num2str(t)])
            keyboard
            ok=false;
        end
        if max(abs(testSetsWSR(t).Q2(:)-Q2(:)))>zeta
            display(['FAIL, wrong Q2, test No. ' num2str(t)])
            keyboard
            ok=false;
        end
    end
catch e
    disp('FAIL, error');
    disp(getReport(e,'extended'))
    ok = false;
end
if ok
    disp('OK')
end