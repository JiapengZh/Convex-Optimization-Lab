function test_c_mimo_mac(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-3;

ok = true;
try
    for t=1:length(testSetsWSR)
        R = c_mimo_mac(testSetsWSR(t).H1,testSetsWSR(t).H2,testSetsWSR(t).P1,testSetsWSR(t).P2,W);
        if max(abs(testSetsWSR(t).R(:)-R(:)))>zeta
            display(['FAIL, wrong rate points, test No. ' num2str(t)])
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