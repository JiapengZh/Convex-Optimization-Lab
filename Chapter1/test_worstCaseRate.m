function test_worstCaseRate(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-2;

ok = true;
disp('Testing worstCaseRate')
try
    for t=1:length(testSetsWCR)
        for i=1:length(testSetsWCR(t).N)
            [R,s,mu] = worstCaseRate(testSetsWCR(t).p,testSetsWCR(t).N(i));
            
            if max(abs(testSetsWCR(t).R(i)-R))>zeta
                display(['FAIL, wrong rate, test No. ' num2str(t) ', with noise power ' num2str(testSetsWCR(t).N(i))])
                keyboard
                ok=false;
            end
            if max(abs(testSetsWCR(t).s(:,i)-s))>zeta
                display(['FAIL, wrong worst case powers, test No. ' num2str(t) ', with noise power ' num2str(testSetsWCR(t).N(i))])
                keyboard
                ok=false;
            end
            if max(abs(testSetsWCR(t).mu(i)-mu))>zeta
                display(['FAIL, wrong Lagrangian multiplier, test No. ' num2str(t) ', with noise power ' num2str(testSetsWCR(t).N(i))])
                keyboard
                ok=false;
            end
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

