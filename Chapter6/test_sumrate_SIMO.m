function test_sumrate_SIMO(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-3;

ok = true;
try
    for sysSize = 1:size(testSets,1)
        nIter = 100;
        if sysSize == 3
            nIter = 5;
        end
        for iter=1:nIter
            for i = 1:size(testSets(sysSize,iter).x,2)
                R=sumrate_SIMO(testSets(sysSize,iter).H,testSets(sysSize,iter).x(:,i),testSets(sysSize,iter).y(:,i));
                if max(abs(R-testSets(sysSize,iter).sumrate_SIMO(i))) > zeta
                    display(['FAIL, wrong result, test No. ' num2str(t) ',' num2str(i)])
                    keyboard
                    ok=false;
                end
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
