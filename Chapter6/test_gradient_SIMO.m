function test_gradient_SIMO(filename)

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
                g=gradient_SIMO(testSets(sysSize,iter).H,testSets(sysSize,iter).x(:,i));
                if max(abs(g-testSets(sysSize,iter).gradient_SIMO(:,i))) > zeta
                    display(['FAIL, wrong result, test No. ' num2str(iter) ',' num2str(i)])
                    
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
