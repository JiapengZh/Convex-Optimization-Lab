function test_sumrate_zf(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-3;

ok = true;
try
    for sysSize = 1:size(testSets,1)
        nIter = 20;
        if sysSize == 3
            nIter = 5;
        end
        for iter=1:nIter
            for i = 1:size(testSets(sysSize,iter).sel,2)
                for p = 1:length(Parr)
                    [R, q] = sumrate_zf(testSets(sysSize,iter).H,Parr(p),testSets(sysSize,iter).sel(:,i));
                    
                    if max(abs(R-testSets(sysSize,iter).sumrate_zf_R(i,p))) > zeta
                        display(['FAIL, wrong result for sumrate, test No. ' num2str(iter) ',' num2str(i) ',' num2str(p)])
                        keyboard
                        ok=false;
                    end
                    
                    if max(abs(q-testSets(sysSize,iter).sumrate_zf_q(:,i,p))) > zeta
                        display(['FAIL, wrong result for power allocation, test No. ' num2str(iter) ',' num2str(i) ',' num2str(p)])
                        keyboard
                        ok=false;
                    end
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
