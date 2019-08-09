function test_srmax_brb(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-2;

ok = true;
try
    for sysSize = 1:size(testSets,1)
        nIter = 20;
        if sysSize == 3
            nIter = 5;
        end
        for iter=1:nIter
            for p = 1:length(Parr)
                [R, q] = srmax_brb(testSets(sysSize,iter).H,Parr(p),zeta);
                
                if max(abs(R-testSets(sysSize,iter).srmax_brb_R(p))) > 2*zeta
                    display(['FAIL, wrong result for sumrate, test No. ' num2str(sysSize) ',' num2str(iter) ',' num2str(p)])
                    keyboard
                    ok=false;
                end
                
                if max(abs(q-testSets(sysSize,iter).srmax_brb_q(:,p))) > 2*zeta
                    display(['FAIL, wrong result for power allocation, test No. ' num2str(sysSize) ',' num2str(iter) ',' num2str(p)])
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
