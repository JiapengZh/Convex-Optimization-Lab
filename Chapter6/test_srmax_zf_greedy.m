function test_srmax_zf_greedy(filename)

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
            for p = 1:length(Parr)
                [R, q] = srmax_zf_greedy(testSets(sysSize,iter).H,Parr(p));
                
                if max(abs(R-testSets(sysSize,iter).srmax_zf_greedy_R(p))) > zeta
                    display(['FAIL, wrong result for sumrate, test No. ' num2str(iter) ',' num2str(p)])
                    keyboard
                    ok=false;
                end
                
                if max(abs(q-testSets(sysSize,iter).srmax_zf_greedy_q(:,p))) > zeta
                    display(['FAIL, wrong result for power allocation, test No. ' num2str(iter) ',' num2str(p)])
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
