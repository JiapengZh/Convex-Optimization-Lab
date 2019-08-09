function test_reduction(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-3;

ok = true;
try
    for t=1:length(as)
        for p = 1:length(Parr)
            [a,b]=reduction(Parr(p),as{t},bs{t});
            
            if max(abs(a(:)-aSol{t,p}(:))) > zeta
                display(['FAIL, wrong result for lower bound, test No. ' num2str(t) ',' num2str(p)])
                keyboard
                ok=false;
            end
            if max(abs(b(:)-bSol{t,p}(:))) > zeta
                display(['FAIL, wrong result for upper bound, test No. ' num2str(t) ',' num2str(p)])
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
