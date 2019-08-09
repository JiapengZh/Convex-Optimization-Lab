function test_projQ(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-7;

ok = true;
try
    for t=1:length(testSets)
            Q = projQ(testSets(t).Q,testSets(t).P);
            if max(abs(testSets(t).projQ(:)-Q(:)))>zeta
                display(['FAIL, wrong projection, test No. ' num2str(t)])
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