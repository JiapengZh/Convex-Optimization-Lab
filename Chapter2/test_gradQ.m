function test_gradQ(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-7;

ok = true;
try
    for t=1:length(testSets)
            g = gradQ(testSets(t).Q,testSets(t).H);
            if max(abs(testSets(t).gradQ(:)-g(:)))>zeta
                display(['FAIL, wrong gradient, test No. ' num2str(t)])
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