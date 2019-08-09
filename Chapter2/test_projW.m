function test_projW(filename)

if nargin<1, filename = 'testData'; end

load(filename)

zeta = 1e-7;

ok = true;
try
    for t=1:length(testSets)
            W = projW(testSets(t).W,testSets(t).P);
            if max(abs(testSets(t).projW(:)-W(:)))>zeta
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