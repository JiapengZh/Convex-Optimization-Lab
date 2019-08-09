function test_PMin(filename)

if nargin<1, filename = 'testData'; end

load(filename)

ok = true;
try
    for t=1:length(testSetsPMIN)
        [P,T,status] = PMinVectorBC_yalmip(testSetsPMIN(t).H,testSetsPMIN(t).rho,testSetsPMIN(t).sigma);
        %[P,T,status] = PMinVectorBC_cvx(testSetsPMIN(t).H,testSetsPMIN(t).rho,testSetsPMIN(t).sigma);
        ok = pmincheck(testSetsPMIN(t).H,testSetsPMIN(t).rho,testSetsPMIN(t).sigma,P,T,status, ...
            testSetsPMIN(t).P,testSetsPMIN(t).T,testSetsPMIN(t).status);
    end
catch e
    disp('FAIL, error');
    disp(getReport(e,'extended'))
    ok = false;
end
if ok
    disp('OK')
end