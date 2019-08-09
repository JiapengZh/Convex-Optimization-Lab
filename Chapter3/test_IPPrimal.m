noTests = 10;

zeta = 1e-2;
err = zeros(noTests,1);

ok = true;

for test = 1:noTests
    m = 200;
    n = 400;
    
    
    %optimal solutions
    x = rand(n,1);
    s = rand(n,1);
    I = rand(n,1)>0.5;
    x(I) = 0;
    s(~I) = 0;
    
    y = randn(m,1);
    
    %parameters
    A = randn(m,n);
    c = s + A'*y;
    b = A*x;
    
    [xh yh] = IPPrimal(c,A,b,0.5,1,1e-6);
    
    sh = max(0,c - A'*yh);
    
    rd = A'*yh + sh - c;
    rp = A*xh - b;
    rs = xh.*sh;
    
    err = max( [max(abs(rd)), max(abs(rp)), max(abs(rs))] );
    
    if err > zeta
        keyboard
        ok = false;
    end
    
end

if ok
    display('OK')
end
