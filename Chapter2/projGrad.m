function [x,val,value,count] = projGrad(fun,grad,proj,x,iterMax,method)

ksi = 10^(-10);
x = proj(x);
value = zeros(iterMax+1,1);
value(1) = fun(x);
count = 1;

if exist('method','var') == 0
      s = 1;
      for i = 1:iterMax
          xOld = x;
          x = xOld + s*grad(xOld);
          x = proj(x);
          value(i+1) = fun(x);
          val = fun(x);
          gApprox = (x - xOld)/s;
          frobNorm = gApprox(:)'*gApprox(:);
          count = count + 1;
          if frobNorm < ksi
            %x_star = x;
              val = fun(x);
        
              break
          end
      end
end
        
if exist('method', 'var')
    %error('nn')
    if method == 1
        s_start=100;
        for i = 1:iterMax
            %sa = 1 / i;
            xOld=x;
            fun_p=@(s) -fun(proj(x+s*grad(x)));
            s = fminbnd(fun_p,0,s_start);
            x = xOld + s*grad(xOld);
            x=proj(x);
            value(i+1) = fun(x);
            gApprox = (x - xOld)/s;
            frobNorm = gApprox(:)'*gApprox(:);
            val = fun(x);
            count = count + 1;
            s_start=2*s_start;
            if frobNorm < ksi
                %x_star = x;
                val = fun(x);
        
                break
            end
        end
        
    end
    if method == 2
        s_zero = 100;
        beta=0.5;
        sigma=0.2;
        fun_p=@(s) fun(proj(x+s*grad(x)));      
        for i = 1:iterMax
            s = s_zero;
            xOld = x;     
            %val = fun(x);
            pp = proj(x + s*grad(x)) - x;
            fanshushu = pp(:)'*pp(:);
            while (fun_p(s) - fun_p(0)) < (sigma*fanshushu)/s
                s = beta * s;
                %s_zero = 2 * s;
                pp = proj(x + s*grad(x)) - x;
                fanshushu = pp(:)'*pp(:);   
            end 
            
            x = xOld + s*grad(xOld);
            x = proj(x);
            value(i+1) = fun(x);
            val = fun(x);   
            count = count + 1;
            gApprox = (x - xOld)/s;
            frobNorm = gApprox(:)'*gApprox(:);
            if frobNorm < ksi
                %x_star = x;
                val = fun(x);      
                break
            end
            s_zero = 1*s;
            %fanshushu = norm(pp,'fro');
            
            
            %fanshushu = pp(:)'*pp(:);
            
            
            %if (fun_p(s) - fun_p(0)) >= (sigma*fanshushu)/s
            %    val = fun(x);
            %    break
            %end
            
        end
    end
    if method == 3
        for i = 1:iterMax
            sa = 1 / i;
            xOld = x;
            x = xOld + sa*grad(xOld);
            x = proj(x);
            value(i+1) = fun(x);
            val = fun(x);
            gApprox = (x - xOld)/sa;
            frobNorm = gApprox(:)'*gApprox(:);
            count = count + 1;
            if frobNorm < ksi
                %x_star = x;
                val = fun(x);
        
                break
            end
        end
        
    end
end
%disp('hanhan')
%disp('hanhan')
end