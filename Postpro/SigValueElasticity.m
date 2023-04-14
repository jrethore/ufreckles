function X1=SigValueElasticity(omega)
% We try to solve the equation sin2aw+asin2w=0 and sin2aw-asin2w=0
f1=@(x) (sin(omega*x))^2-x^2*(sin(omega))^2;

Num = 2;
% The number of singularity values from solving two equations above


Time=10;

a = 0.01;
step = 0.01;
b = a+step;

for i=1:Num
    
    t=cputime;
    while (f1(a)*f1(b)>0)&&(cputime-t<Time)
        b=b+step;
    end;
    
    if (f1(a)*f1(b)<=0)        
        [X1(i),Er1(i),No1(i)]=bisection(f1,a,b);
        if No1(i)==0
            a=b+step;
            b=a+step;
        else
            a=b;
            b=a+step;
        end;
    else
        X1(i)=0; Er1(i)=NaN; No1(i)=0;
    end;
    
end;
end