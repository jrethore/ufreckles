
function [X,Er,No,count]=bisection(f,A,B)
% function [X,Er,No]=bisection(f,a,b)
% f : The function of being sovled
% x0 : The initial value of iteraction
% del: The error of x
% ep : The error of f(x)

% X : The returned value of solution
% Er : The value of function compared with 0
% No : The number of interaction after finishing

% f=2^(-x)+exp(x)+2*cos(x)-6;
 
% a=1;
% b=3;
 
del=10^(-8);   % danh gia gia tri x khi bn-an du gan nhau
ep=10^(-8); % danh gia gia tri f(x) di gan 0
% M=20;           % so lan lap cao nhat de ket thuc thuat toan
X=[]; Er=[]; No=[];
a=A;
b=B;
if (f(a)*f(b)>0)||(a>=b)
    
    disp('We can not use the bisection method');
    
elseif f(a)==0
   
    X=a;
    Er=ep;
    No=0;
    
elseif f(b)==0
    
    X=b;
    Er=ep;
    No=0;
    
else
    
    No = log((b-a)/ep)/log(2);

    
    mp=a+(b-a)/2;
    count=1;
    
    while abs(b-a) >= del && abs(f(mp))>=ep && count<=No
        
        if f(a)*f(b) <= 0
            if(f(a)*f(mp)<=0) 
                b=mp;
            else
                a=mp;
            end
        end
        
        mp=a+(b-a)/2;
        count=count+1;
    end
    
    
    mp=a+(b-a)/2;

    X=mp;
    Er=f(mp);
    No=count;
    
 %   toc
end;