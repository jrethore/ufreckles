function [N]=NURBSBasisFunc(ii,pl,u,u_knotl,nders)
if nargin < 5 , nders=0;end


N=zeros([length(u),pl+1,nders+1]);
ndu=zeros(length(u),pl+1,pl+1);
a=zeros(length(u),pl+1,nders+1);



leftl=zeros(length(u),pl+1);
rightl=zeros(length(u),pl+1);

N(:,1,1) = 1;
ndu(:,1,1)=1;

for j = 1:pl
    leftl(:,j+1) = u - u_knotl(ii+1-j);
    rightl(:,j+1) = u_knotl(ii+j) - u;
    saved = 0;
    for r = 0:(j-1)
        ndu(:,j+1,r+1) = rightl(:,r+2)+leftl(:,j-r+1);
        temp = ndu(:,r+1,j)./ndu(:,j+1,r+1);
        ndu(:,r+1,j+1) = saved + rightl(:,r+2).*temp;
        saved = leftl(:,j-r+1).*temp;
    end
    ndu(:,j+1,j+1) = saved;
end

    N(:,:,1)=ndu(:,:,pl+1);
    
for r=0:pl
    s1=0;
    s2=1;
    a(:,1,1)=1;
    for k=1:nders
        d=zeros(length(u),1,1);
        rk=r-k;
        pk=pl-k;
        if(r>=k)
            a(:,s2+1,1)=a(:,s1+1,1)./ndu(:,pk+2,rk+1);
            d(:,1,1) = a(:,s2+1,1).*ndu(:,rk+1,pk+1);
        end
        if(rk>=-1)
            j1=1;
        else
            j1=-rk;
        end
        if((r-1)<=pk)
            j2=k-1;
        else
            j2=pl-r;
        end
        for j=j1:j2
            a(:,s2+1,j+1)=(a(:,s1+1,j+1) - a(:,s1+1,j))./ndu(:,pk+2,rk+j+1);
            d(:,1,1) = d(:,1,1) + a(:,s2+1,j+1).*ndu(:,rk+j+1,pk+1);
        end
        if(r<=pk)
            a(:,s2+1,k+1) = -a(:,s1+1,k)./ndu(:,pk+2,r+1);
            d(:,1,1) = d(:,1,1) + a(:,s2+1,k+1).*ndu(:,r+1,pk+1);
        end
        N(:,r+1,k+1)=d(:,1,1);
        j=s1;
        s1=s2;
        s2=j;
    end
end

% multiply by correction factor
r=pl;
for k=1:nders
    for j=0:pl
        N(:,j+1,k+1)=N(:,j+1,k+1)*r;
    end
    r=r*(pl-k);
end


end


%       N=zeros(length(u),pl+1);
%       leftl=zeros(length(u),pl+1);
%       rightl=zeros(length(u),pl+1);
%       N(:,1) = 1;
%       for j = 1:pl
%          leftl(:,j+1) = u - u_knotl(ii+1-j);
%          rightl(:,j+1) = u_knotl(ii+j) - u;
%          saved = 0;
%          for r = 0:(j-1)
%             temp = N(:,r+1)./(rightl(:,r+2)+leftl(:,j-r+1));
%             N(:,r+1) = saved + rightl(:,r+2).*temp;
%             saved = leftl(:,j-r+1).*temp;
%          end
%          N(:,j+1) = saved;
%       end
