function [Sloc,R]=GloToLoc(Sglo,Rold,DO)
if numel(Rold)>1
    if nargout>1
        assert(nargin==3);
        assert(numel(DO)>1);
        Ii=sparse(ones(size(DO,1),1));
        Oi=sparse(size(DO,1),1);
        I=[Ii,Oi,Oi,Oi,Ii,Oi,Oi,Oi,Ii];
        
%        full(DO(size(DO,1),:))
        
        DR=[Ii,-0.5*DO(:,1),-0.5*DO(:,2),...
            0.5*DO(:,1),Ii,-0.5*DO(:,3),...
            0.5*DO(:,2), 0.5*DO(:,3),Ii];
        
%        reshape(full(DR(size(DO,1),:)),3,3)
        
%        inv(reshape(full(DR(size(DO,1),:)),3,3))
        
        iDR=inverse(DR);
        
%        reshape(full(iDR(size(DO,1),:)),3,3)
        
        Q=mult(iDR,-2*(DR-I));
        
        Q=Q+I;
        R=mult(Q,Rold);
        
%        reshape(full(R(size(DO,1),:)),3,3)
        
    else
        R=Rold;
        if nargin<3,DO=0;end
        if DO
            R=R(:,[1,4,7,2,5,8,3,6,9]);
        end
    end
    if size(Sglo,2)<10
        if size(Sglo,2)==6,symflag=1;else symflag=0;end
        if symflag
            sym2full=sparse([1,4,6,4,2,5,6,5,3],1:9,1,6,9);
            Sglo=Sglo*sym2full;
        end
        RT=R(:,[1,4,7,2,5,8,3,6,9]);
        Stemp=mult(Sglo,R);
        Sloc=mult(RT,Stemp);
        if symflag
            full2sym=sparse([1,5,9,2,6,3],1:6,1,9,6);
            Sloc=Sloc*full2sym;
        end
    else
        indi1=repmat((1:3),6,2);
        indj1=repmat([1;2;3;1;2;1],6,1);
        indi2=repmat([1,2,3,2,3,1],6,1);
        indj2=repmat([1,2,3,2,3,3],1,6);
        R=Rold(:,indj1(:)+3*(indi1(:)-1)).*Rold(:,indj2(:)+3*(indi2(:)-1));
        indi1=repmat((1:3),3,2);
        indj1=repmat([1;2;1],6,1);
        indi2=repmat([1,2,3,2,3,1],3,1);
        indj2=[repmat([2;3;3],3,1);repmat([1;2;1],3,1)];
        indi=repmat(4:6,1,6);
        indj=repmat(1:6,3,1);

        complete=sparse(1:18,indi(:)+6*(indj(:)-1),1,18,36);
        R=R+(Rold(:,indj1(:)+3*(indi1(:)-1)).*Rold(:,indj2(:)+3*(indi2(:)-1)))*complete;
        trans=reshape(1:36,6,6);
        trans=trans';
        if DO
            R=R(:,trans(:));
        end
        Stemp=mult(Sglo,R);
        RT=R(:,trans(:));
        Sloc=mult(RT,Stemp);
    end
else
    Sloc=Sglo;
    R=Rold;
end
end
function B=inverse(A)
dete=A(:,1).*A(:,5).*A(:,9)+A(:,3).*A(:,4).*A(:,8)+A(:,2).*A(:,6).*A(:,7)...
    -A(:,3).*A(:,5).*A(:,7)-A(:,1).*A(:,6).*A(:,8)-A(:,2).*A(:,4).*A(:,9);
B=[(A(:,5).*A(:,9)-A(:,6).*A(:,8))./dete,...
    -(A(:,2).*A(:,9)-A(:,3).*A(:,8))./dete,...
    (A(:,2).*A(:,6)-A(:,3).*A(:,5))./dete,...
    -(A(:,4).*A(:,9)-A(:,6).*A(:,7))./dete,...
    (A(:,1).*A(:,9)-A(:,3).*A(:,7))./dete,...
    -(A(:,1).*A(:,6)-A(:,3).*A(:,4))./dete,...
    (A(:,4).*A(:,8)-A(:,5).*A(:,7))./dete,...
    -(A(:,1).*A(:,8)-A(:,2).*A(:,7))./dete,...
    (A(:,1).*A(:,5)-A(:,2).*A(:,4))./dete];

end
function C=mult(A,B)
nel=sqrt(size(A,2));
[indj,indi]=meshgrid(1:nel,1:nel);
C=sparse(1,1);
for k=1:nel
    C=C+A(:,k+(indi(:)-1)*nel).*B(:,indj(:)+(k-1)*nel);
end

end
