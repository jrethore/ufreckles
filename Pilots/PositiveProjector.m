function [Pxx,Pyy,Pxy]=PositiveProjector(exx,eyy,exy,P,opt)
if nargin<5,opt=1;end
switch opt
    case 1
    delta = 1e-6*max(abs(exx));
    Pxx=PositiveStrain(exx+delta,eyy,exy);
    Pxx=(Pxx-P)/delta;
    delta = 1e-6*max(abs(eyy));
    Pyy=PositiveStrain(exx,eyy+delta,exy);
    Pyy=(Pyy-P)/delta;
    delta = 1e-6*max(abs(exy));
    Pxy=PositiveStrain(exx,eyy,exy+delta);
    Pxy=(Pxy-P)/delta;
    case 2
    delta=sqrt((exx-eyy).^2+(exy).^2);
    m=(-exx+eyy+delta)./(exy);
    m(isnan(m))=0;
    m(isinf(m))=0;
    nn=(1)./sqrt(1+m.^2);
    lambdas=[exx+0.5*exy.*m,eyy-0.5*exy.*m];
    [~,ids]=sort(abs(lambdas),2);
    lambdas=[lambdas(:,1).*(ids(:,1)==1)+lambdas(:,2).*(ids(:,1)==2),lambdas(:,1).*(ids(:,2)==1)+lambdas(:,2).*(ids(:,2)==2)];
    nx=[nn,-nn.*m];
    nx=[nx(:,1).*(ids(:,1)==1)+nx(:,2).*(ids(:,1)==2),nx(:,1).*(ids(:,2)==1)+nx(:,2).*(ids(:,2)==2)];
    ny=[nn.*m,nn];
    ny=[ny(:,1).*(ids(:,1)==1)+ny(:,2).*(ids(:,1)==2),ny(:,1).*(ids(:,2)==1)+ny(:,2).*(ids(:,2)==2)];
    detV=1-2*((nx(:,1).*ny(:,2)-ny(:,1).*nx(:,2))>0);
    n1=[detV.*nx(:,1),detV.*ny(:,1)];
    n2=[nx(:,2),ny(:,2)];

    
    M1=[n1(:,1).*n1(:,1),n1(:,2).*n1(:,2),n1(:,1).*n1(:,2)];
    M2=[n2(:,1).*n2(:,1),n2(:,2).*n2(:,2),n2(:,1).*n2(:,2)];
    M12=[2*(n1(:,1).*n2(:,1)),2*(n1(:,2).*n2(:,2)),n1(:,1).*n2(:,2)+n2(:,1).*n1(:,2)];
    g1=(lambdas(:,1)+abs(lambdas(:,1)))/2;
    g2=(lambdas(:,2)+abs(lambdas(:,2)))/2;
    beta12=0.95<lambdas(:,1)./lambdas(:,2)<1.05;
    beta12=beta12.*((lambdas(:,1)>0)*0.5)+(1-beta12).*(0.5*(g1-g2)./(lambdas(:,1)-lambdas(:,2)));
    beta12(isnan(beta12))=0.5;
%    beta12=(0.5*(g1-g2)./(lambda1-lambda2));
    H1=ComputeKron(M1,M1);
    H2=ComputeKron(M2,M2);
    H12=ComputeKron(M12,M12);
    Pxx=diag(sparse(lambdas(:,1)>0))*H1+diag(sparse(lambdas(:,2)>0))*H2+diag(sparse(beta12))*H12;
    if any(isnan(Pxx(:)))
        keyboard
    end
    Pyy=full(Pxx(:,3+(1:3))*diag(sparse([1,1,1])));
    Pxy=full(Pxx(:,6+(1:3))*diag(sparse([2,2,2])));
    Pxx=full(Pxx(:,1:3)*diag(sparse([1,1,1])));
    case 3
        Pyy=zeros(size(exx,1),3);
        Pxy=zeros(size(exx,1),3);
        Pxx=zeros(size(exx,1),3);
        for ip=1:size(exx,1)
            [Pi]=TensorProjection2D(exx(ip),eyy(ip),exy(ip));
            Pxx(ip,:)=Pi(1,:);
            Pyy(ip,:)=Pi(2,:);
            Pxy(ip,:)=Pi(3,:);
        end
        
        
        
end
    function H=ComputeKron(M,N)
        H=[M(:,1).*N(:,1),M(:,1).*N(:,2),M(:,1).*N(:,3),...
            M(:,2).*N(:,1),M(:,2).*N(:,2),M(:,2).*N(:,3),...
            M(:,3).*N(:,1),M(:,3).*N(:,2),M(:,3).*N(:,3),];
    end
end