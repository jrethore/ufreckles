function [Nplus, Nminus]=ComputeNabla(lset,trupix)

siz=size(lset);

Dxplus=(lset(3:siz(1),2:siz(2)-1)-lset(2:siz(1)-1,2:siz(2)-1))/trupix;
Dxpp=Dxplus.*(repmat(1,size(Dxplus))+sign(Dxplus))/2;
Dxpm=Dxplus.*(repmat(1,size(Dxplus))-sign(Dxplus))/2;
Dxminus=(lset(2:siz(1)-1,2:siz(2)-1)-lset(1:siz(1)-2,2:siz(2)-1))/trupix;
Dxmp=Dxminus.*(repmat(1,size(Dxminus))+sign(Dxminus))/2;
Dxmm=Dxminus.*(repmat(1,size(Dxminus))-sign(Dxminus))/2;

Dyplus=(lset(2:siz(1)-1,3:siz(2))-lset(2:siz(1)-1,2:siz(2)-1))/trupix;
Dypp=Dyplus.*(repmat(1,size(Dyplus))+sign(Dyplus))/2;
Dypm=Dyplus.*(repmat(1,size(Dyplus))-sign(Dyplus))/2;
Dyminus=(lset(2:siz(1)-1,2:siz(2)-1)-lset(2:siz(1)-1,1:siz(2)-2))/trupix;
Dymp=Dyminus.*(repmat(1,size(Dyminus))+sign(Dyminus))/2;
Dymm=Dyminus.*(repmat(1,size(Dyminus))-sign(Dyminus))/2;

Nplus=repmat(1,siz);
Nminus=Nplus;

Nplus(2:siz(1)-1,2:siz(2)-1)=sqrt(max(Dxmp.*Dxmp,Dxpm.*Dxpm)+max(Dymp.*Dymp,Dypm.*Dypm));
Nminus(2:siz(1)-1,2:siz(2)-1)=sqrt(max(Dxmm.*Dxmm,Dxpp.*Dxpp)+max(Dymm.*Dymm,Dypp.*Dypp));
Nplus(1,:)=Nplus(2,:);
Nminus(1,:)=Nminus(2,:);
Nplus(siz(1),:)=Nplus(siz(1)-1,:);
Nminus(siz(1),:)=Nminus(siz(1)-1,:);
Nplus(:,1)=Nplus(:,2);
Nminus(:,1)=Nminus(:,2);
Nplus(:,siz(2))=Nplus(:,siz(2)-1);
Nminus(:,siz(2))=Nminus(:,siz(2)-1);
