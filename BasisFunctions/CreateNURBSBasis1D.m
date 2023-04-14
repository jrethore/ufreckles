function [phix,phiy]=CreateNURBSBasis1D(mesh_file,sizeim,p,pscale,timosh,type_nurbs)

if nargin < 4 , pscale=1;end
if nargin < 5 , timosh=0;end
if nargin < 6 , type_nurbs=1;end
load(mesh_file,'xo','yo','uo','vo','Nnodes','Nelems','Smesh');
xo=(xo-0.5)*pscale+1;
yo=(yo-0.5)*pscale+1;
xo=reshape(xo,Nnodes);
yo=reshape(yo,Nnodes);
xo=xo(1:Nnodes(1),1)';
yo=yo(1,1:Nnodes(2));
uo=(uo-1)*pscale+1;
vo=(vo-1)*pscale+1;


indp=zeros((p+1)*Smesh(2)*pscale,1);
indn=zeros((p+1)*Smesh(2)*pscale,1);
val=zeros((p+1)*Smesh(2)*pscale,1);
dval=zeros((p+1)*Smesh(2)*pscale,1);
if type_nurbs
Nny=length(yo)+p-1;
else
Nny=length(yo)*p-(p-1);    
end
nel=0;
toto=1;
for iy=1:Nelems(2)
    yp=(yo(iy):(yo(iy+1)-1));

	if type_nurbs
    [f,df]=NURBSDersBasisFunc(iy+p,p,yp',vo-0.5);
    else
    toto=toto+p;
    [f,df]=NURBSDersBasisFunc(toto,p,yp',vo-0.5);
    end

	Sel=length(yp);
    for ip=1:(p+1)
        if type_nurbs
        indn(nel+(1:Sel))=iy+ip-1;
        else
        indn(nel+(1:Sel))=iy*p-p+ip;
        end            
        indp(nel+(1:Sel))=yp;
		val(nel+(1:Sel))=f(:,ip);
		dval(nel+(1:Sel))=df(:,ip);
        nel=nel+Sel;
    end
end
phi1y=sparse(indp,indn,val,sizeim(2)*pscale,Nny);
dphi1y=sparse(indp,indn,dval,sizeim(2)*pscale,Nny);


%%
% figure1 = figure('XVisual',...
%     '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%     'PaperSize',[20.98 29.68]);
% axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
%     'FontName','Times');
% box('on');
% hold('all');
% 
% 
% xi=1:size(phi1y)-1;
% 
%     p1=plot(2*(xi-1)/xi(length(xi))-1,phi1y(xi,:),'LineWidth',2);
% 
% ylim([0,1]);
% 
% 
% %xlabel('Element size [pixel]','FontSize',20,'FontName','Times','Interpreter','latex','Interpreter','latex');
% %ylabel('Mean error [pixel]','FontSize',20,'FontName','Times','Interpreter','latex','Interpreter','latex');
% print ('-djpeg', ['p12-nurbs']);
% print ('-depsc', ['p12-nurbs']);
%     p1=plot(2*(vo-1)/xi(length(xi))-1,0*vo,'ko','LineWidth',2,'LineStyle','None','MarkerSize',10);
% print ('-djpeg', ['p12-nurbs-knot']);
% print ('-depsc', ['p12-nurbs-knot']);
% 2*(vo-1.)/xi(length(xi))-1
%%













% figure
% plot(phi1y)
% figure
% plot(dphi1y)
% [gy,gx]=gradient(phi1y);
% figure
% plot(gx)
phi1x=sparse(1:Smesh(1),1,1,sizeim(1)*pscale,1);
%dphi1x=sparse(1:Smesh(1),1,(-(0:Smesh(1)-1)/(Smesh(1)-1)+0.5)*2*Sel,sizeim(1)*pscale,1);
dphi1x=sparse(1:Smesh(1),1,-(0:Smesh(1)-1)+0.5*Smesh(1),sizeim(1)*pscale,1);

phix1=phi1x(:)*(phi1y(:)');
[indx,indy,val]=find(phix1);
indp=mod(indx,sizeim(1)*pscale)+sizeim(1)*pscale*(mod(indy,sizeim(2)*pscale)-1);
indn=ceil(indx/(sizeim(1)*pscale))+(ceil(indy/(sizeim(2)*pscale))-1);
phix1=sparse(indp,indn,val,prod(sizeim*pscale),Nny);

phiy1=dphi1x(:)*(dphi1y(:)');
[indx,indy,val]=find(phiy1);
indp=mod(indx,sizeim(1)*pscale)+sizeim(1)*pscale*(mod(indy,sizeim(2)*pscale)-1);
indn=ceil(indx/(sizeim(1)*pscale))+(ceil(indy/(sizeim(2)*pscale))-1);
phiy1=sparse(indp,indn,val,prod(sizeim*pscale),Nny);

% figure
% imagesc(reshape(phix1(:,1),sizeim))
% figure
% imagesc(reshape(phiy1(:,1),sizeim))

yp=1:Smesh(2);
f1=(yp-0.5)/(Smesh(2)+1);
f2=1-(yp-0.5)/(Smesh(2)+1);
f1=sparse(1:Smesh(2),1,f1',sizeim(2)*pscale,1);
f2=sparse(1:Smesh(2),1,f2',sizeim(2)*pscale,1);
phi2y=[f1,f2];
phiy2=phi1x(:)*(phi2y(:)');

[indx,indy,val]=find(phiy2);
indp=mod(indx,sizeim(1)*pscale)+sizeim(1)*pscale*(mod(indy,sizeim(2)*pscale)-1);
indn=ceil(indx/(sizeim(1)*pscale))+(ceil(indy/(sizeim(2)*pscale))-1);
phiy2=sparse(indp,indn,val,prod(sizeim*pscale),2);
phix2=sparse(size(phiy2,1),size(phiy2,2));

% figure
% imagesc(reshape(phix2(:,1),sizeim))
% figure
% imagesc(reshape(phiy2(:,1),sizeim))


% indp=zeros((1+1)*Smesh(2)*pscale,1);
% indn=zeros((1+1)*Smesh(2)*pscale,1);
% val=zeros((1+1)*Smesh(2)*pscale,1);
% Nny=length(yo)+1-1;
% nel=0;
% 					vo1=[yo(1),yo,yo(length(yo))];
% for iy=1:Nelems(2)
%     yp=(yo(iy):(yo(iy+1)-1));
% 
%     [f]=NURBSBasisFunc(iy+1,1,yp',vo1-0.5);
% 
% 	Sel=length(yp);
%     for ip=1:(1+1)
%         indn(nel+(1:Sel))=iy+ip-1;
%         indp(nel+(1:Sel))=yp;
% 		val(nel+(1:Sel))=f(:,ip);
%         nel=nel+Sel;
%     end
% end
% phi2y=sparse(indp,indn,val,sizeim(2)*pscale,Nny);
% phiy2=phi1x(:)*(phi2y(:)');
%  figure
%  plot(phi2y)
% 
% [indx,indy,val]=find(phiy2);
% indp=mod(indx,sizeim(1)*pscale)+sizeim(1)*pscale*(mod(indy,sizeim(2)*pscale)-1);
% indn=ceil(indx/(sizeim(1)*pscale))+(ceil(indy/(sizeim(2)*pscale))-1);
% phiy2=sparse(indp,indn,val,prod(sizeim*pscale),Nny);
% phix2=sparse(size(phiy2,1),size(phiy2,2));




if ~timosh
phix=[phix2,phix1];
phiy=[phiy2,phiy1];
else
phix=[phix2,phix1,0*phix1];
phiy=[phiy2,0*phiy1,phiy1];
    
end


end