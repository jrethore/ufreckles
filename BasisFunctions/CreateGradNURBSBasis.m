function [dphidx,dphidy,Xi,Yi,wi]=CreateGradNURBSBasis(mesh_file,p,GaussPts,pscale,frame)
if nargin < 3 , GaussPts='pixels';end
if nargin < 4 , pscale=1;end
if nargin < 5 , frame='parametric';end
load(mesh_file,'Px','Py','uo','vo','Nbselems');
uo=(uo-0.5)*pscale+0.5;
vo=(vo-0.5)*pscale+0.5;
Px=(Px-0.5)*pscale+0.5;
Py=(Py-0.5)*pscale+0.5;


switch GaussPts
    case 'pixels'
        [Yi,Xi]=meshgrid((vo(1+p(2))+0.5):(vo(length(vo)-p(2))-0.5),(uo(1+p(1))+0.5):(uo(length(uo)-p(1))-0.5));
        Xi=Xi(:);Yi=Yi(:);
        wi=1;
    case 'knot'
        [Yi,Xi]=meshgrid(vo,uo);
        Xi=Xi(:);Yi=Yi(:);
        wi=1;
    case 'nodes'
        load(mesh_file,'ui','vi');
        Yi=vi;
        Xi=ui;
        wi=1;
    case 'Gauss_points'
        [xgi,wxgi]=GetGaussPointsLine(1+p(1));
        [ygi,wygi]=GetGaussPointsLine(1+p(2));
        [yg,xg]=meshgrid(ygi,xgi);
        wg=wxgi*wygi';
        xg=xg(:);yg=yg(:);wg=wg(:);
        N=[0.25*(1-xg).*(1-yg),0.25*(1+xg).*(1-yg),0.25*(1+xg).*(1+yg),0.25*(1-xg).*(1+yg)];
        N_r=[-0.25*(1-yg),0.25*(1-yg),0.25*(1+yg),-0.25*(1+yg)];
        N_s=[-0.25*(1-xg),-0.25*(1+xg),0.25*(1+xg),0.25*(1-xg)];
        Xi=zeros(prod(Nbselems)*prod(1+p),1);
        Yi=zeros(prod(Nbselems)*prod(1+p),1);
        wi=zeros(prod(Nbselems)*prod(1+p),1);
        Sel=prod(1+p);
        np=0;
        for ix=1:Nbselems(1)
            for iy=1:Nbselems(2)
                inods=ix+p(1)+[0,1,1,0];
                jnods=iy+p(2)+[0,0,1,1];
                ui=uo(inods(:))';
                vi=vo(jnods(:))';
                dxdr=N_r*ui;
                dydr=N_r*vi;
                dxds=N_s*ui;
                dyds=N_s*vi;
                detJ=(dxdr.*dyds-dydr.*dxds);
                indp=np+(1:Sel);
                Xi(indp)=N*ui;
                Yi(indp)=N*vi;
                wi(indp)=detJ.*wg;
                np=np+Sel;
            end
        end
end
npt=numel(Xi);

indp=zeros((p(1)+1)*npt,1);
indn=zeros((p(1)+1)*npt,1);
val=zeros((p(1)+1)*npt,1);
dval=zeros((p(1)+1)*npt,1);

Nnx=Nbselems(1)+p(1);


nel=0;
for ix=1:Nbselems(1)
    if ix==Nbselems(1)
        found=find((Xi>=uo(ix+p(1)))&(Xi<=uo(ix+p(1)+1)));
    else
        found=find((Xi>=uo(ix+p(1)))&(Xi<uo(ix+p(1)+1)));
    end
    xp=Xi(found);
    [f]=NURBSBasisFunc(ix+p(1),p(1),xp',uo,1);
    
    Sel=length(xp);
    for ip=1:(p(1)+1)
        indn(nel+(1:Sel))=ix+ip-1;
        indp(nel+(1:Sel))=found;
        dval(nel+(1:Sel))=f(:,ip,2);
        val(nel+(1:Sel))=f(:,ip,1);
        nel=nel+Sel;
    end
end
phix=sparse(indp,indn,val,npt,Nnx);
dphix=sparse(indp,indn,dval,npt,Nnx);

indp=zeros((p(2)+1)*npt,1);
indn=zeros((p(2)+1)*npt,1);
val=zeros((p(2)+1)*npt,1);
dval=zeros((p(2)+1)*npt,1);

Nny=Nbselems(2)+p(2);
nel=0;
for iy=1:Nbselems(2)
    if iy==Nbselems(2)
        found=find((Yi>=vo(iy+p(2)))&(Yi<=vo(iy+p(2)+1)));
    else
        found=find((Yi>=vo(iy+p(2)))&(Yi<vo(iy+p(2)+1)));
    end
    yp=Yi(found);
    [f]=NURBSBasisFunc(iy+p(2),p(2),yp',vo,1);
    
    Sel=length(yp);
    for ip=1:(p(2)+1)
        indn(nel+(1:Sel))=iy+ip-1;
        indp(nel+(1:Sel))=found;
        dval(nel+(1:Sel))=f(:,ip,2);
        val(nel+(1:Sel))=f(:,ip,1);
        nel=nel+Sel;
    end
end
phiy=sparse(indp,indn,val,npt,Nny);
dphiy=sparse(indp,indn,dval,npt,Nny);


[indj,indi]=meshgrid(1:Nny,1:Nnx);
phi=phix(:,indi(:)).*phiy(:,indj(:));
dphidx=dphix(:,indi(:)).*phiy(:,indj(:));
dphidy=phix(:,indi(:)).*dphiy(:,indj(:));
if strcmp(frame,'physical')
    
                  dxdr=dphidx*Px(:);
                  dydr=dphidx*Py(:);
                  dxds=dphidy*Px(:);
                  dyds=dphidy*Py(:);
                  detJ=(dxdr.*dyds-dydr.*dxds);

    N_x=diag(sparse(dyds./detJ))*dphidx+diag(sparse(-dydr./detJ))*dphidy;
    N_y=diag(sparse(-dxds./detJ))*dphidx+diag(sparse(dxdr./detJ))*dphidy;
    dphidx=N_x;
    dphidy=N_y;
    wi=wi.*detJ;
end
Xi=phi*Px(:);
Yi=phi*Py(:);
wi=diag(sparse(wi));
end