function [phi,Xi,Yi,Zi,wgi]=CreateNURBSBasis3D(mesh_file,p,GaussPts,pscale)
if nargin < 3 , GaussPts='knot';end
if nargin < 4 , pscale=1;end
load(mesh_file,'Px','Py','Pz','uo','vo','wo','Nbselems');
uo=(uo-0.5)*pscale+0.5;
vo=(vo-0.5)*pscale+0.5;
wo=(wo-0.5)*pscale+0.5;
Px=(Px-0.5)*pscale+0.5;
Py=(Py-0.5)*pscale+0.5;
Pz=(Pz-0.5)*pscale+0.5;


switch GaussPts
    case 'pixels'
        error('Not allowed in 3D');
    case 'knot'
        [Yi,Xi,Zi]=meshgrid(vo,uo,wo);
        Xi=Xi(:);Yi=Yi(:);Zi=Zi(:);
        wgi=1;
    case 'nodes'
        load(mesh_file,'ui','vi','wi');
        Yi=vi;
        Xi=ui;
        Zi=wi;
        wgi=1;
    case 'Gauss_points'
        [xgi,wxgi]=GetGaussPointsLine(2*p(1));
        [ygi,wygi]=GetGaussPointsLine(2*p(2));
        [zgi,wzgi]=GetGaussPointsLine(2*p(3));
        [yg,xg,zg]=meshgrid(ygi,xgi,zgi);
        wg=wxgi*wygi';
        wg=wg(:)*wzgi';
        xg=xg(:);yg=yg(:);zg=zg(:);wg=wg(:);
        N=[0.125*(1-xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1+yg).*(1-zg),0.125*(1-xg).*(1+yg).*(1-zg),...
            0.125*(1-xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1+yg).*(1+zg),0.125*(1-xg).*(1+yg).*(1+zg)];
        N_r=[-0.125*(1-yg).*(1-zg),0.125*(1-yg).*(1-zg),0.125*(1+yg).*(1-zg),-0.125*(1+yg).*(1-zg),...
            -0.125*(1-yg).*(1+zg),0.125*(1-yg).*(1+zg),0.125*(1+yg).*(1+zg),-0.125*(1+yg).*(1+zg)];
        N_s=[-0.125*(1-xg).*(1-zg),-0.125*(1+xg).*(1-zg),0.125*(1+xg).*(1-zg),0.125*(1-xg).*(1-zg),...
            -0.125*(1-xg).*(1+zg),-0.125*(1+xg).*(1+zg),0.125*(1+xg).*(1+zg),0.125*(1-xg).*(1+zg)];
        N_t=[-0.125*(1-xg).*(1-yg),-0.125*(1+xg).*(1-yg),-0.125*(1+xg).*(1+yg),-0.125*(1-xg).*(1+yg),...
            0.125*(1-xg).*(1-yg),0.125*(1+xg).*(1-yg),0.125*(1+xg).*(1+yg),0.125*(1-xg).*(1+yg)];
        Xi=zeros(prod(Nbselems)*8*prod(p),1);
        Yi=zeros(prod(Nbselems)*8*prod(p),1);
        Zi=zeros(prod(Nbselems)*8*prod(p),1);
        wgi=zeros(prod(Nbselems)*8*prod(p),1);
        Sel=8*prod(p);
        np=0;
        for ix=1:Nbselems(1)
            for iy=1:Nbselems(2)
                for iz=1:Nbselems(3)
                    inods=ix+p(1)+[0,1,1,0,0,1,1,0];
                    jnods=iy+p(2)+[0,0,1,1,0,0,1,1];
                    knods=iz+p(3)+[0,0,0,0,1,1,1,1];
                    ui=uo(inods(:))';
                    vi=vo(jnods(:))';
                    wi=wo(knods(:))';
                    dxdr=N_r*ui;
                    dydr=N_r*vi;
                    dzdr=N_r*wi;
                    dxds=N_s*ui;
                    dyds=N_s*vi;
                    dzds=N_s*wi;
                    dxdt=N_t*ui;
                    dydt=N_t*vi;
                    dzdt=N_t*wi;
                    detJ =dxdr.*dyds.*dzdt + dxdt .*dydr.*dzds +...
                        dxds.*dydt.*dzdr - dxdt .*dyds.*dzdr -...
                        dxdr.*dydt.*dzds - dxds .*dydr.*dzdt;
                    indp=np+(1:Sel);
                    Xi(indp)=N*ui;
                    Yi(indp)=N*vi;
                    Zi(indp)=N*wi;
                    wgi(indp)=detJ.*wg;
                    np=np+Sel;
                end
            end
        end
end
npt=numel(Xi);
indp=zeros((p(1)+1)*npt,1);
indn=zeros((p(1)+1)*npt,1);
val=zeros((p(1)+1)*npt,1);

Nnx=Nbselems(1)+p(1);


nel=0;
for ix=1:Nbselems(1)
    if ix==Nbselems(1)
        found=find((Xi>=uo(ix+p(1)))&(Xi<=uo(ix+p(1)+1)));
    else
        found=find((Xi>=uo(ix+p(1)))&(Xi<uo(ix+p(1)+1)));
    end
    xp=Xi(found);
    [f]=NURBSBasisFunc(ix+p(1),p(1),xp',uo);
    
    Sel=length(xp);
    for ip=1:(p(1)+1)
        indn(nel+(1:Sel))=ix+ip-1;
        indp(nel+(1:Sel))=found;
        val(nel+(1:Sel))=f(:,ip);
        nel=nel+Sel;
    end
end
phix=sparse(indp,indn,val,npt,Nnx);

indp=zeros((p(2)+1)*npt,1);
indn=zeros((p(2)+1)*npt,1);
val=zeros((p(2)+1)*npt,1);

Nny=Nbselems(2)+p(2);
nel=0;
for iy=1:Nbselems(2)
    if iy==Nbselems(2)
        found=find((Yi>=vo(iy+p(2)))&(Yi<=vo(iy+p(2)+1)));
    else
        found=find((Yi>=vo(iy+p(2)))&(Yi<vo(iy+p(2)+1)));
    end
    yp=Yi(found);
    [f]=NURBSBasisFunc(iy+p(2),p(2),yp',vo);
    
    Sel=length(yp);
    for ip=1:(p(2)+1)
        indn(nel+(1:Sel))=iy+ip-1;
        indp(nel+(1:Sel))=found;
        val(nel+(1:Sel))=f(:,ip);
        nel=nel+Sel;
    end
end
phiy=sparse(indp,indn,val,npt,Nny);

indp=zeros((p(3)+1)*npt,1);
indn=zeros((p(3)+1)*npt,1);
val=zeros((p(3)+1)*npt,1);

Nnz=Nbselems(3)+p(3);
nel=0;
for iz=1:Nbselems(3)
    if iz==Nbselems(3)
        found=find((Zi>=wo(iz+p(3)))&(Zi<=wo(iz+p(3)+1)));
    else
        found=find((Zi>=wo(iz+p(3)))&(Zi<wo(iz+p(3)+1)));
    end
    zp=Zi(found);
    [f]=NURBSBasisFunc(iz+p(3),p(3),zp',wo);
    
    Sel=length(zp);
    for ip=1:(p(3)+1)
        indn(nel+(1:Sel))=iz+ip-1;
        indp(nel+(1:Sel))=found;
        val(nel+(1:Sel))=f(:,ip);
        nel=nel+Sel;
    end
end
phiz=sparse(indp,indn,val,npt,Nnz);

[indj,indi,indk]=meshgrid(1:Nny,1:Nnx,1:Nnz);
phi=phix(:,indi(:)).*phiy(:,indj(:)).*phiz(:,indk(:));

Xi=phi*Px(:);
Yi=phi*Py(:);
Zi=phi*Pz(:);

wgi=diag(sparse(wgi));
end