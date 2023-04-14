function [maskp,inde]=GetWeigthDetJ(mesh_file,sizeim,pscale,GaussPts,selected_elts,deg)
zo=1;rint=0;
[pp,filname,ext]=fileparts(mesh_file);
if isempty(ext)
    mesh_file=[mesh_file,'.mat'];
end
load(mesh_file,'-mat','rint','xo','yo','zo','Nnodes','Nelems','Smesh');
if nargin < 3 , pscale=1;end
if nargin < 4 , GaussPts='pixels';end
if nargin < 5 , selected_elts=[];end
if nargin < 6 , deg=1;end
Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
if deg>1, rint=0;end

xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
if length(zo)>1
    zo=(zo-0.5)*pscale+0.5;
end
ns=[8,8];ng=1;

load(mesh_file,'-mat','conn','ng','ns','elt');
if length(zo)>1
    elt(elt==4)=5;
end
if strcmp(GaussPts,'pixels') ,ngb=ng*prod(ns);ngt=ng*(1-0.75*(ng==4))*prod(ns)*2;ngq=ng*prod(ns);ngte=ng*prod(ns)*2;ngp=ng*prod(ns)*2;ngh=ng*prod(ns);st=2;
elseif strcmp(GaussPts,'Gauss_points'), ngb=2;ngt=1+6*(deg==2);ngq=1+3*(~rint)+12*(deg==2);ngte=1;ngp=1+(~rint);ngh=1+7*(~rint);ng=1;ns=1;st=1;
elseif strcmp(GaussPts,'sub_cells'), ngb=2*prod(ns);ngt=(1+6*(deg==2))*prod(ns)*2;ngq=4*prod(ns);ngte=1*prod(ns)*2;ngp=2*prod(ns)*2;ngh=1*prod(ns);ng=1;st=2;
end

if ng>0
    nx=ngb*sum(elt==2)+ngt*sum(elt==3)+ngq*sum(elt==4)+ngte*sum(elt==5)+ngp*sum(elt==6)+ngh*sum(elt==8);
else
    nx=prod(sizeim*pscale);
end

spformat=~isempty(selected_elts);
if spformat
    npts=ngb*sum(elt(selected_elts)==2)+ngt*sum(elt(selected_elts)==3)+ngq*sum(elt(selected_elts)==4)+ngte*sum(elt(selected_elts)==5)+ngp*sum(elt(selected_elts)==6)+ngh*sum(elt(selected_elts)==8);
    maskp=zeros(npts,1);
    indp=zeros(npts,1);
    inde=zeros(npts,1);
    nind=0;
else
    maskp=zeros(nx,1);
    inde=zeros(nx,1);
end
if ng>0
    if any(elt==2)
        [xgb,wgb]=GetGaussPointsLine(ngb);
        Selb=length(xgb);
        Nb_r=[-0.5+0*xgb,0.5+0*xgb];
    end
    if any(elt==3)
        [xgt,ygt,wgt]=GetGaussPointsTriangle(ngt/(st*prod(ns)),ns);
        Selt=length(xgt);
        Nt_r=[-1+0*xgt,1+0*xgt,0*ygt];
        Nt_s=[-1+0*ygt,0*xgt,1+0*ygt];
    end
    if any(elt==5)
        [xgte,ygte,zgte,wgte]=GetGaussPointsTetrahedron(ngte,ns);
        Selte=length(xgte);
        Nte_r=[-1+0*xgte,1+0*xgte,0*ygte,0*zgte];
        Nte_s=[-1+0*xgte,0*xgte,1+0*ygte,0*zgte];
        Nte_t=[-1+0*xgte,0*xgte,0*ygte,1+0*zgte];
    end
    if any(elt==4)
        [xgq,ygq,wgq]=GetGaussPointsQuadrangle(ngq/(prod(ns)),ns);
        Selq=length(xgq);
        Nq_r=[-0.25*(1-ygq),0.25*(1-ygq),0.25*(1+ygq),-0.25*(1+ygq)];
        Nq_s=[-0.25*(1-xgq),-0.25*(1+xgq),0.25*(1+xgq),0.25*(1-xgq)];
    end
    if any(elt==6)
        [xgp,ygp,zgp,wgp]=GetGaussPointsWedge(ngp/(prod(ns)),ns);
        Selp=length(xgp);
        Np_r=[-0.5*(1-zgp),0.5*(1-zgp),(0*ygp),...
            -0.5*(1+zgp),0.5*(1+zgp),(0*ygp)];
        Np_s=[-0.5*(1-zgp),(0*xgp),0.5*(1-zgp),...
            -0.5*(1-zgp),(0*xgp),0.5*(1-zgp)];
        Np_t=[-0.5*(1-xgp-ygp),-0.5*xgp,-0.5*ygp,...
            0.5*(1-xgp-ygp),0.5*xgp,0.5*ygp];
    end
    if any(elt==8)
        [xgh,ygh,zgh,wgh]=GetGaussPointsHexaedron(ngh/prod(ns),ns);
        Selh=length(xgh);
        Nh_r=[-0.125*(1-ygh).*(1-zgh),0.125*(1-ygh).*(1-zgh),0.125*(1+ygh).*(1-zgh),-0.125*(1+ygh).*(1-zgh),...
            -0.125*(1-ygh).*(1+zgh),0.125*(1-ygh).*(1+zgh),0.125*(1+ygh).*(1+zgh),-0.125*(1+ygh).*(1+zgh)];
        Nh_s=[-0.125*(1-xgh).*(1-zgh),-0.125*(1+xgh).*(1-zgh),0.125*(1+xgh).*(1-zgh),0.125*(1-xgh).*(1-zgh),...
            -0.125*(1-xgh).*(1+zgh),-0.125*(1+xgh).*(1+zgh),0.125*(1+xgh).*(1+zgh),0.125*(1-xgh).*(1+zgh)];
        Nh_t=[-0.125*(1-xgh).*(1-ygh),-0.125*(1+xgh).*(1-ygh),-0.125*(1+xgh).*(1+ygh),-0.125*(1-xgh).*(1+ygh),...
            0.125*(1-xgh).*(1-ygh),0.125*(1+xgh).*(1-ygh),0.125*(1+xgh).*(1+ygh),0.125*(1-xgh).*(1+ygh)];
    end
    npix=0;
    for i1=1:prod(Nelems)
        inods=conn(i1,1:elt(i1));
        go=isempty(selected_elts);
        if ~go
            go=any(selected_elts==i1);
        end
        if go
            switch elt(i1)
                case 2
                    dxdr=Nb_r*xo(inods);
                    dydr=Nb_r*yo(inods);
                    wg=wgb;
                    detJ=sqrt(dxdr.^2+dydr.^2);
                    Sel=Selb;
                    
                case 3
                    dxdr=Nt_r*xo(inods);
                    dydr=Nt_r*yo(inods);
                    dxds=Nt_s*xo(inods);
                    dyds=Nt_s*yo(inods);
                    wg=wgt;
                    detJ=(dxdr.*dyds-dydr.*dxds);
                    Sel=Selt;
                case 4
                    dxdr=Nq_r*xo(inods);
                    dydr=Nq_r*yo(inods);
                    dxds=Nq_s*xo(inods);
                    dyds=Nq_s*yo(inods);
                    wg=wgq;
                    detJ=(dxdr.*dyds-dydr.*dxds);
                    Sel=Selq;
                case 5
                    inods=inods(1:4);
                    dxdr=Nte_r*xo(inods);
                    dydr=Nte_r*yo(inods);
                    dzdr=Nte_r*zo(inods);
                    dxds=Nte_s*xo(inods);
                    dyds=Nte_s*yo(inods);
                    dzds=Nte_s*zo(inods);
                    dxdt=Nte_t*xo(inods);
                    dydt=Nte_t*yo(inods);
                    dzdt=Nte_t*zo(inods);
                    detJ=dxdr.*dyds.*dzdt+dxds.*dydt.*dzdr+dxdt.*dydr.*dzds...
                        -dzdr.*dyds.*dxdt-dzds.*dydt.*dxdr-dzdt.*dydr.*dxds;
                    wg=wgte;
                    Sel=Selte;
                case 6
                    dxdr=Np_r*xo(inods);
                    dydr=Np_r*yo(inods);
                    dzdr=Np_r*zo(inods);
                    dxds=Np_s*xo(inods);
                    dyds=Np_s*yo(inods);
                    dzds=Np_s*zo(inods);
                    dxdt=Np_t*xo(inods);
                    dydt=Np_t*yo(inods);
                    dzdt=Np_t*zo(inods);
                    detJ=dxdr.*dyds.*dzdt+dxds.*dydt.*dzdr+dxdt.*dydr.*dzds...
                        -dzdr.*dyds.*dxdt-dzds.*dydt.*dxdr-dzdt.*dydr.*dxds;
                    wg=wgp;
                    Sel=Selp;
                case 8
                    dxdr=Nh_r*xo(inods);
                    dydr=Nh_r*yo(inods);
                    dzdr=Nh_r*zo(inods);
                    dxds=Nh_s*xo(inods);
                    dyds=Nh_s*yo(inods);
                    dzds=Nh_s*zo(inods);
                    dxdt=Nh_t*xo(inods);
                    dydt=Nh_t*yo(inods);
                    dzdt=Nh_t*zo(inods);
                    detJ=dxdr.*dyds.*dzdt+dxds.*dydt.*dzdr+dxdt.*dydr.*dzds...
                        -dzdr.*dyds.*dxdt-dzds.*dydt.*dxdr-dzdt.*dydr.*dxds;
                    wg=wgh;
                    Sel=Selh;
                    %                    npix=npix+Selh;
                    %             for ig=1:ngh
                    %                 J=[dxdr(ig),dxds(ig),dxdt(ig);...
                    %                     dydr(ig),dyds(ig),dydt(ig);...
                    %                     dzdr(ig),dzds(ig),dzdt(ig)];
                    %                 detJ=det(J);
                    %                 wg=wgh(ig);
                    %                 maskp(npix+1)=detJ*wg;
                    %                 npix=npix+1;
                    %                 if detJ<0
                    %                     error('DetJ <0');
                    %                 end
                    %             end
            end

            if spformat
                indp(nind+(1:Sel))=npix+(1:Sel);
                maskp(nind+(1:Sel))=detJ.*wg;
                inde(nind+(1:Sel))=i1;
                nind=nind+Sel;
            else
                maskp(npix+(1:Sel))=detJ.*wg;
                inde(npix+(1:Sel))=i1;
            end
            if any(detJ)<0
                error('DetJ <0');
            end
        end
        switch elt(i1)

            case 2
                npix=npix+Selb;
            case 3
                npix=npix+Selt;
            case 4
                npix=npix+Selq;
            case 5
                npix=npix+Selte;
            case 6
                npix=npix+Selp;
            case 8
                npix=npix+Selh;
        end
    end
    if spformat
        maskp=sparse(indp,1,maskp,nx,1);
        inde=sparse(indp,1,inde,nx,1);
    end
else
    maskp=1;
end


maskp=diag(sparse(maskp));



end
