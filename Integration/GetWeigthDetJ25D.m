function [maskp,inde]=GetWeigthDetJ25D(mesh_file,sizeim,pscale,GaussPts,selected_elts)

load(mesh_file,'rint','Xo','Yo','Zo','Nnodes','Nelems','Smesh');
if nargin < 3 , pscale=1;end
if nargin < 4 , GaussPts='pixels';end
if nargin < 5 , selected_elts=[];end


load(mesh_file,'conn','ng','ns','elt');

if strcmp(GaussPts,'pixels') ,ngt=ng*(1-0.75*(ng==4))*prod(ns)*2;ngq=ng*prod(ns);ngp=ng*prod(ns)*2;ngh=ng*prod(ns);st=1+(prod(ns)==1);
elseif strcmp(GaussPts,'Gauss_points'), ngt=1;ngq=1+3*(~rint);ngp=1+(~rint);ngh=1+7*(~rint);ng=1;ns=1;st=1;
elseif strcmp(GaussPts,'sub_cells'), ngt=1*prod(ns)*2;ngq=4*prod(ns);ngp=2*prod(ns)*2;ngh=1*prod(ns);ng=1;st=1+(prod(ns)==1);
end

if ng>0
    nx=ngt*sum(elt==3)/st+ngq*sum(elt==4)+ngp*sum(elt==6)/st+ngh*sum(elt==8);
else
    nx=prod(sizeim*pscale);
end

spformat=~isempty(selected_elts);
if spformat
    npts=ngt*sum(elt(selected_elts)==3)/st+ngq*sum(elt(selected_elts)==4)+ngp*sum(elt(selected_elts)==6)/st+ngh*sum(elt(selected_elts)==8);
    maskp=zeros(npts,1);
    indp=zeros(npts,1);
    inde=zeros(npts,1);
    nind=0;
else
    maskp=zeros(nx,1);
    inde=zeros(nx,1);
end
if ng>0
    if any(elt==3)
        [xgt,ygt,wgt]=GetGaussPointsTriangle(max(1,ngt/(2*prod(ns))),ns);
        Selt=length(xgt);
        Nt_r=[-1+0*xgt,1+0*xgt,0*ygt];
        Nt_s=[-1+0*ygt,0*xgt,1+0*ygt];
    end
    if any(elt==4)
        [xgq,ygq,wgq]=GetGaussPointsQuadrangle(ngq/(prod(ns)),ns);
        Selq=length(xgq);
        Nq_r=[-0.25*(1-ygq),0.25*(1-ygq),0.25*(1+ygq),-0.25*(1+ygq)];
        Nq_s=[-0.25*(1-xgq),-0.25*(1+xgq),0.25*(1+xgq),0.25*(1-xgq)];
    end
    if any(elt==6)
        [xgp,ygp,zgp,wgp]=GetGaussPointsWedge(ngp/(prod(ns)/2),ns);
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
                case 3
                    dxdr=Nt_r*Xo(inods);
                    dydr=Nt_r*Yo(inods);
                    dzdr=Nt_r*Zo(inods);
                    dxds=Nt_s*Xo(inods);
                    dyds=Nt_s*Yo(inods);
                    dzds=Nt_s*Zo(inods);
                    wg=wgt;
                    d3 = dxdr.*dyds - dydr.*dxds;
                    d2 = dzdr.*dxds - dxdr.*dzds;
                    d1 = dydr.*dzds - dzdr.*dyds;

                    detJ = sqrt ( d1.*d1 + d2.*d2 + d3.*d3 );

                    Sel=Selt;
                case 4
                    dxdr=Nq_r*Xo(inods);
                    dydr=Nq_r*Yo(inods);
                    dzdr=Nq_r*Zo(inods);
                    dxds=Nq_s*Xo(inods);
                    dyds=Nq_s*Yo(inods);
                    dzds=Nq_s*Zo(inods);
                    wg=wgq;
                    d3 = dxdr.*dyds - dydr.*dxds;
                    d2 = dzdr.*dxds - dxdr.*dzds;
                    d1 = dydr.*dzds - dzdr.*dyds;

                    detJ = sqrt ( d1.*d1 + d2.*d2 + d3.*d3 );

                    Sel=Selq;
                case 6
                    dxdr=Np_r*Xo(inods);
                    dydr=Np_r*Yo(inods);
                    dzdr=Np_r*Zo(inods);
                    dxds=Np_s*Xo(inods);
                    dyds=Np_s*Yo(inods);
                    dzds=Np_s*Zo(inods);
                    dxdt=Np_t*Xo(inods);
                    dydt=Np_t*Yo(inods);
                    dzdt=Np_t*Zo(inods);
                    detJ=dxdr.*dyds.*dzdt+dxds.*dydt.*dzdr+dxdt.*dydr.*dzds...
                        -dzdr.*dyds.*dxdt-dzds.*dydt.*dxdr-dzdt.*dydr.*dxds;
                    wg=wgp;
                    Sel=Selp;
                case 8
                    dxdr=Nh_r*Xo(inods);
                    dydr=Nh_r*Yo(inods);
                    dzdr=Nh_r*Zo(inods);
                    dxds=Nh_s*Xo(inods);
                    dyds=Nh_s*Yo(inods);
                    dzds=Nh_s*Zo(inods);
                    dxdt=Nh_t*Xo(inods);
                    dydt=Nh_t*Yo(inods);
                    dzdt=Nh_t*Zo(inods);
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

            case 3
                npix=npix+Selt;
            case 4
                npix=npix+Selq;
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
