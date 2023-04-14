function [nxp,nyp,inde]=GetNormals(mesh_file,selected_elts)

load(mesh_file,'rint','xo','yo','zo','Nnodes','Nelems','Smesh');
pscale=1;
if nargin < 2 , selected_elts=[];end

xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
if length(zo)>1
    zo=(zo-0.5)*pscale+0.5;
end

load(mesh_file,'conn','ng','ns','elt');
if length(zo)>1
    elt(elt==4)=5;
end
ngb=2;ngt=1;ngq=1+3*(~rint);ngte=1;ngp=1+(~rint);ngh=1+7*(~rint);ng=1;ns=1;st=1;


    nx=ngb*sum(elt==2)+ngt*sum(elt==3)/st+ngq*sum(elt==4)+ngte*sum(elt==5)/st+ngp*sum(elt==6)/st+ngh*sum(elt==8);

spformat=~isempty(selected_elts);
if spformat
    npts=ngb*sum(elt(selected_elts)==2)+ngt*sum(elt(selected_elts)==3)/st+ngq*sum(elt(selected_elts)==4)+ngte*sum(elt(selected_elts)==5)/st++ngp*sum(elt(selected_elts)==6)/st+ngh*sum(elt(selected_elts)==8);
    nxp=zeros(npts,1);
    nyp=zeros(npts,1);
    indp=zeros(npts,1);
    inde=zeros(npts,1);
    nind=0;
else
    nxp=zeros(nx,1);
    nyp=zeros(nx,1);
    inde=zeros(nx,1);
end
    if any(elt==2)
        [xgb,wgb]=GetGaussPointsLine(ngb);
        Selb=length(xgb);
        Nb_r=[-0.5+0*xgb,0.5+0*xgb];
    end
    if any(elt==3)
        [xgt,ygt,wgt]=GetGaussPointsTriangle(max(1,ngt/(2*prod(ns))),ns);
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
                case 2
                    xn=xo(inods);
                    yn=yo(inods);
                    tx=diff(xn);
                    ty=diff(yn);
                         tnorm=abs(tx+1i*ty);
                        tx=(tx/tnorm);
                        ty=(ty/tnorm);
                        nye=-tx;
                        nxe=ty;
                   
                    Sel=Selb;
                                        
                otherwise
                    error('not coded yet');
            end

            if spformat
                indp(nind+(1:Sel))=npix+(1:Sel);
                nxp(nind+(1:Sel))=nxe;
                nyp(nind+(1:Sel))=nye;
                inde(nind+(1:Sel))=i1;
                nind=nind+Sel;
            else
                nxp(npix+(1:Sel))=nxe;
                nyp(npix+(1:Sel))=nye;
                inde(npix+(1:Sel))=i1;
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
        nxp=sparse(indp,1,nxp,nx,1);
        nyp=sparse(indp,1,nyp,nx,1);
        inde=sparse(indp,1,inde,nx,1);
    end

end
