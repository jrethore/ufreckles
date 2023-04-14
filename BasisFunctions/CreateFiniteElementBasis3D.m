function [phi]=CreateFiniteElementBasis3D(mesh_file,sizeim,pscale,selected_nodes,GaussPts,full_size,selected_elts)

load(mesh_file,'rint','xo','yo','zo','Nnodes','Nelems','Smesh');
if nargin < 3 , pscale=1;end
if nargin < 4 , selected_nodes=[];end
if nargin < 5 , GaussPts='pixels';end
if nargin < 6 , full_size=false;end
if nargin < 7 , selected_elts=[];end

xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
zo=(zo-0.5)*pscale+0.5;
load(mesh_file,'conn','elt','ng','ns');

incp(1)=1;incp(2)=incp(1)*sizeim(1)*pscale;incp(3)=incp(2)*sizeim(2)*pscale;
if strcmp(GaussPts,'rpixels')
    GaussPts='pixels';
    ng=0;
end
if strcmp(GaussPts,'pixels') ,ngte=ng*prod(ns)*2;ngt=ng*prod(ns)*2;ngq=ng*prod(ns);st=1+(prod(ns)==1);
elseif strcmp(GaussPts,'Gauss_points'),ngte=1;ngt=1+(~rint);ngq=1+7*(~rint);ng=1;ns=1;st=1;
elseif strcmp(GaussPts,'sub_cells'), ngte=1*prod(ns)*2;ngt=2*prod(ns)*2;ngq=1*prod(ns);ng=1;st=1+(prod(ns)==1);
end
if ng>0
    if isempty(selected_nodes)
        if isempty(selected_elts)
        npts=4*ngte*sum(elt==4)/st+6*ngt*sum(elt==6)/st+8*ngq*sum(elt==8);
        else
        npts=4*ngte*sum(elt(selected_elts)==4)/st+6*ngt*sum(elt(selected_elts)==6)/st+8*ngq*sum(elt(selected_elts)==8);
        end
    else
        if isempty(selected_elts)
        selected_elts=GetEltsFromNodes(conn,elt,selected_nodes);
        end
        npts=4*ngte*sum(elt(selected_elts)==4)/st+6*ngt*sum(elt(selected_elts)==6)/st+8*ngq*sum(elt(selected_elts)==8);       
    end
    nx=ngte*sum(elt==4)/st+ngt*sum(elt==6)/st+ngq*sum(elt==8);
else
    npts=ceil(mean(elt))*ceil(prod(Smesh))*pscale^3;nx=prod(sizeim*pscale);
end
if isempty(selected_nodes)||full_size
    ny=prod(Nnodes);
else
    ny=length(selected_nodes);
end
indn=zeros(npts,1);
indp=zeros(npts,1);
val=zeros(npts,1);
if ng>0
    if any(elt==4)
        [xgte,ygte,zgte,wgte]=GetGaussPointsTetrahedron(ngte,ns);
        Selte=length(xgte);
        Nte=[1-xgte-ygte-zgte,xgte,ygte,zgte];
    end
    if any(elt==6)
        [xgt,ygt,zgt,wgt]=GetGaussPointsWedge(ngt/(2*prod(ns)),ns);
        Selt=length(xgt);
        Nt=[0.5*(1-xgt-ygt).*(1-zgt),0.5*(xgt).*(1-zgt),0.5*(ygt).*(1-zgt),...
            0.5*(1-xgt-ygt).*(1+zgt),0.5*(xgt).*(1+zgt),0.5*(ygt).*(1+zgt)];
    end
    if any(elt==8)
        [xgq,ygq,zgq,wgq]=GetGaussPointsHexaedron(ngq/prod(ns),ns);
        Selq=length(xgq);
        Nq=[0.125*(1-xgq).*(1-ygq).*(1-zgq),0.125*(1+xgq).*(1-ygq).*(1-zgq),0.125*(1+xgq).*(1+ygq).*(1-zgq),0.125*(1-xgq).*(1+ygq).*(1-zgq),...
            0.125*(1-xgq).*(1-ygq).*(1+zgq),0.125*(1+xgq).*(1-ygq).*(1+zgq),0.125*(1+xgq).*(1+ygq).*(1+zgq),0.125*(1-xgq).*(1+ygq).*(1+zgq)];
    end

end
nel=0;
npix=0;
for i1=1:prod(Nelems)
    if isempty(selected_nodes)
        if isempty(selected_elts)
            go=1;
        else
            go=any(selected_elts==i1);
        end
    else
        if any(selected_elts==i1)
            go=1;
        else
            go=0;
        end
    end
    if go
    inods=conn(i1,1:elt(i1));
        if ng==0
            xn=xo(inods);
            yn=yo(inods);
            zn=zo(inods);
            [ypix,xpix,zpix]=meshgrid(ceil(min(yn)):floor(max(yn)),ceil(min(xn)):floor(max(xn)),ceil(min(zn)):floor(max(zn)));
            if elt(i1)==4
                if pscale==1
                [xg,yg,zg,wg]=GetGaussPointsTetrahedron(ngte,ns,xn,yn,zn,xpix(:),ypix(:),zpix(:));
                else
                    xg=(xpix-min(xn))/(max(xn)-min(xn));
                    yg=(ypix-min(yn))/(max(yn)-min(yn));
                    zg=(zpix-min(zn))/(max(zn)-min(zn));
                            wg=~(abs(xg)>1|abs(yg)>1|abs(zg)>1|(1-xg-yg-zg<0));

                end
                N=[1-xg-yg-zg,xg,yg,zg];
            elseif elt(i1)==6
                if pscale==1
                [xg,yg,zg,wg]=GetGaussPointsWedge(ngt,ns,xn,yn,zn,xpix(:),ypix(:),zpix(:));
                else
                    xg=-1+2*(xpix-min(xn))/(max(xn)-min(xn));
                    yg=-1+2*(ypix-min(yn))/(max(yn)-min(yn));
                    zg=-1+2*(zpix-min(zn))/(max(zn)-min(zn));
                            wg=~(abs(xg)>1|abs(yg)>1|abs(zg)>1);

                end
                N=[0.5*(1-xg-yg).*(1-zg),0.5*(xg).*(1-zg),0.5*(yg).*(1-zg),...
                    0.5*(1-xg-yg).*(1+zg),0.5*(xg).*(1+zg),0.5*(yg).*(1+zg)];
            elseif elt(i1)==8
                if pscale==1
                [xg,yg,zg,wg]=GetGaussPointsHexaedron(ngt,ns,xn,yn,zn,xpix(:),ypix(:),zpix(:));
                else
                     xg=-1+2*(xpix(:)-min(xn))/(max(xn)-min(xn));
                    yg=-1+2*(ypix(:)-min(yn))/(max(yn)-min(yn));
                    zg=-1+2*(zpix(:)-min(zn))/(max(zn)-min(zn));
                           wg=~(abs(xg)>1|abs(yg)>1|abs(zg)>1);
                end
                N=[0.125*(1-xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1+yg).*(1-zg),0.125*(1-xg).*(1+yg).*(1-zg),...
                    0.125*(1-xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1+yg).*(1+zg),0.125*(1-xg).*(1+yg).*(1+zg)];
            end
            Sel=length(xg);
            npix=1+incp(1)*(xpix-1)+incp(2)*(ypix-1)+incp(3)*(zpix-1);
            for in=1:elt(i1)
                if isempty(selected_nodes)||full_size
                    id=inods(in);
                else
                    id=find(selected_nodes==inods(in));
                end

                if ~isempty(id)
                    indn(nel+(1:Sel))=id;
                    indp(nel+(1:Sel))=npix;
                    val(nel+(1:Sel))=wg.*N(:,in);
                    nel=nel+Sel;
                end
            end
        else
            if elt(i1)==4
                for in=1:4
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)
                        indn(nel+(1:Selte))=id;
                        indp(nel+(1:Selte))=npix+(1:Selte);
                        val(nel+(1:Selte))=Nte(:,in);
                        nel=nel+Selte;
                    end
                end
            elseif elt(i1)==6
                for in=1:6
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)
                        indn(nel+(1:Selt))=id;
                        indp(nel+(1:Selt))=npix+(1:Selt);
                        val(nel+(1:Selt))=Nt(:,in);
                        nel=nel+Selt;
                    end
                end
%                npix=npix+Selt;
            elseif elt(i1)==8
                for in=1:8
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)

                        indn(nel+(1:Selq))=id;
                        indp(nel+(1:Selq))=npix+(1:Selq);
                        val(nel+(1:Selq))=Nq(:,in);
                        nel=nel+Selq;
                    end
                end
%                npix=npix+Selq;

            end
        end
    end
    if ng>0
        switch elt(i1)
            case 4
                npix=npix+Selte;
            case 6
                npix=npix+Selt;
            case 8
                npix=npix+Selq;
        end
    end
end


if nel<length(indn)
    indn=indn(1:nel);
    indp=indp(1:nel);
    val=val(1:nel);
end
phi=sparse(indp,indn,val,nx,ny);


end
