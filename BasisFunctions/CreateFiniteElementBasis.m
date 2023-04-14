function [phi]=CreateFiniteElementBasis(mesh_file,sizeim,pscale,selected_nodes,GaussPts,full_size,deg)
[pp,filname,ext]=fileparts(mesh_file);
if isempty(ext)
    mesh_file=[mesh_file,'.mat'];
end
rflag=1;
load(mesh_file,'-mat','rflag','rint','xo','yo','Nnodes','Nelems');
if nargin < 3 , pscale=1;end
if nargin < 4 , selected_nodes=[];end
if nargin < 5 , GaussPts='pixels';end
if nargin < 6 , full_size=false;end
if nargin < 7 , deg=1;end
if deg>1, rint=0;end
Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
        
xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
ns=[8,8];ng=1;
load(mesh_file,'-mat','conn','elt','ng','ns');

incp(1)=1;incp(2)=incp(1)*sizeim(1)*pscale;
if strcmp(GaussPts,'rpixels')
    GaussPts='pixels';
    ng=0;
end
if strcmp(GaussPts,'pixels') ,ngb=ng*prod(ns);ngt=ng*(1-0.75*(ng==4))*prod(ns)*2;ngq=ng*prod(ns);ngp=ng*prod(ns)*2;ngh=ng*prod(ns)/2;st=2;
elseif strcmp(GaussPts,'Gauss_points'),ngb=2;  ngt=1+6*(deg==2);ngq=1+3*(~rint)+12*(deg==2);ngp=1+(~rint);ngh=1+7*(~rint);ng=1;ns=1;st=1;
elseif strcmp(GaussPts,'sub_cells'),ngb=2*prod(ns);  ngt=(1+6*(deg==2))*prod(ns)*2;ngq=4*prod(ns);ngp=2*prod(ns)*2;ngh=8*prod(ns);ng=1;st=2;
end
if ng>0
    npts=2*ngb*sum(elt==2)+3*ngt*sum(elt==3)+4*ngq*sum(elt==4);
    nx=ngb*sum(elt==2)+ngt*sum(elt==3)+ngq*sum(elt==4);
else
    npts=ceil(mean(elt))*ceil(prod(Smesh))*pscale^2;
    nx=prod(sizeim*pscale);
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
    if any(elt==2)
        [xgb,wgb]=GetGaussPointsLine(ngb);
        Selb=length(xgb);
        Nb=[0.5*(1-xgb),0.5*(1+xgb)];
    end
    if any(elt==3)
        [xgt,ygt,wgt]=GetGaussPointsTriangle(ngt/(st*prod(ns)),ns);
        Selt=length(xgt);
        Nt=[1-xgt-ygt,xgt,ygt];
    end
    if any(elt==4)
        [xgq,ygq,wgq]=GetGaussPointsQuadrangle(ngq/(prod(ns)),ns);
        Selq=length(xgq);
        Nq=[0.25*(1-xgq).*(1-ygq),0.25*(1+xgq).*(1-ygq),0.25*(1+xgq).*(1+ygq),0.25*(1-xgq).*(1+ygq)];
    end
end
nel=0;
npix=0;

for i1=1:prod(Nelems)
    inods=conn(i1,1:elt(i1));
    go=0;
    if isempty(selected_nodes)
        go=1;
    else
        for in=1:elt(i1)
            if any(selected_nodes==conn(i1,in))
                go=1;
            end
        end
    end
    if go
        if ng==0
            xn=xo(inods);
            yn=yo(inods);

                ipix=max(1,ceil(min(xn)-2)):min(sizeim(1),floor(max(xn)+2));
                jpix=max(1,ceil(min(yn)-2)):min(sizeim(2),floor(max(yn)+2));
                [ypix,xpix]=meshgrid(jpix,ipix);
                if rflag&&elt(i1)==4
                      xg=-1+2*(xpix(:)-min(xn))/(max(xn)-min(xn));
                    yg=-1+2*(ypix(:)-min(yn))/(max(yn)-min(yn));
                            wg=~(abs(xg)>1|abs(yg)>1);
                else
                     [xg,yg,wg]=GetGaussPointsPixels(elt(i1),xn,yn,xpix(:),ypix(:));
                   
                end
        
            
%                [ypix,xpix]=meshgrid(ceil(min(yn)):floor(max(yn)),ceil(min(xn)):floor(max(xn)));
%            [xg,yg,wg]=GetGaussPointsPixels(elt(i1),xn,yn,xpix(:),ypix(:));
            

            
                ind=find(wg);
                xg=xg(ind);yg=yg(ind);
            assert(elt(i1)>2,'Not Coded yet')
                N=GetFiniteElementShapeFunctions(elt(i1),xg,yg);
%             if elt(i1)==3
%                 figure
%             plot(xn,yn,'bx')
%             hold on
%             plot(xpix,ypix,'k+')
%               plot(N(:,1:elt(i1))*xn,N(:,1:elt(i1))*yn,'ro')
%             pause
%             end
%             if elt(i1)==3
%                 [xg,yg,wg]=GetGaussPointsTriangle(ngt/(prod(ns)/2),ns,xn,yn,xpix(:),ypix(:));
%                 ind=find(wg);
%                 xg=xg(ind);yg=yg(ind);
%                 N=[1-xg-yg,xg,yg];
%             elseif elt(i1)==4
%                 [xg,yg,wg]=GetGaussPointsQuadrangle(ngq/(prod(ns)),ns,xn,yn,xpix(:),ypix(:));
%                 ind=find(wg);
%                 xg=xg(ind);yg=yg(ind);
%                 N=[0.25*(1-xg).*(1-yg),0.25*(1+xg).*(1-yg),0.25*(1+xg).*(1+yg),0.25*(1-xg).*(1+yg)];
%             end
            Sel=length(xg);
            npix=1+incp(1)*(xpix(ind)-1)+incp(2)*(ypix(ind)-1);
            for in=1:elt(i1)
                if isempty(selected_nodes)||full_size
                    id=inods(in);
                else
                    id=find(selected_nodes==inods(in));
                end
                
                if ~isempty(id)
                    indn(nel+(1:Sel))=id;
                    indp(nel+(1:Sel))=npix;
                    val(nel+(1:Sel))=N(:,in);
                    nel=nel+Sel;
                end
            end
        else
            if elt(i1)==2
                for in=1:2
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)
                        indn(nel+(1:Selb))=id;
                        indp(nel+(1:Selb))=npix+(1:Selb);
                        val(nel+(1:Selb))=Nb(:,in);
                        nel=nel+Selb;
                    end
                end
                npix=npix+Selb;
            elseif elt(i1)==3
                for in=1:3
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
                npix=npix+Selt;
            elseif elt(i1)==4
                for in=1:4
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
                npix=npix+Selq;
                
            end
        end
    else
        if ng>0
            if elt(i1)==2
                npix=npix+Selb;
            elseif elt(i1)==3
                npix=npix+Selt;
            elseif elt(i1)==4
                npix=npix+Selq;
            end
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
