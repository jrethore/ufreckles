function [dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis25D(mesh_file,sizeim,pscale,selected_nodes,GaussPts,full_size,selected_elts)
rint=false;
[pp,filname,ext]=fileparts(mesh_file);
if isempty(ext)
    mesh_file=[mesh_file,'.mat'];
end

load(mesh_file,'-mat','rint','Xo','Yo','Zo','Nnodes','Nelems','Smesh');

if nargin < 3 , pscale=1;end
if nargin < 4 , selected_nodes=[];end
if nargin < 5 , GaussPts='pixels';end
if nargin < 6 , full_size=false;end
if nargin < 7 , selected_elts=[];end

assert(~((pscale>1)&(strcmp(GaussPts,'pixels'))));

load(mesh_file,'-mat','conn','elt','ng','ns');
incp(1)=1;incp(2)=incp(1)*sizeim(1)*pscale;
if strcmp(GaussPts,'pixels') ,ngt=ng*(1-0.75*(ng==4))*prod(ns)*2;ngq=ng*prod(ns);ngp=ng*prod(ns)*2;ngh=ng*prod(ns)/2;st=1+(prod(ns)==1);
elseif strcmp(GaussPts,'Gauss_points'), ngt=1;ngq=1+3*(~rint);ngp=1+(~rint);ngh=1+7*(~rint);ng=1;ns=1;st=1;
elseif strcmp(GaussPts,'sub_cells'), ngt=1*prod(ns)*2;ngq=4*prod(ns);ngp=2*prod(ns)*2;ngh=8*prod(ns);ng=1;st=1+(prod(ns)==1);
end

assert(ng>0)
if ng>0
    npts=3*ngt*sum(elt==3)/st+4*ngq*sum(elt==4);nx=ngt*sum(elt==3)/st+ngq*sum(elt==4);
else
    npts=mean(elt)*ceil(prod(Smesh))*pscale^2;nx=prod(sizeim*pscale);
end
if isempty(selected_nodes)||full_size
    ny=prod(Nnodes);
else
    ny=length(selected_nodes);
end
indn=zeros(npts,1);
indp=zeros(npts,1);
valx=zeros(npts,1);
valy=zeros(npts,1);
valz=zeros(npts,1);
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
    
end
nel=0;
npix=0;
for i1=1:prod(Nelems)
    inods=conn(i1,1:elt(i1));
    go=0;
    if isempty(selected_nodes)
        go=1;
    else
        go=0;
        for in=1:elt(i1)
            if any(selected_nodes==conn(i1,in))
                go=1;
            end
        end
    end
    if go
        if elt(i1)==3
            dxdr=Nt_r*Xo(inods);
            dydr=Nt_r*Yo(inods);
            dzdr=Nt_r*Zo(inods);
            dxds=Nt_s*Xo(inods);
            dyds=Nt_s*Yo(inods);
            dzds=Nt_s*Zo(inods);
            
            d3 = dxdr.*dyds - dydr.*dxds;
            d2 = dzdr.*dxds - dxdr.*dzds;
            d1 = dydr.*dzds - dzdr.*dyds;
            
            DetJ = sqrt ( d1.*d1 + d2.*d2 + d3.*d3 );
            
            
            dxdt = d1./DetJ;
            dydt = d2./DetJ;
            dzdt = d3./DetJ;
            
            
            invJ(:,1) = (dyds .* dzdt - dydt .* dzds) ./ DetJ;
            invJ(:,4) = -(dydr .* dzdt - dydt .* dzdr) ./ DetJ;
            invJ(:,7)  = (dydr .* dzds - dyds .* dzdr) ./ DetJ;
            
            invJ(:,2)  = -(dxds .* dzdt - dxdt .* dzds) ./ DetJ;
            invJ(:,5)  = (dxdr .* dzdt - dxdt .* dzdr) ./ DetJ;
            invJ(:,8)  = -(dxdr .* dzds - dxds .* dzdr) ./ DetJ;
            
            invJ(:,3)  = (dxds .* dydt - dxdt .* dyds) ./ DetJ;
            invJ(:,6)  = -(dxdr .* dydt - dxdt .* dydr) ./ DetJ;
            invJ(:,9)  = (dxdr .* dyds - dxds .* dydr) ./ DetJ;
            for in=1:3
                if isempty(selected_nodes)||full_size
                    id=inods(in);
                else
                    id=find(selected_nodes==inods(in));
                end
                if ~isempty(id)
                    N_x=Nt_r(:,in).*invJ(:,1)+Nt_s(:,in).*invJ(:,4);
                    N_y=Nt_r(:,in).*invJ(:,2)+Nt_s(:,in).*invJ(:,5);
                    N_z=Nt_r(:,in).*invJ(:,3)+Nt_s(:,in).*invJ(:,6);
                    
                    indn(nel+(1:Selt))=id;
                    indp(nel+(1:Selt))=npix+(1:Selt);
                    valx(nel+(1:Selt))=N_x(:);
                    valy(nel+(1:Selt))=N_y(:);
                    valz(nel+(1:Selt))=N_z(:);
                    nel=nel+Selt;
                end
            end
            npix=npix+Selt;
        elseif elt(i1)==4
            dxdr=Nq_r*Xo(inods);
            dydr=Nq_r*Yo(inods);
            dzdr=Nq_r*Zo(inods);
            dxds=Nq_s*Xo(inods);
            dyds=Nq_s*Yo(inods);
            dzds=Nq_s*Zo(inods);
            
            d3 = dxdr.*dyds - dydr.*dxds;
            d2 = dzdr.*dxds - dxdr.*dzds;
            d1 = dydr.*dzds - dzdr.*dyds;
            
            DetJ = sqrt ( d1.*d1 + d2.*d2 + d3.*d3 );
            
            
            dxdt = d1./DetJ;
            dydt = d2./DetJ;
            dzdt = d3./DetJ;
            
            
            
            invJ(:,1) = (dyds .* dzdt - dydt .* dzds) ./ DetJ;
            invJ(:,4) = -(dydr .* dzdt - dydt .* dzdr) ./ DetJ;
            invJ(:,7)  = (dydr .* dzds - dyds .* dzdr) ./ DetJ;
            
            invJ(:,2)  = -(dxds .* dzdt - dxdt .* dzds) ./ DetJ;
            invJ(:,5)  = (dxdr .* dzdt - dxdt .* dzdr) ./ DetJ;
            invJ(:,8)  = -(dxdr .* dzds - dxds .* dzdr) ./ DetJ;
            
            invJ(:,3)  = (dxds .* dydt - dxdt .* dyds) ./ DetJ;
            invJ(:,6)  = -(dxdr .* dydt - dxdt .* dydr) ./ DetJ;
            invJ(:,9)  = (dxdr .* dyds - dxds .* dydr) ./ DetJ;
            for in=1:4
                if isempty(selected_nodes)||full_size
                    id=inods(in);
                else
                    id=find(selected_nodes==inods(in));
                end
                if ~isempty(id)
                    N_x=Nq_r(:,in).*invJ(:,1)+Nq_s(:,in).*invJ(:,4);
                    N_y=Nq_r(:,in).*invJ(:,2)+Nq_s(:,in).*invJ(:,5);
                    N_z=Nq_r(:,in).*invJ(:,3)+Nq_s(:,in).*invJ(:,6);
                    indn(nel+(1:Selq))=id;
                    indp(nel+(1:Selq))=npix+(1:Selq);
                    valx(nel+(1:Selq))=N_x(:);
                    valy(nel+(1:Selq))=N_y(:);
                    valz(nel+(1:Selq))=N_z(:);
                    nel=nel+Selq;
                end
            end
            npix=npix+Selq;
            
        end
    else
        if ng>0
            if elt(i1)==3
                npix=npix+Selt;
            elseif elt(i1)==4
                npix=npix+Selq;
                
            end
        end
    end
    clear invJ
end

if nel<length(indn)
    indn=indn(1:nel);
    indp=indp(1:nel);
    valx=valx(1:nel);
    valy=valy(1:nel);
    valz=valz(1:nel);
end

dphidx=sparse(indp,indn,valx,nx,ny);
dphidy=sparse(indp,indn,valy,nx,ny);
dphidz=sparse(indp,indn,valz,nx,ny);


end