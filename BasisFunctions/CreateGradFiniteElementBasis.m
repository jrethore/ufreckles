function [dphidx,dphidy]=CreateGradFiniteElementBasis(mesh_file,sizeim,pscale,selected_nodes,GaussPts,full_size,deg)
rint=false;
[pp,filname,ext]=fileparts(mesh_file);
if isempty(ext)
    mesh_file=[mesh_file,'.mat'];
end
load(mesh_file,'-mat','rint','xo','yo','Nnodes','Nelems','Smesh');

if nargin < 3 , pscale=1;end
if nargin < 4 , selected_nodes=[];end
if nargin < 5 , GaussPts='pixels';end
if nargin < 6 , full_size=false;end
if nargin < 7 , deg=1;end
if deg>1, rint=0;end


xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
ng=1;ns=8*ones(2,1);
load(mesh_file,'-mat','conn','elt','ng','ns');
incp(1)=1;incp(2)=incp(1)*sizeim(1)*pscale;
if strcmp(GaussPts,'pixels') ,ngb=ng*prod(ns);ngt=ng*(1-0.75*(ng==4))*prod(ns)*2;ngq=ng*prod(ns);ngp=ng*prod(ns)*2;ngh=ng*prod(ns)/2;st=2;
elseif strcmp(GaussPts,'Gauss_points'),ngb=2; ngt=1+6*(deg==2);ngq=1+3*(~rint)+12*(deg==2);ngp=1+(~rint);ngh=1+7*(~rint);ng=1;ns=1;st=1;
elseif strcmp(GaussPts,'sub_cells'),ngb=2*prod(ns); ngt=(1+6*(deg==2))*prod(ns)*2;ngq=4*prod(ns);ngp=2*prod(ns)*2;ngh=8*prod(ns);ng=1;st=2;
end
if ng>0
    npts=2*ngb*sum(elt==2)+3*ngt*sum(elt==3)+4*ngq*sum(elt==4);
    nx=ngb*sum(elt==2)+ngt*sum(elt==3)+ngq*sum(elt==4);
else
    npts=mean(elt)*ceil(prod(Smesh))*pscale^2;
    nx=prod(sizeim*pscale);
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
        if ng==0
            xn=xo(inods);
            yn=yo(inods);
            [ypix,xpix]=meshgrid(ceil(min(yn)):floor(max(yn)),ceil(min(xn)):floor(max(xn)));
            
            if elt(i1)==3
                [xg,yg,wg]=GetGaussPointsTriangle(ngt/(prod(ns)/2),ns,xn,yn,xpix(:),ypix(:));
                N_r=[-1+0*xg,1+0*xg,0*yg];
                N_s=[-1+0*yg,0*xg,1+0*yg];
                
            elseif elt(i1)==2
                error('Not Coded yet');
            elseif elt(i1)==4
                [xg,yg,wg]=GetGaussPointsQuadrangle(ngq/(prod(ns)),ns,xn,yn,xpix(:),ypix(:));
                N_r=[-0.25*(1-yg),0.25*(1-yg),0.25*(1+yg),-0.25*(1+yg)];
                N_s=[-0.25*(1-xg),-0.25*(1+xg),0.25*(1+xg),0.25*(1-xg)];
            end
            Sel=length(xg);
            npix=1+incp(1)*(xpix-1)+incp(2)*(ypix-1);
            dxdr=N_r*xo(inods);
            dydr=N_r*yo(inods);
            dxds=N_s*xo(inods);
            dyds=N_s*yo(inods);
            detJ=(dxdr.*dyds-dydr.*dxds);
            invJ=[dyds./detJ,-dxds./detJ,-dydr./detJ,dxdr./detJ];
            for in=1:elt(i1)
                if isempty(selected_nodes)||full_size
                    id=inods(in);
                else
                    id=find(selected_nodes==inods(in));
                end
                if ~isempty(id)
                    N_x=N_r(:,in).*invJ(:,1)+N_s(:,in).*invJ(:,3);
                    N_y=N_r(:,in).*invJ(:,2)+N_s(:,in).*invJ(:,4);
                    indn(nel+(1:Sel))=id;
                    indp(nel+(1:Sel))=npix;
                    valx(nel+(1:Sel))=wg.*N_x(:);
                    valy(nel+(1:Sel))=wg.*N_y(:);
                    nel=nel+Sel;
                end
            end
            
            
        else
            if elt(i1)==2
                dxdr=Nb_r*xo(inods);
                dydr=Nb_r*yo(inods);
                detJ=(dxdr.^2+dydr.^2);
                invJ=[dxdr./detJ,dydr./detJ];
%                sg=-[1,-1];
                for in=1:2
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)
                         N_x=Nb_r(:,in).*invJ(:,1);
                         N_y=Nb_r(:,in).*invJ(:,2);
% dx=diff(xo(inods));
% dy=diff(yo(inods));
% if abs(dx)<eps,dx=Inf;end
% if abs(dy)<eps,dy=Inf;end
% 
%                         N_x=sg(in)*[ones(ngb,1)]/dx;
%                         N_y=sg(in)*[ones(ngb,1)]/dy;
                        indn(nel+(1:Selb))=id;
                        indp(nel+(1:Selb))=npix+(1:Selb);
                        valx(nel+(1:Selb))=N_x(:);
                        valy(nel+(1:Selb))=N_y(:);
                        nel=nel+Selb;
                    end
                end
                npix=npix+Selb;
            elseif elt(i1)==3
                dxdr=Nt_r*xo(inods);
                dydr=Nt_r*yo(inods);
                dxds=Nt_s*xo(inods);
                dyds=Nt_s*yo(inods);
                % if strcmp(GaussPts,'sub_cells')
                %     figure
                % plot(xo(inods),yo(inods),'+')
                %         Nt=[1-xgt-ygt,xgt,ygt];
                % hold on
                % plot(Nt*xo(inods),Nt*yo(inods),'o')
                % pause
                %
                % end
                detJ=(dxdr.*dyds-dydr.*dxds);
                invJ=[dyds./detJ,-dxds./detJ,-dydr./detJ,dxdr./detJ];
                for in=1:3
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)
                        N_x=Nt_r(:,in).*invJ(:,1)+Nt_s(:,in).*invJ(:,3);
                        N_y=Nt_r(:,in).*invJ(:,2)+Nt_s(:,in).*invJ(:,4);
                        indn(nel+(1:Selt))=id;
                        indp(nel+(1:Selt))=npix+(1:Selt);
                        valx(nel+(1:Selt))=N_x(:);
                        valy(nel+(1:Selt))=N_y(:);
                        nel=nel+Selt;
                    end
                end
                npix=npix+Selt;
            elseif elt(i1)==4
                dxdr=Nq_r*xo(inods);
                dydr=Nq_r*yo(inods);
                dxds=Nq_s*xo(inods);
                dyds=Nq_s*yo(inods);
                
                detJ=(dxdr.*dyds-dydr.*dxds);
                invJ=[dyds./detJ,-dxds./detJ,-dydr./detJ,dxdr./detJ];
                for in=1:4
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)
                        N_x=Nq_r(:,in).*invJ(:,1)+Nq_s(:,in).*invJ(:,3);
                        N_y=Nq_r(:,in).*invJ(:,2)+Nq_s(:,in).*invJ(:,4);
                        indn(nel+(1:Selq))=id;
                        indp(nel+(1:Selq))=npix+(1:Selq);
                        valx(nel+(1:Selq))=N_x(:);
                        valy(nel+(1:Selq))=N_y(:);
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
    valx=valx(1:nel);
    valy=valy(1:nel);
end

dphidx=sparse(indp,indn,valx,nx,ny);
dphidy=sparse(indp,indn,valy,nx,ny);



end
