function [dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis3D(mesh_file,sizeim,pscale,selected_nodes,GaussPts,full_size,selected_elts)
rint=false;
[pp,filname,ext]=fileparts(mesh_file);
if isempty(ext)
    mesh_file=[mesh_file,'.mat'];
end
load(mesh_file,'-mat','rint','xo','yo','zo','Nnodes','Nelems','Smesh');

if nargin < 3 , pscale=1;end
if nargin < 4 , selected_nodes=[];end
if nargin < 5 , GaussPts='pixels';end
if nargin < 6 , full_size=false;end
if nargin < 7 , selected_elts=[];end


xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
zo=(zo-0.5)*pscale+0.5;

load(mesh_file,'-mat','conn','elt','ng','ns');
incp(1)=1;incp(2)=incp(1)*sizeim(1)*pscale;

if strcmp(GaussPts,'pixels') ,ngte=ng*prod(ns)*2;ngt=ng*prod(ns)*2;ngq=ng*prod(ns);st=1+(prod(ns)==1);
elseif strcmp(GaussPts,'Gauss_points'), ngte=1;ngt=1+(~rint);ngq=1+7*(~rint);ng=1;ns=1;st=1;
elseif strcmp(GaussPts,'sub_cells'), ngte=1*prod(ns)*2;ngt=2*prod(ns)*2;ngq=1*prod(ns);ng=1;st=1+(prod(ns)==1);
end
assert(ng>0)
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
    npts=mean(elt)*prod(Smesh)*pscale^2;nx=prod(sizeim*pscale);
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
    if any(elt==4)
        [xgte,ygte,zgte,wgte]=GetGaussPointsTetrahedron(ngte,ns);
        Selte=length(xgte);
        Nte_r=[-1+0*xgte,1+0*xgte,0*ygte,0*zgte];
        Nte_s=[-1+0*xgte,0*xgte,1+0*ygte,0*zgte];
        Nte_t=[-1+0*xgte,0*xgte,0*ygte,1+0*zgte];
%         Nte_r=Nte_r(:,[3,4,1,2]);
%         Nte_s=Nte_s(:,[3,4,1,2]);
%         Nte_t=Nte_t(:,[3,4,1,2]);
    end
    if any(elt==6)
        [xgt,ygt,zgt,wgt]=GetGaussPointsWedge(ngt/(2*prod(ns)),ns);
        Selt=length(xgt);
        Nt_r=[-0.5*(1-zgt),0.5*(1-zgt),(0*ygt),...
            -0.5*(1+zgt),0.5*(1+zgt),(0*ygt)];
        Nt_s=[-0.5*(1-zgt),(0*xgt),0.5*(1-zgt),...
            -0.5*(1-zgt),(0*xgt),0.5*(1-zgt)];
        Nt_t=[-0.5*(1-xgt-ygt),-0.5*xgt,-0.5*ygt,...
            0.5*(1-xgt-ygt),0.5*xgt,0.5*ygt];
    end
    if any(elt==8)
        [xgq,ygq,zgq,wgq]=GetGaussPointsHexaedron(ngq/prod(ns),ns);
        Selq=length(xgq);
        Nq_r=[-0.125*(1-ygq).*(1-zgq),0.125*(1-ygq).*(1-zgq),0.125*(1+ygq).*(1-zgq),-0.125*(1+ygq).*(1-zgq),...
            -0.125*(1-ygq).*(1+zgq),0.125*(1-ygq).*(1+zgq),0.125*(1+ygq).*(1+zgq),-0.125*(1+ygq).*(1+zgq)];
        Nq_s=[-0.125*(1-xgq).*(1-zgq),-0.125*(1+xgq).*(1-zgq),0.125*(1+xgq).*(1-zgq),0.125*(1-xgq).*(1-zgq),...
            -0.125*(1-xgq).*(1+zgq),-0.125*(1+xgq).*(1+zgq),0.125*(1+xgq).*(1+zgq),0.125*(1-xgq).*(1+zgq)];
        Nq_t=[-0.125*(1-xgq).*(1-ygq),-0.125*(1+xgq).*(1-ygq),-0.125*(1+xgq).*(1+ygq),-0.125*(1-xgq).*(1+ygq),...
            0.125*(1-xgq).*(1-ygq),0.125*(1+xgq).*(1-ygq),0.125*(1+xgq).*(1+ygq),0.125*(1-xgq).*(1+ygq)];
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
        if elt(i1)==4
            dxdr=Nte_r*xo(inods);
            dydr=Nte_r*yo(inods);
            dzdr=Nte_r*zo(inods);
            dxds=Nte_s*xo(inods);
            dyds=Nte_s*yo(inods);
            dzds=Nte_s*zo(inods);
            dxdt=Nte_t*xo(inods);
            dydt=Nte_t*yo(inods);
            dzdt=Nte_t*zo(inods);
            
           
            DetJ =dxdr.*dyds.*dzdt + dxdt .*dydr.*dzds +...
                dxds.*dydt.*dzdr - dxdt .*dyds.*dzdr -...
                dxdr.*dydt.*dzds - dxds .*dydr.*dzdt;
            invJ(:,1) = (dyds .* dzdt - dydt .* dzds) ./ DetJ;
            invJ(:,4) = -(dydr .* dzdt - dydt .* dzdr) ./ DetJ;
            invJ(:,7)  = (dydr .* dzds - dyds .* dzdr) ./ DetJ;

            invJ(:,2)  = -(dxds .* dzdt - dxdt .* dzds) ./ DetJ;
            invJ(:,5)  = (dxdr .* dzdt - dxdt .* dzdr) ./ DetJ;
            invJ(:,8)  = -(dxdr .* dzds - dxds .* dzdr) ./ DetJ;

            invJ(:,3)  = (dxds .* dydt - dxdt .* dyds) ./ DetJ;
            invJ(:,6)  = -(dxdr .* dydt - dxdt .* dydr) ./ DetJ;
            invJ(:,9)  = (dxdr .* dyds - dxds .* dydr) ./ DetJ;
            % tic
            for in=1:4
                if isempty(selected_nodes)||full_size
                    id=inods(in);
                else
                    id=find(selected_nodes==inods(in));
                end
                if ~isempty(id)
                    N_x=Nte_r(:,in).*invJ(:,1)+Nte_s(:,in).*invJ(:,4)+Nte_t(:,in).*invJ(:,7);
                    N_y=Nte_r(:,in).*invJ(:,2)+Nte_s(:,in).*invJ(:,5)+Nte_t(:,in).*invJ(:,8);
                    N_z=Nte_r(:,in).*invJ(:,3)+Nte_s(:,in).*invJ(:,6)+Nte_t(:,in).*invJ(:,9);
                    ii=nel+(1:Selte);
                    indn(ii)=id;
                    indp(ii)=npix+(1:Selte);
                    valx(ii)=N_x(:);
                    valy(ii)=N_y(:);
                    valz(ii)=N_z(:);
                    nel=nel+Selte;
                end
            end
        elseif elt(i1)==6
            dxdr=Nt_r*xo(inods);
            dydr=Nt_r*yo(inods);
            dzdr=Nt_r*zo(inods);
            dxds=Nt_s*xo(inods);
            dyds=Nt_s*yo(inods);
            dzds=Nt_s*zo(inods);
            dxdt=Nt_t*xo(inods);
            dydt=Nt_t*yo(inods);
            dzdt=Nt_t*zo(inods);
            for ig=1:ngt
                J=[dxdr(ig),dxds(ig),dxdt(ig);...
                    dydr(ig),dyds(ig),dydt(ig);...
                    dzdr(ig),dzds(ig),dzdt(ig)];
                for in=1:6
                    if isempty(selected_nodes)||full_size
                        id=inods(in);
                    else
                        id=find(selected_nodes==inods(in));
                    end
                    if ~isempty(id)
                        N_rst=[Nt_r(ig,in);Nt_s(ig,in);Nt_t(ig,in)];
                        N_xyz=J'\N_rst;
                        indn(nel+1)=id;
                        indp(nel+1)=npix+ig;
                        valx(nel+1)=N_xyz(1);
                        valy(nel+1)=N_xyz(2);
                        valz(nel+1)=N_xyz(3);
                        nel=nel+1;
                    end
                end
            end
            %            npix=npix+Selt;
        elseif elt(i1)==8
            dxdr=Nq_r*xo(inods);
            dydr=Nq_r*yo(inods);
            dzdr=Nq_r*zo(inods);
            dxds=Nq_s*xo(inods);
            dyds=Nq_s*yo(inods);
            dzds=Nq_s*zo(inods);
            dxdt=Nq_t*xo(inods);
            dydt=Nq_t*yo(inods);
            dzdt=Nq_t*zo(inods);

            DetJ =dxdr.*dyds.*dzdt + dxdt .*dydr.*dzds +...
                dxds.*dydt.*dzdr - dxdt .*dyds.*dzdr -...
                dxdr.*dydt.*dzds - dxds .*dydr.*dzdt;
            invJ(:,1) = (dyds .* dzdt - dydt .* dzds) ./ DetJ;
            invJ(:,4) = -(dydr .* dzdt - dydt .* dzdr) ./ DetJ;
            invJ(:,7)  = (dydr .* dzds - dyds .* dzdr) ./ DetJ;

            invJ(:,2)  = -(dxds .* dzdt - dxdt .* dzds) ./ DetJ;
            invJ(:,5)  = (dxdr .* dzdt - dxdt .* dzdr) ./ DetJ;
            invJ(:,8)  = -(dxdr .* dzds - dxds .* dzdr) ./ DetJ;

            invJ(:,3)  = (dxds .* dydt - dxdt .* dyds) ./ DetJ;
            invJ(:,6)  = -(dxdr .* dydt - dxdt .* dydr) ./ DetJ;
            invJ(:,9)  = (dxdr .* dyds - dxds .* dydr) ./ DetJ;
            % tic
            for in=1:8
                if isempty(selected_nodes)||full_size
                    id=inods(in);
                else
                    id=find(selected_nodes==inods(in));
                end
                if ~isempty(id)
                    N_x=Nq_r(:,in).*invJ(:,1)+Nq_s(:,in).*invJ(:,4)+Nq_t(:,in).*invJ(:,7);
                    N_y=Nq_r(:,in).*invJ(:,2)+Nq_s(:,in).*invJ(:,5)+Nq_t(:,in).*invJ(:,8);
                    N_z=Nq_r(:,in).*invJ(:,3)+Nq_s(:,in).*invJ(:,6)+Nq_t(:,in).*invJ(:,9);
                    ii=nel+(1:Selq);
                    indn(ii)=id;
                    indp(ii)=npix+(1:Selq);
                    valx(ii)=N_x(:);
                    valy(ii)=N_y(:);
                    valz(ii)=N_z(:);
                    nel=nel+Selq;
                end
            end


            %             toc
            %             tic

            %             for ig=1:ngq
            %                 J=[dxdr(ig),dxds(ig),dxdt(ig);...
            %                    dydr(ig),dyds(ig),dydt(ig);...
            %                    dzdr(ig),dzds(ig),dzdt(ig)];
            %                 for in=1:8
            %                 if isempty(selected_nodes)||full_size
            %                     id=inods(in);
            %                 else
            %                     id=find(selected_nodes==inods(in));
            %                 end
            %                 if ~isempty(id)
            %                     N_rst=[Nq_r(ig,in);Nq_s(ig,in);Nq_t(ig,in)];
            %                     N_xyz=J\N_rst;
            %                     indn(nel+1)=id;
            %                     indp(nel+1)=npix+ig;
            %                     valx(nel+1)=N_xyz(1);
            %                     valy(nel+1)=N_xyz(2);
            %                     valz(nel+1)=N_xyz(3);
            %                     nel=nel+1;
            %                 end
            %                 end
            %             end
            %            toc
            %            npix=npix+Selq;
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
    valx=valx(1:nel);
    valy=valy(1:nel);
    valz=valz(1:nel);
end

dphidx=sparse(indp,indn,valx,nx,ny);
dphidy=sparse(indp,indn,valy,nx,ny);
dphidz=sparse(indp,indn,valz,nx,ny);



end
