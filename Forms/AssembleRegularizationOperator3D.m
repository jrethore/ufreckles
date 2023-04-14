function R=AssembleRegularizationOperator3D(nmod,iscale,select)
tic();
disp(sprintf('    Regularization operator...'));
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'selected','xo','yo','zo','conn','elt','Nnodes','Nelems');
load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,iscale-1)),'Nddl_tot');
if select
    disp(sprintf('    Tikhonov boundary regularization...'));
    Nddl=prod(Nnodes);
    ngt=1;ngq=4;ns=ones(2,1);
    if any(elt==6)
error('NOT CODED YET')
    end
    if any(elt==8)
%         [xgq,ygq,zgq,wgq]=GetGaussPointsHexaedron(ngq,ns);
%         Nq_r=[-0.125*(1-ygq).*(1-zgq),0.125*(1-ygq).*(1-zgq),0.125*(1+ygq).*(1-zgq),-0.125*(1+ygq).*(1-zgq),...
%             -0.125*(1-ygq).*(1+zgq),0.125*(1-ygq).*(1+zgq),0.125*(1+ygq).*(1+zgq),-0.125*(1+ygq).*(1+zgq)];
%         Nq_s=[-0.125*(1-xgq).*(1-zgq),-0.125*(1+xgq).*(1-zgq),0.125*(1+xgq).*(1-zgq),0.125*(1-xgq).*(1-zgq),...
%             -0.125*(1-xgq).*(1+zgq),-0.125*(1+xgq).*(1+zgq),0.125*(1+xgq).*(1+zgq),0.125*(1-xgq).*(1+zgq)];
%         Nq_t=[-0.125*(1-xgq).*(1-ygq),-0.125*(1+xgq).*(1-ygq),-0.125*(1+xgq).*(1+ygq),-0.125*(1-xgq).*(1+ygq),...
%             0.125*(1-xgq).*(1-ygq),0.125*(1+xgq).*(1-ygq),0.125*(1+xgq).*(1+ygq),0.125*(1-xgq).*(1+ygq)];
    [xgq,ygq,wgq]=GetGaussPointsQuadrangle(ngq);
    Selq=length(xgq);
    Nq_r=[-0.25*(1-ygq),0.25*(1-ygq),0.25*(1+ygq),-0.25*(1+ygq)];
    Nq_s=[-0.25*(1-xgq),-0.25*(1+xgq),0.25*(1+xgq),0.25*(1-xgq)];
    end
    ind1=[];
    ind2=[];
    val=[];

    for i1=1:prod(Nelems)
        inods=conn(i1,1:elt(i1));
        go=sum(~selected(inods));
        if go==4
            jnods=inods(~selected(inods));
            if elt(i1)==6
                error('NOT CODED YET');
            elseif elt(i1)==8
                dxdr=Nq_r*xo(jnods);
                dydr=Nq_r*yo(jnods);
                dzdr=Nq_r*zo(jnods);
                dxds=Nq_s*xo(jnods);
                dyds=Nq_s*yo(jnods);
                dzds=Nq_s*zo(jnods);
            end
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

            
            
            
            N_x=zeros(numel(DetJ),1);
            N_y=zeros(numel(DetJ),1);
            N_z=zeros(numel(DetJ),1);
            if elt(i1)==6
 error('NOT CODED YET');
            elseif elt(i1)==8
                for in=1:4
                    N_x(:,in)=Nq_r(:,in).*invJ(:,1)+Nq_s(:,in).*invJ(:,4);
                    N_y(:,in)=Nq_r(:,in).*invJ(:,2)+Nq_s(:,in).*invJ(:,5);
                    N_z(:,in)=Nq_r(:,in).*invJ(:,3)+Nq_s(:,in).*invJ(:,6);
                end

                wdetJo=diag(sparse(wgq.*DetJ));
            end
            Kel=N_x'*wdetJo*N_x+N_y'*wdetJo*N_y+N_z'*wdetJo*N_z;
            [indj,indi]=meshgrid(jnods,jnods);
            ind1=[ind1;indi(:)];
            ind2=[ind2;indj(:)];
            val=[val;Kel(:)];
            ind1=[ind1;indi(:)+Nddl];
            ind2=[ind2;indj(:)+Nddl];
            val=[val;Kel(:)];
            ind1=[ind1;indi(:)+2*Nddl];
            ind2=[ind2;indj(:)+2*Nddl];
            val=[val;Kel(:)];
            clear Ny Nx Nz
        end
    end
    R=sparse(ind1,ind2,val,Nddl_tot,Nddl_tot);
%     disp(sprintf('    Tikhonov boundary regularization...'));
%     Nddl=prod(Nnodes);
%     ngt=2;ngq=8;ns=ones(3,1);
%     if any(elt==6)
%         [xgt,ygt,zgt,wgt]=GetGaussPointsWedge(ngt,ns);
%         Nt_r=[-0.5*(1-zgt),0.5*(1-zgt),(0*ygt),...
%             -0.5*(1+zgt),0.5*(1+zgt),(0*ygt)];
%         Nt_s=[-0.5*(1-zgt),(0*xgt),0.5*(1-zgt),...
%             -0.5*(1-zgt),(0*xgt),0.5*(1-zgt)];
%         Nt_t=[-0.5*(1-xgt-ygt),-0.5*xgt,-0.5*ygt,...
%             0.5*(1-xgt-ygt),0.5*xgt,0.5*ygt];
%     end
%     if any(elt==8)
%         [xgq,ygq,zgq,wgq]=GetGaussPointsHexaedron(ngq,ns);
%         Nq_r=[-0.125*(1-ygq).*(1-zgq),0.125*(1-ygq).*(1-zgq),0.125*(1+ygq).*(1-zgq),-0.125*(1+ygq).*(1-zgq),...
%             -0.125*(1-ygq).*(1+zgq),0.125*(1-ygq).*(1+zgq),0.125*(1+ygq).*(1+zgq),-0.125*(1+ygq).*(1+zgq)];
%         Nq_s=[-0.125*(1-xgq).*(1-zgq),-0.125*(1+xgq).*(1-zgq),0.125*(1+xgq).*(1-zgq),0.125*(1-xgq).*(1-zgq),...
%             -0.125*(1-xgq).*(1+zgq),-0.125*(1+xgq).*(1+zgq),0.125*(1+xgq).*(1+zgq),0.125*(1-xgq).*(1+zgq)];
%         Nq_t=[-0.125*(1-xgq).*(1-ygq),-0.125*(1+xgq).*(1-ygq),-0.125*(1+xgq).*(1+ygq),-0.125*(1-xgq).*(1+ygq),...
%             0.125*(1-xgq).*(1-ygq),0.125*(1+xgq).*(1-ygq),0.125*(1+xgq).*(1+ygq),0.125*(1-xgq).*(1+ygq)];
%     end
%     ind1=[];
%     ind2=[];
%     val=[];
% 
%     for i1=1:prod(Nelems)
%         inods=conn(i1,1:elt(i1));
%         go=any(~selected(inods));
%         if go
%             if elt(i1)==6
%                 dxdr=Nt_r*xo(inods);
%                 dydr=Nt_r*yo(inods);
%                 dzdr=Nt_r*zo(inods);
%                 dxds=Nt_s*xo(inods);
%                 dyds=Nt_s*yo(inods);
%                 dzds=Nt_s*zo(inods);
%                 dxdt=Nt_t*xo(inods);
%                 dydt=Nt_t*yo(inods);
%                 dzdt=Nt_t*zo(inods);
%             elseif elt(i1)==8
%                 dxdr=Nq_r*xo(inods);
%                 dydr=Nq_r*yo(inods);
%                 dzdr=Nq_r*zo(inods);
%                 dxds=Nq_s*xo(inods);
%                 dyds=Nq_s*yo(inods);
%                 dzds=Nq_s*zo(inods);
%                 dxdt=Nq_t*xo(inods);
%                 dydt=Nq_t*yo(inods);
%                 dzdt=Nq_t*zo(inods);
%             end
%             DetJ =dxdr.*dyds.*dzdt + dxdt .*dydr.*dzds +...
%                 dxds.*dydt.*dzdr - dxdt .*dyds.*dzdr -...
%                 dxdr.*dydt.*dzds - dxds .*dydr.*dzdt;
%             invJ(:,1) = (dyds .* dzdt - dydt .* dzds) ./ DetJ;
%             invJ(:,4) = -(dydr .* dzdt - dydt .* dzdr) ./ DetJ;
%             invJ(:,7)  = (dydr .* dzds - dyds .* dzdr) ./ DetJ;
% 
%             invJ(:,2)  = -(dxds .* dzdt - dxdt .* dzds) ./ DetJ;
%             invJ(:,5)  = (dxdr .* dzdt - dxdt .* dzdr) ./ DetJ;
%             invJ(:,8)  = -(dxdr .* dzds - dxds .* dzdr) ./ DetJ;
% 
%             invJ(:,3)  = (dxds .* dydt - dxdt .* dyds) ./ DetJ;
%             invJ(:,6)  = -(dxdr .* dydt - dxdt .* dydr) ./ DetJ;
%             invJ(:,9)  = (dxdr .* dyds - dxds .* dydr) ./ DetJ;
%             N_x=zeros(numel(DetJ),1);
%             N_y=zeros(numel(DetJ),1);
%             N_z=zeros(numel(DetJ),1);
%             if elt(i1)==6
%                 for in=1:6
%                     N_x(:,in)=Nt_r(:,in).*invJ(:,1)+Nt_s(:,in).*invJ(:,4)+Nt_t(:,in).*invJ(:,7);
%                     N_y(:,in)=Nt_r(:,in).*invJ(:,2)+Nt_s(:,in).*invJ(:,5)+Nt_t(:,in).*invJ(:,8);
%                     N_z(:,in)=Nt_r(:,in).*invJ(:,3)+Nt_s(:,in).*invJ(:,6)+Nt_t(:,in).*invJ(:,9);
%                 end
% 
%                 wdetJo=diag(sparse(wgt.*DetJ));
%             elseif elt(i1)==8
%                 for in=1:8
%                     N_x(:,in)=Nq_r(:,in).*invJ(:,1)+Nq_s(:,in).*invJ(:,4)+Nq_t(:,in).*invJ(:,7);
%                     N_y(:,in)=Nq_r(:,in).*invJ(:,2)+Nq_s(:,in).*invJ(:,5)+Nq_t(:,in).*invJ(:,8);
%                     N_z(:,in)=Nq_r(:,in).*invJ(:,3)+Nq_s(:,in).*invJ(:,6)+Nq_t(:,in).*invJ(:,9);
%                 end
% 
%                 wdetJo=diag(sparse(wgq.*DetJ));
%             end
%             Kel=N_x'*wdetJo*N_x+N_y'*wdetJo*N_y+N_z'*wdetJo*N_z;
%             [indj,indi]=meshgrid(inods,inods);
%             ind1=[ind1;indi(:)];
%             ind2=[ind2;indj(:)];
%             val=[val;Kel(:)];
%             ind1=[ind1;indi(:)+Nddl];
%             ind2=[ind2;indj(:)+Nddl];
%             val=[val;Kel(:)];
%             ind1=[ind1;indi(:)+2*Nddl];
%             ind2=[ind2;indj(:)+2*Nddl];
%             val=[val;Kel(:)];
%             clear Ny Nx Nz
%         end
%     end
%     R=sparse(ind1,ind2,val,Nddl_tot,Nddl_tot);
else
    disp(sprintf('    Tikhonov regularization...'));
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'maskg');

    load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
    load(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,iscale-1)),'epsyy');
    load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,iscale-1)),'Uxy','Uyx');
    load(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,iscale-1)),'epszz');
    load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,iscale-1)),'Uxz','Uzx');
    load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,iscale-1)),'Uyz','Uzy');

    R= epsxx'*maskg*wdetJ*epsxx+Uxy'*maskg*wdetJ*Uxy+Uxz'*maskg*wdetJ*Uxz...
        +epsyy'*maskg*wdetJ*epsyy+Uyx'*maskg*wdetJ*Uyx+Uyz'*maskg*wdetJ*Uyz...
        +epszz'*maskg*wdetJ*epszz+Uzy'*maskg*wdetJ*Uzy+Uzx'*maskg*wdetJ*Uzx;
    if isfield(param,'enrichment')&&(iscale==1)
        load(fullfile('TMP',sprintf('%d_3d_xepsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
        load(fullfile('TMP',sprintf('%d_3d_xepsyy_%d',nmod,iscale-1)),'epsyy');
        load(fullfile('TMP',sprintf('%d_3d_xepsxy_%d',nmod,iscale-1)),'Uxy','Uyx');
        load(fullfile('TMP',sprintf('%d_3d_xepszz_%d',nmod,iscale-1)),'epszz');
        load(fullfile('TMP',sprintf('%d_3d_xepsxz_%d',nmod,iscale-1)),'Uzx','Uxz');
        load(fullfile('TMP',sprintf('%d_3d_xepsyz_%d',nmod,iscale-1)),'Uyz','Uzy');
        load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'masks');

        Rx= epsxx'*masks*wdetJ*epsxx+Uxy'*masks*wdetJ*Uxy+Uxz'*masks*wdetJ*Uxz...
            +epsyy'*masks*wdetJ*epsyy+Uyx'*masks*wdetJ*Uyx+Uyz'*masks*wdetJ*Uyz...
            +epszz'*masks*wdetJ*epszz+Uzy'*masks*wdetJ*Uzy+Uzx'*masks*wdetJ*Uzx;

        [indi,indj,val]=find(Rx);
        Nddls=size(R,1);
        keep=~((indi<=Nddls)&(indj<=Nddls));
        Rx=sparse(indi,indj,val.*keep,Nddl_tot,Nddl_tot);
        Ro=sparse(Nddl_tot-size(R,1),Nddl_tot-size(R,2));
        R=blkdiag(R,Ro);
        R=R+Rx;


    end
end
end
