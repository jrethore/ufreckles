function R=AssembleRegularizationOperator(cut_nodes,iscale,nmod,cond_type)
if nargin<4, cond_type=1;end
load(fullfile('TMP','params'));
param0=param;
tic();
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
disp(sprintf('    Regularization operator...'));
load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','conn','elt','Nnodes','Nelems');
Nddl_tot=2*prod(Nnodes);
face_elt=[];
if isfield(param,'enrichment')&&(iscale==1)
    load(fullfile('TMP',sprintf('%d_phix_%d',nmod,iscale-1)),'Nddl_tot');
    if ~iscell(param0.levelset_file)
        nbfis=1;
    else
        nbfis=numel(param0.levelset_file);
    end
    for ic=1:nbfis
        load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'face_elts');
        face_elt=[face_elt;face_elts];
    end
end
preallocate=isempty(cut_nodes);
npts=2*(3*3*sum(elt==3)+4*4*sum(elt==4));
if preallocate
    uncut=repmat(0,Nnodes);
    ind1=zeros(npts,1);
    ind2=zeros(npts,1);
    val=zeros(npts,1);
else
    uncut=repmat(1,Nnodes);
    uncut(cut_nodes)=0;
    ind1=[];
    ind2=[];
    val=[];
end
nind=0;
if isempty(unmasked_nodes)
    unmasked_nodes=1:prod(Nnodes);
else
    tmp=zeros(prod(Nnodes),1);
    tmp(unmasked_nodes)=1:numel(unmasked_nodes);
    unmasked_nodes=tmp;
    clear tmp
end
Nddl=max(unmasked_nodes);
ngt=1;ngq=4;ng=1;
if any(elt==3)
    [xgt,ygt,wgt]=GetGaussPointsTriangle(ngt);
    Selt=length(xgt);
    Nt_r=[-1+0*xgt,1+0*xgt,0*ygt];
    Nt_s=[-1+0*ygt,0*xgt,1+0*ygt];
end
if any(elt==4)
    [xgq,ygq,wgq]=GetGaussPointsQuadrangle(ngq);
    Selq=length(xgq);
    Nq_r=[-0.25*(1-ygq),0.25*(1-ygq),0.25*(1+ygq),-0.25*(1+ygq)];
    Nq_s=[-0.25*(1-xgq),-0.25*(1+xgq),0.25*(1+xgq),0.25*(1-xgq)];
end


for i1=1:prod(Nelems)
    if ~any(face_elt==i1)
        inods=conn(i1,1:elt(i1));
        switch cond_type
            case 1
                cond= (~any(uncut(inods)))&(mean(unmasked_nodes(inods)>0)==1);
            case 2
                cond= (any(~uncut(inods)))&(mean(unmasked_nodes(inods)>0)==1);
        end
        if cond
            if elt(i1)==3
                dxdr=Nt_r*xo(inods);
                dydr=Nt_r*yo(inods);
                dxds=Nt_s*xo(inods);
                dyds=Nt_s*yo(inods);
                
                detJ=(dxdr.*dyds-dydr.*dxds);
                invJ=[dyds./detJ,-dxds./detJ,-dydr./detJ,dxdr./detJ];
                N_x=zeros(size(Nt_r));
                N_y=zeros(size(Nt_r));
                for in=1:3
                    N_x(:,in)=Nt_r(:,in).*invJ(:,1)+Nt_s(:,in).*invJ(:,3);
                    N_y(:,in)=Nt_r(:,in).*invJ(:,2)+Nt_s(:,in).*invJ(:,4);
                end
                wdetJo=diag(sparse(wgt.*detJ));
            elseif elt(i1)==4
                dxdr=Nq_r*xo(inods);
                dydr=Nq_r*yo(inods);
                dxds=Nq_s*xo(inods);
                dyds=Nq_s*yo(inods);
                
                detJ=(dxdr.*dyds-dydr.*dxds);
                invJ=[dyds./detJ,-dxds./detJ,-dydr./detJ,dxdr./detJ];
                N_x=zeros(size(Nq_r));
                N_y=zeros(size(Nq_r));
                for in=1:4
                    N_x(:,in)=Nq_r(:,in).*invJ(:,1)+Nq_s(:,in).*invJ(:,3);
                    N_y(:,in)=Nq_r(:,in).*invJ(:,2)+Nq_s(:,in).*invJ(:,4);
                end
                wdetJo=diag(sparse(wgq.*detJ));
            end
            if 0
                eo=zeros(size(N_x));
                exx=[N_x,eo];eyy=[eo,N_y];exy=0.5*[N_y,N_x];
                Kel=exx'*(wdetJo*exx)+eyy'*(wdetJo*eyy)+2*exy'*(wdetJo*exy);
%                 ids=repmat(0,1,elt(i1));
%                 for in=1:elt(i1)
%                     ids(in)=find(unmasked_nodes==inods(in));
%                 end
                ids=unmasked_nodes(inods);
                ids=[ids,ids+Nddl];
                [indj,indi]=meshgrid(ids,ids);
                if preallocate
                    ipt=nind+(1:numel(indi));
                    nval=numel(ipt);
                    ind1(ipt)=indi(:);
                    ind2(ipt)=indj(:);
                    val(ipt)=Kel(:);
                    nind=nind+nval;
                else
                    ind1=[ind1;indi(:)];
                    ind2=[ind2;indj(:)];
                    val=[val;Kel(:)];
                    ind1=[ind1;indi(:)+Nddl];
                    ind2=[ind2;indj(:)+Nddl];
                    val=[val;Kel(:)];
                end
            else
                Kel=N_x'*wdetJo*N_x+N_y'*wdetJo*N_y;
%                 ids=repmat(0,1,elt(i1));
%                 for in=1:elt(i1)
%                     ids(in)=find(unmasked_nodes==inods(in));
%                 end
                ids=unmasked_nodes(inods);
                [indj,indi]=meshgrid(ids,ids);
                if preallocate
                    ipt=nind+(1:numel(indi));
                    nval=numel(ipt);
                    ind1(ipt)=indi(:);
                    ind2(ipt)=indj(:);
                    val(ipt)=Kel(:);
                    ind1(ipt+nval)=indi(:)+Nddl;
                    ind2(ipt+nval)=indj(:)+Nddl;
                    val(ipt+nval)=Kel(:);
                    nind=nind+2*nval;
                else
                    ind1=[ind1;indi(:)];
                    ind2=[ind2;indj(:)];
                    val=[val;Kel(:)];
                    ind1=[ind1;indi(:)+Nddl];
                    ind2=[ind2;indj(:)+Nddl];
                    val=[val;Kel(:)];
                end
            end
        end
    end
end
if (nind<npts)&&preallocate
    ind1=ind1(1:nind);
    ind2=ind2(1:nind);
    val=val(1:nind);
end
R=sparse(ind1,ind2,val,Nddl_tot,Nddl_tot);
end
