function [lso,ls1,nmesh]=ComputeLevelSetFromPoints(nmod,iscale,xon,yon,nxon,nyon,refine)
check=1;
if nargin<7,refine=0;end
load(fullfile('TMP','params'),'param');
param0=param;
tau=param0.transition_length;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
thickness=0;
if isfield(param,'line_thickness')
    thickness=round(0.5*(param.line_thickness));
end
roi=param0.roi;
sizeim=[roi(2)-roi(1)+1,roi(4)-roi(3)+1];
mesh_size=round(length(xon)/max(20,param.Nelems))*ones(1,2);
xo=round(0:mesh_size(1):sizeim(1)-1)+1;
yo=round(0:mesh_size(2):sizeim(2)-1)+1;
xo=xo-floor(mean(xo)-sizeim(1)/2);
yo=yo-floor(mean(yo)-sizeim(2)/2);
Nnodes=[length(xo),length(yo),1];
Nelems=max(Nnodes-1,1);
[yo,xo]=meshgrid(yo,xo);
yo=yo(:)-0.5;
xo=xo(:)-0.5;
incn=[0,1,Nnodes(1)+1,Nnodes(1)];
nroot=reshape(1:prod(Nnodes),Nnodes);
nroot=nroot(1:max(1,Nnodes(1)-1),:);
nroot=nroot(:,1:max(1,Nnodes(2)-1));
conn=repmat(nroot(:),1,4)+repmat(incn,prod(Nelems),1);
elt=repmat(4,prod(Nelems),1);
dX=repmat(xo',length(xon),1)-repmat(xon,1,length(xo));
dY=repmat(yo',length(xon),1)-repmat(yon,1,length(xo));
[lson,ind]=min(abs(dX+1i*dY));
lson=lson';
if check
    figure
    plot(xo,yo,'b+')
    hold on
end
nbandn=lson<(tau*2^(iscale-1)+thickness+mean(mesh_size));
elts=GetEltsFromNodes(conn,elt,find(nbandn));
elt=elt(elts);
conn=conn(elts,:);
[nbandn,nconn]=GetNodesFromElts(conn,elt);
xo=xo(nbandn);
yo=yo(nbandn);
lson=lson(nbandn);
ind=ind(nbandn);
conn=nconn;
Nnodes=[length(xo),1,1];
Nelems=[length(elt),1,1];

[elt,conn,xo,yo,Nnodes,Nelems]=BuiltRefinedConnectivity(elt,conn,xo,yo,ones(length(xo),1),refine);
xo=xo+0.1*(abs(xo-round(xo))==0);
yo=yo+0.1*(abs(yo-round(yo))==0);
mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
ng=0;
ns=round(mesh_size/2);
Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
rint=false;
zo=1;
save(mesh_file,'xo','yo','zo','Nnodes','Nelems','elt','conn','ng','ns','Smesh','rint');
if refine
    dX=repmat(xo',length(xon),1)-repmat(xon,1,length(xo));
    dY=repmat(yo',length(xon),1)-repmat(yon,1,length(xo));
    [lson,ind]=min(abs(dX+1i*dY));
    lson=lson';
    
    
end

phip=CreateFiniteElementBasis(mesh_file,sizeim,1,[],'pixels');
nmesh=sum(phip,2)>0;
phi=CreateFiniteElementBasis(mesh_file,sizeim,1,[],'Gauss_points');
[dphidx,dphidy]=CreateGradFiniteElementBasis(mesh_file,sizeim,1,[],'Gauss_points');
[wdetJ]=GetWeigthDetJ(mesh_file,sizeim,1,'Gauss_points');
xg=phi*xo;
yg=phi*yo;

if check
    plot(xo,yo,'ro')
%     figure
%     imagesc(reshape(sum(phip,2),sizeim)')
%     invgray=colormap('gray');
%     invgray=flipud(invgray);
%     colormap(invgray);
%     axis off;
%     axis xy;
%     axis image;
%     hold on
%     plot(xon,yon,'w.','MarkerSize',0.25)
%     plot(xo,yo,'g.')
%     print ('-depsc2', [param0.result_file,'-mesh_adapt.eps']);
    
    figure
    hold on
    for ie=1:length(elt)
        inods=conn(ie,1:elt(ie));
        inods=[inods,inods(1)];
        plot(xo(inods),yo(inods),'-r','LineWidth',2)
    end
    plot(xon,yon,'k.','MarkerSize',5)
    axis off;
    axis xy;
    axis image;
    print ('-depsc2', [param0.result_file,'-mesh_adapt.eps']);
    
    figure
    hold on
    imagesc(reshape(phip*lson,sizeim)')
    axis off;
    axis xy;
    axis image;
    figure
    hold on
    imagesc(reshape(phip*(1+0*lson),sizeim)')
    axis off;
    axis xy;
    axis image;
end
lson=-lson.*sign((xo-xon(ind)).*nxon(ind)+(yo-yon(ind)).*nyon(ind));
lso=reshape(phip*lson,sizeim);

% lso=LSReinit(lso,20,1);
% Mp=phip'*phip;
% Fp=phip'*lso(:);
% lson=Mp\Fp;
if check
    mask=sum(phip,2)==0;
    lsp=phip*lson;
    lsp(mask)=max(lsp(:));
    figure
    cb=colormap(hot);
    cb=cb(:,[3,2,1]);
    colormap(cb);
    hold on
    imagesc(reshape(lsp,sizeim)')
    axis off;
    axis xy;
    axis image;
    plot(xon,yon,'w.','MarkerSize',0.5)
%    colorbar('FontSize',20)
    print ('-depsc2', [param0.result_file,'-ls0.eps']);
end
nxg=dphidx*lson;
nyg=dphidy*lson;
nn=abs(nxg+i*nyg);
nxg=nxg./nn;
nyg=nyg./nn;
%     Nx=phi\nxg;
%     Ny=phi\nyg;
%     R=min(max(dphidx*Nx+dphidy*Ny,1/(10*max(sizeim))),100/min(sizeim));
% R=(dphidx*Nx+dphidy*Ny);

enx=-nyg;
eny=nxg;
M=dphidx'*wdetJ*dphidx+dphidy'*wdetJ*dphidy;
%    F=dphidx'*wdetJ*(enx.*R)+dphidy'*wdetJ*(eny.*R);
F=dphidx'*wdetJ*(enx)+dphidy'*wdetJ*(eny);

if param.closed
    
    xg=mean(xo);yg=mean(yg);
    if isfield(param,'starting_point')
        xc=param.starting_point(1)-roi(1)+1;
        yc=param.starting_point(2)-roi(3)+1;
        tfissx=xg-xc;tfissy=yg-yc;
        tnorm=abs(tfissx+i*tfissy);
        tfissx=tfissx/tnorm;tfissy=tfissy/tnorm;
    else
        xc=xon(1);
        yc=yon(1);
        tfissx=nxon(1);
        tfissy=nyon(1);
    end
    inv=sign(tfissx*(xg-xc)+tfissy*(yg-yc));
    if inv<0
        tfissx=-tfissx;
        tfissy=-tfissy;
    end
    nfissx=-tfissy;
    nfissy=tfissx;
    xg=xc+tfissx*(max(lso(:))+2*mean(mesh_size));
    yg=yc+tfissy*(max(lso(:))+2*mean(mesh_size));
    crack=(xo-xg)*nfissx+(yo-yg)*nfissy;
    front=(xo-xg)*tfissx+(yo-yg)*tfissy;
    hnodes=2*double(crack>=0)-1;
    %         figure
    %         imagesc(reshape(phip*crack,sizeim))
    %         figure
    %         imagesc(reshape(phip*front,sizeim))
    %%
    enriched=zeros(length(xo),1);
    face_elts=[];
    for ie=1:length(elt)
        inods=conn(ie,1:elt(ie));
        infront_behind_0=mean(front(inods)>=0);
        above_below=abs(mean(hnodes(inods)));
        if (above_below<1)&&(infront_behind_0<1)
            enriched(inods)=1;
            face_elts=[face_elts,ie];
        end
    end
    face_nodes=find(enriched);
    if check
        figure
        plot(xo,yo,'b+')
        hold on
        plot(xo(face_nodes),yo(face_nodes),'ro')
    end
    
    wdetJs=GetWeigthDetJ(mesh_file,sizeim,1,'sub_cells',face_elts);
    phipe=CreateFiniteElementBasis(mesh_file,sizeim,1,face_nodes,'pixels',true);
    phis=CreateFiniteElementBasis(mesh_file,sizeim,1,face_nodes,'sub_cells',true);
    [dphidxs,dphidys]=CreateGradFiniteElementBasis(mesh_file,sizeim,1,face_nodes,'sub_cells',true);
    
    nxs=dphidxs*lson;
    nys=dphidys*lson;
    nn=abs(nxs+1i*nys);
    nxs=nxs./nn;
    nys=nys./nn;
    %     Rs=min(max(dphidxs*Nx+dphidys*Ny,1/(10*max(sizeim))),100/min(sizeim));
    %      Rs=(dphidxs*Nx+dphidys*Ny);
    
    
    enx=nys;
    eny=-nxs;
    
    hn=double(crack>=0);
    hn=diag(sparse(hn(face_nodes)));
    heaviside=double(phis*crack>=0);
    heaviside=diag(sparse(heaviside));
    dphihdx=heaviside*dphidxs(:,face_nodes)-dphidxs(:,face_nodes)*hn;
    dphihdy=heaviside*dphidys(:,face_nodes)-dphidys(:,face_nodes)*hn;
    dphidxs=[dphidxs,dphihdx];
    dphidys=[dphidys,dphihdy];
    
    
    
    heaviside=double(phip*crack>=0);
    heaviside=diag(sparse(heaviside));
    phipe=heaviside*phipe(:,face_nodes)-phipe(:,face_nodes)*hn;
    
    Ms=dphidxs'*wdetJs*dphidxs+dphidys'*wdetJs*dphidys;
    %        Fs=dphidxs'*wdetJs*(enx.*Rs)+dphidys'*wdetJs*(eny.*Rs);
    Fs=dphidxs'*wdetJs*(enx)+dphidys'*wdetJs*(eny);
    
    [indi,indj,val]=find(Ms);
    Nddl_tot=size(Ms,1);
    Nddls=size(M,1);
    keep=~((indi<=Nddls)&(indj<=Nddls));
    Ms=sparse(indi,indj,val.*keep,Nddl_tot,Nddl_tot);
    keep=~((1:Nddl_tot)'<=Nddls);
    Fs=keep.*Fs;
    if Nddl_tot>size(M,2)
        Mo=sparse(Nddl_tot-size(M,1),Nddl_tot-size(M,2));
        M=blkdiag(M,Mo);
        Fo=zeros(Nddl_tot-size(F,1),1);
        F=[F;Fo];
    end
    M=M+Ms;
    F=F+Fs;
    phip=[phip,phipe];
end

indn=1;
indi=[indn];
indj=1:length(indi);
C=sparse(indi,indj,1,size(M,1),length(indi));
Mt=[M,C;C',sparse(size(C,2),size(C,2))];
Lt=Mt\[F;zeros(size(C,2),1)];
ls1n=Lt(1:length(F));
ls1=reshape(phip*ls1n,sizeim);
ls1=ls1-min(ls1(:));
on=sub2ind(sizeim,round(xon),round(yon));
ls1=min(max(ls1,min(ls1((on)))),max(ls1((on))));
if check
    lsp=ls1;
    lsp(mask)=max(lsp(:));
    figure
    cb=colormap(hot);
    cb=cb(:,[3,2,1]);
    colormap(cb);
    hold on
    imagesc(reshape(lsp,sizeim)')
    axis off;
    axis xy;
    axis image;
    plot(xon,yon,'w.','MarkerSize',0.5)
%    colorbar('FontSize',20)
    print ('-depsc2',[param0.result_file, '-ls1.eps']);
end

end