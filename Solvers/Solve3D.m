function [U]=Solve3D(Uini,nmod,iscale)
clear MedianFilter
if nargin<3,iscale=1;end
pscale=2^(iscale-1);
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
[~,filreso,~]=fileparts(param0.result_file);
save(fullfile('TMP',sprintf('%s-scale-%d-restart',filreso,iscale)),'param');
roi=param0.roi;
if iscell(param0.deformed_image)
    nim=length(param0.deformed_image);
else
    nim=1;
end
if isfield(param0,'restart')
    restart=param0.restart;
else
    restart=1;
end
if restart
    nmax=nim;
else
    if iscale==1
        nmax=nim;
    else
        if nim<2
            nmax=nim;
        else
            nmax=2;
        end
    end
end
inorm=false;
if isfield(param0,'normalize_grey_level')
    inorm=param0.normalize_grey_level;
end
reg_type='none';
mfilter=0;
if isfield(param0,'regularization_type')
    reg_type=param0.regularization_type;
    if iscale>1&&~(strcmp(reg_type,'none'))
        reg_type='tiko';
    end
    if strcmp(param0.regularization_type,'median')
        reg_type='median';
        mfilter=1;
        lmed=1;
        if isfield(param0,'regularization_parameter')
            lmed=param0.regularization_parameter;
        end
    elseif strcmp(param0.regularization_type,'tiko')
        mfilter=1;lmed=1;
    end
end
maxiter=param0.iter_max;
conv=param0.convergance_limit;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
disp(sprintf('Starting resolution for scale %d...',iscale));
mesh_file=fullfile('TMP',sprintf('%d_3d_mesh_%d.mat',nmod,iscale-1));
load(mesh_file,'rflag','rint','xo','yo','zo','Nnodes','Nelems','Smesh','conn','elt','ng','ns');
invmap=~rflag;
pstep=1;
if isfield(param0,'psample')
    pstep=param0.psample;
end
face_elts=[];
face_nodes=[];
if isfield(param,'enrichment')&&(iscale==1)
    if iscell(param.levelset_file)
        error('not coded yet');
    else
        ncrack=1;
    end
    for ic=1:ncrack
        if ncrack==1,lvl7file=param.levelset_file;end
        load(lvl7file,'crack','front','zone');
        [iface_nodes,iface_elts,crackn]=TreatmentofEnrichment3D(ic,nmod,crack,front,zone);
        clear front
        face_elts=[face_elts,iface_elts];
        face_nodes=[face_nodes,iface_nodes];
    end
    if size(Uini,1)==(3*prod(Nnodes))
        Uini=[Uini;zeros(3*numel(face_nodes),size(Uini,2))];
    end
    
end
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
switch param0.stack_format
    case 'bin'
        fidi=fopen(fullfile('TMP',sprintf('dsample0_%d',iscale-1)));
        im0=fread(fidi,prod(sizeim));
        fclose(fidi);
        im0=reshape(im0,sizeim);
    case 'mat'
        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
end
[mean0,std0]=mexImScalling3D(int64(im0));
mean1=0;std1=1;
dynamic=max(im0(:))-min(im0(:));
F=zeros(3*prod(Nnodes)+3*numel(face_nodes),1);
M=0;
im1=0;
restart_file='';
if isfield(param0,'restart_file')
    restart_file=param0.restart_file;
end
if isempty(restart_file)
    imin=1;
    U=Uini;
else
    load(param0.restart_file,'U','ijm')
    imin=ijm;
    Uini=U;
end
if ~strcmp(reg_type,'none')&&~strcmp(reg_type,'median')
    lc=(min(sizeim)/10);
    lm=lc;
    if isfield(param0,'regularization_parameter')
        lm=param0.regularization_parameter/pscale;
    end
    ki=1/lc;
    V=repmat(cos(2*pi*(xo)*ki).*cos(2*pi*(yo)*ki).*cos(2*pi*(zo)*ki),3,1);
    if isfield(param,'enrichment')&&(iscale==1)
        V=[V;zeros(3*numel(face_nodes),1)];
    end
        hh=min(sqrt(diff(xo(conn(:,1:2)),[],2).^+diff(yo(conn(:,1:2)),[],2).^2+diff(zo(conn(:,1:2)),[],2).^2));
    if param0.regularization_parameter>3*hh
        mfilter=0;
    end
end
fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
if dynamic>255
    fwrite(fid,255);
else
    fwrite(fid,dynamic);
end
model=param;
save(fullfile('TMP',sprintf('%s-scale-%d-restart',filreso,iscale)),'model','rint','Nnodes','Nelems','xo','yo','zo','conn','elt','ng','ns','-append');
for ijm=imin:nmax
    
    display(sprintf('Image %d/%d',ijm,nim));
    res=1;
    ii=1;
    if ~restart&&ijm>2
        Ui=2*U(:,ijm-1)-U(:,ijm-2);
    else
        Ui=Uini(:,ijm);
    end
    merrorp=Inf;
    while ( res>conv && ii< maxiter)
        [merror,disc]=Assemble(ijm);
        if (ii==1)&&(ijm==imin)
            switch reg_type
                case 'tiko'
                    CreateGradBasisFunction3D(iscale,nmod);
                    R=AssembleRegularizationOperator3D(nmod,iscale,isfield(param,'coupling_parameter'));
                    a=(V'*M*V)/(V'*R*V);
                    a=a*(2*lm/lc)^2;
                    R=a*R;
                    P=sparse(size(M,1),size(M,2));
                case 'equilibrium_gap'
                    CreateGradBasisFunction3D(iscale,nmod);
                    Rb=AssembleRegularizationOperator3D(nmod,iscale,1);
                    D=diag(abs(diag(Rb))>0);
                    a=(V'*D*M*V)/(V'*Rb*V);
                    a=a*(2*lm/lc)^2;
                    Rb=a*Rb;
                    disp(sprintf('    Equilibrium Gap regularization...'));
                    LoadMat(nmod);
                    AssembleMechanicalOperator(iscale,nmod,ijm)
                    AssembleEquilibriumGapOperator(iscale,nmod,1)
                    load(fullfile('TMP',sprintf('%d_kk_operator_%d',nmod,iscale-1)),'K');
                    R=K;
                    clear K
                    a=(V'*M*V)/(V'*R*V);
                    a=a*(2*lm/lc)^4;
                    R=a*R;
                    R=R+Rb;
                    P=sparse(size(M,1),size(M,2));
                case 'constitutive_gap'
                    disp(sprintf('    Constitutive Gap regularization...'));
                    CreateGradBasisFunction3D(iscale,nmod);
                    LoadMat(nmod);
                    AssembleMechanicalOperator(iscale,nmod,ijm)
                    AssembleEquilibriumGapOperator(iscale,nmod,0)
                    load(fullfile('TMP',sprintf('%d_k_operator_%d',nmod,iscale-1)),'K');
                    load(fullfile('TMP',sprintf('%d_kk_operator_%d',nmod,iscale-1)),'select');
                    R=K;
                    clear K
                    P=diag(1.e9*sparse((~diag(select))).*diag(R));
                    a=(V'*M*V)/(V'*R*V);
                    a=a*(2*lm/lc)^2;
                    R=a*R;
                    P=a*P;
                otherwise
                    R=sparse(size(M,1),size(M,2));
                    P=sparse(size(M,1),size(M,2));
                    
            end
            M=M+R+P;
            try
                M1 = ichol(M,struct('type','ict','droptol',1e-3));
                M2 = M1';
            catch err
                disp([err.identifier ' : ' err.message ' => using diagonal compensation'])
                diagcomp=((sum(abs(M),2)./diag(M))-2)/100;
                try
                    M1 = ichol(M,struct('type','ict','droptol',1e-3,'diagcomp',max(diagcomp)));
                    M2 = M1';
                catch err
                    disp([err.identifier ' : ' err.message ' => removing zero diag'])
                    ido=isnan(diagcomp);
                    M(ido,ido)=1;
                    disp([err.identifier ' : ' err.message ' => using lumped operator'])
                    M=full(sum(M,2));
                    M1=[];
                    M2=[];
                end
            end
        end
        F=F-R*Ui;
        if issparse(M)
            F=F+M*Ui;
            dU=pcg(M,F,1e-6,1000,M1,M2,Ui);
            dU=dU-Ui;
        else
            dU=F./M;
        end
        disp(sprintf('At iteration # %d',ii));
        disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100/dynamic));
        nU0=max(1,norm(Ui+dU(:)));
        res=norm(dU(:))/nU0;
        disp(sprintf('|dU|=%f',res));
        if mfilter
            res=min(res,abs(merrorp-merror)/abs(merror));
        end
        Ui=Ui+dU;
        ii=ii+1;
        if mfilter
            Ui=MedianFilter(Ui,mesh_file,lmed);
        end
        merrorp=merror;
        U(:,ijm)=Ui;
        save(fullfile('TMP',sprintf('%s-scale-%d-restart',filreso,iscale)),'U','ijm','-append');
    end
%     if mfilter
%         Uii=zeros(size(Ui));
%         Nnn=prod(Nnodes);
%         for ii1=1:prod(Nelems)
%             if mfilter
%                 for iin=1:elt(ii1)
%                     nod=conn(ii1,iin);
%                     if abs(Uii(nod))==0
%                         ee=sum(conn==nod,2)>0;
%                         ee=conn(ee,:);
%                         ee=unique(ee(:));
%                         Uii(nod)=median(Ui(ee));
%                         Uii(nod+Nnn)=median(Ui(ee+Nnn));
%                         Uii(nod+2*Nnn)=median(Ui(ee+2*Nnn));
%                     end
%                 end
%                 Ui=Uii;
%                 U(:,ijm)=Ui;
%                 save(fullfile('TMP',sprintf('%s-scale-%d-restart',filreso,iscale)),'U','-append');
%             end
%         end
%     end
    if dynamic>255
        disc=255*disc/dynamic;
    end
    fwrite(fid,abs(disc));
    
end
disp(sprintf('Enlapsed time for resolution = %6.2f s',cputime -ttic));
fclose(fid);
    function [mdisc,disc]=Assemble(kkk)
        mdisc=0;
        F=0*F;
        %        Fels=zeros(size(conn,1),3*size(conn,2));
        roip=(roi-1)/pscale+1;
%        roip(1:2:end)=floor(roip(1:2:end));
%        roip(2:2:end)=ceil(roip(2:2:end));
        if (kkk==imin)&&(ii==1)
            do_matrix=true;
            nind=(3*6)^2*sum(elt==6)+(3*8)^2*(sum(elt==8)-numel(face_elts))+(2*3*8)^2*numel(face_elts);
            indig=zeros(nind,1);
            indjg=zeros(nind,1);
            val=zeros(nind,1);
            disc=zeros(length(elt),1);
            nel=0;
        else
            do_matrix=false;
        end
        Nn=prod(Nnodes);
        nenr=numel(face_nodes);
        if ii==1
            if nim==1
                fildef=param0.deformed_image;
            else
                fildef=param0.deformed_image{kkk};
            end
            clear im1
            if isfield(param0,'stack_size')
                imsiz0=param0.stack_size;
                fidi=fopen(fildef,'r');
                im1=fread(fidi,prod(imsiz0));
                im1=reshape(im1,imsiz0);
                fclose(fidi);
            else
                [~, ~, ext] = fileparts(fildef);
                switch ext
                    case '.mat'
                        load(fildef,'jm3');
                        im1=(jm3);
                        clear jm3
                    case {'.tif','.tiff'}
                        im1=readTIFFasRAW(fildef);
                end
            end
            if iscale>1
                for ip=1:(iscale-1)
                    scale=2;
                    imsiz0=size(im1);
                    imsiz1=floor(imsiz0/2);
                    nn=2*imsiz1;
                    if nn(1)<imsiz0(1)
                        im1((nn(1)+1):end,:,:)=[];
                    end
                    if nn(2)<imsiz0(2)
                        im1(:,(nn(2)+1):end,:)=[];
                    end
                    if nn(3)<imsiz0(3)
                        im1(:,:,(nn(3)+1):end)=[];
                    end
                    
                    im1=reshape(im1,scale,prod(nn)/scale);
                    im1=mean(im1,1);
                    nn(1)=nn(1)/scale;
                    im1=reshape(im1,nn);
                    
                    im1=permute(im1,[2,3,1]);
                    im1=reshape(im1,scale,prod(nn)/scale);
                    im1=mean(im1,1);
                    nn(2)=nn(2)/scale;
                    im1=reshape(im1,nn([2,3,1]));
                    im1=permute(im1,[3,1,2]);
                    
                    im1=permute(im1,[3,1,2]);
                    im1=reshape(im1,scale,prod(nn)/scale);
                    im1=mean(im1,1);
                    nn(3)=nn(3)/scale;
                    im1=reshape(im1,nn([3,1,2]));
                    im1=permute(im1,[2,3,1]);
                end
            end
            Urbti=round(mean(Ui(0*Nn+(1:Nn))));
            Vrbti=round(mean(Ui(1*Nn+(1:Nn))));
            Wrbti=round(mean(Ui(2*Nn+(1:Nn))));
            [mean1,std1]=mexImScalling3D(int64(im1(max(1,Urbti+roip(1)):min(size(im1,1),Urbti+roip(2)),...
                max(1,Vrbti+roip(3)):min(size(im1,2),Vrbti+roip(4)),...
                max(1,Wrbti+roip(5)):min(size(im1,3),Wrbti+roip(6)))));
            
            
            
        end
        if ng>0
            if any(elt==6)
                [xgt,ygt,zgt,wgt]=GetGaussPointsWedge(ng,ns);
                Nt=[0.5*(1-xgt-ygt).*(1-zgt),0.5*(xgt).*(1-zgt),0.5*(ygt).*(1-zgt),...
                    0.5*(1-xgt-ygt).*(1+zgt),0.5*(xgt).*(1+zgt),0.5*(ygt).*(1+zgt)];
                Nt_r=[-0.5*(1-zgt),0.5*(1-zgt),(0*ygt),...
                    -0.5*(1+zgt),0.5*(1+zgt),(0*ygt)];
                Nt_s=[-0.5*(1-zgt),(0*xgt),0.5*(1-zgt),...
                    -0.5*(1-zgt),(0*xgt),0.5*(1-zgt)];
                Nt_t=[-0.5*(1-xgt-ygt),-0.5*xgt,-0.5*ygt,...
                    0.5*(1-xgt-ygt),0.5*xgt,0.5*ygt];
            end
            if any(elt==8)
                [xgq,ygq,zgq,wgq]=GetGaussPointsHexaedron(ng,ns);
                Nq=[0.125*(1-xgq).*(1-ygq).*(1-zgq),0.125*(1+xgq).*(1-ygq).*(1-zgq),0.125*(1+xgq).*(1+ygq).*(1-zgq),0.125*(1-xgq).*(1+ygq).*(1-zgq),...
                    0.125*(1-xgq).*(1-ygq).*(1+zgq),0.125*(1+xgq).*(1-ygq).*(1+zgq),0.125*(1+xgq).*(1+ygq).*(1+zgq),0.125*(1-xgq).*(1+ygq).*(1+zgq)];
                Nq_r=[-0.125*(1-ygq).*(1-zgq),0.125*(1-ygq).*(1-zgq),0.125*(1+ygq).*(1-zgq),-0.125*(1+ygq).*(1-zgq),...
                    -0.125*(1-ygq).*(1+zgq),0.125*(1-ygq).*(1+zgq),0.125*(1+ygq).*(1+zgq),-0.125*(1+ygq).*(1+zgq)];
                Nq_s=[-0.125*(1-xgq).*(1-zgq),-0.125*(1+xgq).*(1-zgq),0.125*(1+xgq).*(1-zgq),0.125*(1-xgq).*(1-zgq),...
                    -0.125*(1-xgq).*(1+zgq),-0.125*(1+xgq).*(1+zgq),0.125*(1+xgq).*(1+zgq),0.125*(1-xgq).*(1+zgq)];
                Nq_t=[-0.125*(1-xgq).*(1-ygq),-0.125*(1+xgq).*(1-ygq),-0.125*(1+xgq).*(1+ygq),-0.125*(1-xgq).*(1+ygq),...
                    0.125*(1-xgq).*(1-ygq),0.125*(1+xgq).*(1-ygq),0.125*(1+xgq).*(1+ygq),0.125*(1-xgq).*(1+ygq)];
            end
            
        end
        for i1=1:prod(Nelems)
            %                        display(sprintf('elemen %d /%d\n',i1,prod(Nelems)));
            inods=conn(i1,1:elt(i1));
%             if mfilter
%                 Un=zeros(numel(inods),1);
%                 Vn=zeros(numel(inods),1);
%                 Wn=zeros(numel(inods),1);
%                 for in=1:elt(i1)
%                     ee=sum(conn==inods(in),2)>0;
%                     ee=conn(ee,:);
%                     ee=unique(ee(:));
%                     Un(in)=median(Ui(ee));
%                     Vn(in)=median(Ui(ee+1*Nn));
%                     Wn(in)=median(Ui(ee+2*Nn));
%                 end
%             else
                Un=Ui(inods+0*Nn);
                Vn=Ui(inods+1*Nn);
                Wn=Ui(inods+2*Nn);
%            end
            xn=xo(inods);
            yn=yo(inods);
            zn=zo(inods);
            inde=[];
            if ng==0||any(face_elts==i1)
                ipix=max(1,ceil(min(xn)-2)):min(sizeim(1),floor(max(xn)+2));
                jpix=max(1,ceil(min(yn)-2)):min(sizeim(2),floor(max(yn)+2));
                kpix=max(1,ceil(min(zn)-2)):min(sizeim(3),floor(max(zn)+2));
                im0e=double(im0(ipix,jpix,kpix));
                if pstep>1
                    %                     ipixo=int32(1:pstep:size(im0e,1))';
                    %                     jpixo=int32(1:pstep:size(im0e,2))';
                    %                     kpixo=int32(1:pstep:size(im0e,3))';
                    %                     [ypix,xpix,zpix]=meshgrid(jpix(jpixo),ipix(ipixo),kpix(kpixo));
                    [ypix,xpix,zpix]=meshgrid(jpix(1:pstep:end),ipix(1:pstep:end),kpix(1:pstep:end));
                else
                    
                    [ypix,xpix,zpix]=meshgrid(jpix,ipix,kpix);
                end
                if invmap||~(elt(i1)==8)
                    [xg,yg,zg,wg]=GetGaussPointsVoxels(elt(i1),xn,yn,zn,xpix(:),ypix(:),zpix(:));
                else
                    xg=-1+2*(xpix(:)-min(xn))/(max(xn)-min(xn));
                    yg=-1+2*(ypix(:)-min(yn))/(max(yn)-min(yn));
                    zg=-1+2*(zpix(:)-min(zn))/(max(zn)-min(zn));
                    switch elt(i1)
                        case 6
                            wg=~(xg<0|yg<0|1-xg-yg<0|abs(zg)>1);
                        case 4
                            wg=~(xg<0|yg<0|zg<0|1-xg-yg-zg<0);
                        case 8
                            wg=~(abs(xg)>1|abs(yg)>1|abs(zg)>1);
                    end
                end
                N=GetFiniteElementShapeFunctions3D(elt(i1),xg(:),yg(:),zg(:));
                if elt(i1)<8
                    N=N(:,1:elt(i1));
                end
                
                if any(face_elts==i1)
                    heaviside=double(crack(ipix+roi(1)-zone(1),jpix+roi(3)-zone(3),kpix+roi(5)-zone(5))>=0);
                    
                    for in=1:elt(i1)
                        ienr=find(face_nodes==inods(in));
                        if ~isempty(ienr)
                            hn=double(crackn(ienr)>=0);
                            N=[N,N(:,in).*(heaviside(:)-hn)];
                            inde=[inde,ienr];
                            Un=[Un;Ui(3*Nn+ienr+0*nenr)];
                            Vn=[Vn;Ui(3*Nn+ienr+1*nenr)];
                            Wn=[Wn;Ui(3*Nn+ienr+2*nenr)];
                        end
                    end
                end
                phidf=mexPhidf3D(im0e,N,pstep);
                if pstep>1
                    im0e=im0e(1:pstep:end,1:pstep:end,1:pstep:end);
                end
            else
                if elt(i1)==6
                    N=Nt;wg=wgt;
                    N_r=Nt_r;N_s=Nt_s;N_t=Nt_t;
                elseif elt(i1)==8
                    N=Nq;wg=wgq;
                    N_r=Nq_r;N_s=Nq_s;N_t=Nq_t;
                end
                xpix=N*xn;
                ypix=N*yn;
                zpix=N*zn;
                dxdr=N_r*xn;
                dydr=N_r*yn;
                dzdr=N_r*zn;
                dxds=N_s*xn;
                dyds=N_s*yn;
                dzds=N_s*zn;
                dxdt=N_t*xn;
                dydt=N_t*yn;
                dzdt=N_t*zn;
                detJ=dxdr.*dyds.*dzdt+dxds.*dydt.*dzdr+dxdt.*dydr.*dzds...
                    -dzdr.*dyds.*dxdt-dzds.*dydt.*dxdr-dzdt.*dydr.*dxds;
                wg=wg.*detJ;
                
                [gx,gy,gz]=mexGradLinear3D(xpix(:),ypix(:),zpix(:),im0);
                gx=diag(sparse(gx(:)));
                gy=diag(sparse(gy(:)));
                gz=diag(sparse(gz(:)));
                phidf=[gx*N,gy*N,gz*N];
                im0e=mexInterpLinear3D(xpix(:),ypix(:),zpix(:),im0);
            end
            xpix=xpix(:)+(N*Un)/pscale+roip(1)-1;
            ypix=ypix(:)+(N*Vn)/pscale+roip(3)-1;
            zpix=zpix(:)+(N*Wn)/pscale+roip(5)-1;
            xmin=max(floor(min(xpix))-2,1);xmax=min(ceil(max(xpix))+2,size(im1,1));
            ymin=max(floor(min(ypix))-2,1);ymax=min(ceil(max(ypix))+2,size(im1,2));
            zmin=max(floor(min(zpix))-2,1);zmax=min(ceil(max(zpix))+2,size(im1,3));
            im1el=double(im1(xmin:xmax,ymin:ymax,zmin:zmax));
            im1e=mexInterpLinear3D(xpix-xmin+1,ypix-ymin+1,zpix-zmin+1,im1el);
            %                         im1e=mexInterpLinear3D(xpix(:)+(N*Un)/pscale+roip(1)-1,...
            %                 ypix(:)+(N*Vn)/pscale+roip(3)-1,...
            %                 zpix(:)+(N*Wn)/pscale+roip(5)-1,...
            %                 im1);
            if inorm
                im0e=im0e-mean(im0e(:));
                im1e=im1e-mean(im1e(:));
                sc=max(1,std(im0e(:)))/max(1,std(im1e(:)));
                im1e=sc*im1e;
            else
                 im0e=im0e-mean0;
                im1e=im1e-mean1;
                sc=max(1,std0)/max(1,std1);
                im1e=sc*im1e;
            end
            im1e=im0e(:)-im1e(:);
            wg=diag(sparse(wg(:)));
            indo=[inods,3*Nn+inde+0*nenr,inods+Nn,3*Nn+inde+1*nenr,inods+2*Nn,3*Nn+inde+2*nenr];
            if do_matrix
                [indj,indi]=meshgrid(indo,indo);
                ne=1:(3*size(N,2))^2;
                Mel=(wg*phidf);
                Mel=phidf'*Mel;
                indig(nel+ne)=indi(:);
                indjg(nel+ne)=indj(:);
                val(nel+ne)=Mel(:);
                nel=nel+numel(Mel);
            end
            
            Fel=wg*(im1e(:));
            Fel=phidf'*Fel;
            %            Fels(i1,[1:elt(i1),size(conn,2)+(1:elt(i1)),2*size(conn,2)+(1:elt(i1))])=Fel;
            F(indo(:))=F(indo(:))+Fel(:);
            mdisce=sum(wg*abs(im1e(:)))/sum(diag(wg));
            mdisc=mdisc+mdisce;
            disc(i1)=mdisce;
        end
        mdisc=mdisc/length(elt);
        if do_matrix
            if nel<length(indig)
                indig(nel+1:length(indig))=[];
                indjg(nel+1:length(indjg))=[];
                val(nel+1:length(val))=[];
            end
            M=sparse(indig,indjg,val,3*(Nn+nenr),3*(Nn+nenr));
        end
        %        F=full(sparse([conn,Nn+conn,2*Nn+conn],ones(size(Fels)),Fels,3*Nn,1));
    end
end
