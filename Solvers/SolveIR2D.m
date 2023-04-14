function [U,PP]=SolveIR2D(Uini,iscale,nmod,LL,relts)
if nargout>1, Fu=0*Uini;end
if nargin<4,cond=0;else cond=1;end
if nargin<5,relts=[];end
pscale=2^(iscale-1);
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
to=param0.ir_calibration_data.to;
tx=param0.ir_calibration_data.tx;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if strcmp(param.basis,'nurbs')&&iscale==1
    error('NOT CODED YET')
end
mfilter=0;
if isfield(param0,'regularization_type')
    if strcmp(param0.regularization_type,'median')
        mfilter=1;
        clear MedianFilter
    end
end
pgd=0;
if isfield(param0,'do_pgd_prediction')
    pgd=param0.do_pgd_prediction;
end
dotopo=(~isfield(param,'topography'))&&isfield(param0,'calibration_data');
incremental=0;
if isfield(param0,'incremental_resolution')
    incremental=param0.incremental_resolution;
    if iscale>1
        incremental=0;
    end
end
roi=param0.roi;
if isfield(param0,'restart')
    restart=param0.restart;
else
    restart=1;
end
if iscell(param0.reference_image)
    error('NOT CODED YET');
    ncamr=length(param0.reference_image);
    docrosscorrelation=0;
    if isfield(param0,'cross_correlation')
        docrosscorrelation=param0.cross_correlation;
    end
else
    ncamr=1;
    docrosscorrelation=0;
end
if ~isfield(param0,'deformed_image')
    reader=VideoReader(param0.reference_image);
    nbf=reader.NumberOfFrames-1;
    if isfield(param0,'number_of_frames')
        nbf=param0.number_of_frames;
    end
    dim=1;
    if isfield(param0,'video_sampling')
        dim=param0.video_sampling;
    end
    frames=2:dim:nbf;
    nim=length(frames);
    ncamd=1;
else
    if iscell(param0.deformed_image)
        nim=size(param0.deformed_image,2);
        ncamd=size(param0.deformed_image,1);
    else
        nim=1;
        ncamd=1;
    end
end
if ncamr==1
    indcam=ncamd:-1:1;
else
    if docrosscorrelation
        indo=(ncamd:-1:1)';
        indcam=indo';
        for ic=1:ncamd-1
            indo=circshift(indo,1);
            indcam=[indcam;indo'];
        end
    else
        indcam=(1:ncamd)';
    end
end

if restart
    nmax=nim;
else
    if iscale==1
        nmax=nim;
    else
        if nim<2+dotopo
            nmax=nim;
        else
            nmax=2+dotopo;
        end
    end
end

reg_type='none';
if isfield(param0,'regularization_type')
    reg_type=param0.regularization_type;
    if iscale>1&&~(strcmp(reg_type,'none'))
        reg_type='tiko';
    end
end

maxiter=param0.iter_max;
conv=param0.convergance_limit;
%invmap=isfield(param,'mesh_file')&&(iscale==1);
disp(sprintf('Starting resolution for scale %d...',iscale));

mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
load(mesh_file,'rflag','xo','yo','Nnodes','Nelems','Smesh','conn','elt','ng','ns','selected');
if ng==0
else
    error('NOT CODED')
end
invmap=~rflag;
%invmap=isfield(param,'mesh_file')&&(iscale==1);
Fc=zeros(2*prod(Nnodes),1);
Ft=zeros(prod(Nnodes),1);

CreateGradBasisFunction(iscale,nmod);
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim');
dynamic=max(im0(:))-min(im0(:));
M=0;
T=0;
im1=0;
U=Uini;
if cond
    PP=zeros([size(LL,2),size(U,2),size(U,3)]);
end
if ~strcmp(reg_type,'none')
    lc=(min(sizeim)/10);
    lm=lc;
    if isfield(param0,'regularization_parameter')
        lm=param0.regularization_parameter/pscale;
        if isfield(param0,'detect')&&iscale>1
            if param0.detect
                lm=0.1*lm;
            end
        end
        
    end
    ki=1/lc;
    V=repmat(cos(2*pi*(xo)*ki).*cos(2*pi*(yo)*ki),2,1);
end

fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
fwrite(fid,1);
if dynamic>255
    fwrite(fid,255);
else
    fwrite(fid,dynamic);
end
icamr=1;
for ijm=1:nmax
    display(sprintf('Image %d/%d',ijm,nmax));
    for icamd=1:size(indcam,2)
        iz=indcam(icamr,icamd)+(icamr-1)*size(indcam,2)*docrosscorrelation;
        res=1;
        ii=1;
        
        
        if restart
            Ui=Uini(:,ijm,iz);
            if ijm==1&&dotopo&&(indcam(icamr,icamd)==icamr)
                Ui=0*Ui;
                res=0;
                disc=zeros(numel(im0),1);
            end
        else
            if ijm==1
                Ui=Uini(:,ijm,iz);
                if dotopo&&(indcam(icamr,icamd)==icamr)
                    Ui=0*Ui;
                    res=0;
                    disc=zeros(numel(im0),1);
                end
            elseif ijm<3+dotopo
                Ui=Uini(:,ijm,iz);
            else
                DU=U(:,ijm-1,iz)-U(:,ijm-2,iz);
                if incremental
                    Ui=DU;
                else
                    Ui=U(:,ijm-1,iz)+DU;
                    if pgd
                        display('Automatic prediction...');
                        L=U(:,max(1,ijm-4):(ijm-1),iz);
                        dL=Inf;
                        ii=1;
                        while norm(dL)>0.001&&ii<10
                            [merror,disc]=Assemble(false);
                            MM=L'*(M-R-P)*L;
                            FF=L'*Fc;
                            dL=MM\FF;
                            Ui=Ui+L*dL;
                            ii=ii+1;
                        end
                        ii=1;
                    end
                end
                
            end
            
        end
        Ti=Ui(2*prod(Nnodes)+(1:prod(Nnodes)));
        Ui=Ui((1:2*prod(Nnodes)));
        merrorp=Inf;
        while ( res>conv && ii< maxiter)
            AssembleT((ii==1)&&(ijm==1)&&(icamd==1));
            Ti=T\Ft;
            if mfilter
                Ti=MedianFilter(Ti,mesh_file,lm);
            end
            [merror,disc]=Assemble(1);
            
            if (ii==1)&&(ijm==1)&&(icamd==1)
                switch reg_type
                    case 'tiko'
                        R=AssembleRegularizationOperator([],iscale,nmod,2);
                        a=(V'*M*V)/(V'*R*V);
                        a=a*(2*lm/lc)^2;
                        R=a*R;
                        P=sparse(size(M,1),size(M,2));
                    case 'equilibrium_gap'
                        selected_nodes=find(~selected(:));
                        Rb=AssembleRegularizationOperator(selected_nodes,iscale,nmod,2);
                        D=diag(abs(diag(Rb))>0);
                        a=(V'*D*M*V)/(V'*Rb*V);
                        a=a*(0.1*2*lm/lc)^2;
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
                    otherwise
                        R=sparse(size(M,1),size(M,2));
                        P=sparse(size(M,1),size(M,2));
                end
                save(fullfile('TMP',sprintf('%d_operator_%d',nmod,iscale-1)),'M','R');
                
                M=M+R+P;
                if cond
                    MM=LL'*(M*LL);
                else
                    try
                        M1 = ichol(M);
                        M2 = M1';
                    catch err
                        disp([err.identifier ' : ' err.message ' => using diagonal compensation'])
                        diagcomp=(max(sum(abs(M),2)./diag(M))-2)/100;
                        try
                            M1 = ichol(M,struct('type','ict','droptol',1e-3,'diagcomp',diagcomp));
                            M2 = M1';
                        catch err
                            disp([err.identifier ' : ' err.message ' => no compensation'])
                            ido=isnan(diagcomp);
                            M(ido,ido)=1;
                            M1=[];
                            M2=[];
                        end
                    end
                end
            end
            
            F=Fc-R*Ui+M*Ui;
            if cond
                F=LL'*F;
                dP=MM\F;
                dU=LL*dP-Ui;
                PP(:,ijm,iz)=dP;
            else
                dU=pcg(M,F,1e-6,1000,M1,M2,Ui);
                dU=dU-Ui;
            end
            disp(sprintf('At iteration # %d',ii));
            disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100/dynamic));
            if ii==1
                nU0=norm(Ui+dU(:));
            else
                res=norm(dU(:))/nU0;
                disp(sprintf('|dU|=%f',res));
                if mfilter
                    res=min(res,abs(merrorp-merror)/abs(merror));
                end
                
            end
            Ui=Ui+dU;
            ii=ii+1;
            if mfilter
                Ui=MedianFilter(Ui,mesh_file,lm);
            end
            merrorp=merror;
        end
        Ui=[Ui;Ti];
        if incremental
            UpdateMesh(Ui,iscale,nmod,ijm);
            if ijm>1
                U(:,ijm,iz)=U(:,ijm-1,iz)+Ui;
            else
                U(:,ijm,iz)=Ui;
            end
            load(mesh_file,'xo','yo');
        else
            U(:,ijm,iz)=Ui;
        end
        if nargout>1
            Fu(:,ijm)=Fc;
        end
        if dynamic>255
            disc=255*disc/dynamic;
        end
        fwrite(fid,abs(disc));
    end
end
disp(sprintf('Enlapsed time for resolution = %6.2f s',cputime -ttic));
fclose(fid);
    function [mdisc,disc]=Assemble(do_matrix)
        mdisc=0;
        Fc=0*Fc;
        roip=(roi-1)/pscale+1;
        if do_matrix
            nind=(2*3)^2*sum(elt==3)+(2*4)^2*sum(elt==4);
            indig=zeros(nind,1);
            indjg=zeros(nind,1);
            val=zeros(nind,1);
            disc=zeros(length(elt),1);
            nel=0;
        end
        if ii==1
            if ~isfield(param0,'deformed_image')
                im1=double(readim(reader,frames(ijm)));
            else
                if (nim==1)&&(ncamd==1)
                    fildef=param0.deformed_image;
                else
                    fildef=param0.deformed_image{indcam(icamr,icamd),ijm};
                end
                im1=double(readim(fildef));
            end
            if length(size(im1))==3
                im1=mean(im1,3);
            end
            if iscale>1
                for ip=1:(iscale-1)
                    scale=2;
                    imsiz0=size(im1);
                    imsiz1=floor(imsiz0/2);
                    nn=2*imsiz1;
                    im1=im1(1:nn(1),1:nn(2));
                    
                    im1=reshape(im1,scale,prod(nn)/scale);
                    im1=mean(im1,1);
                    nn(1)=nn(1)/scale;
                    im1=reshape(im1,nn);
                    
                    im1=im1';
                    im1=reshape(im1,scale,prod(nn)/scale);
                    im1=mean(im1,1);
                    nn(2)=nn(2)/scale;
                    im1=reshape(im1,nn([2,1]));
                    im1=im1';
                    
                end
            end
        end
        Nn=prod(Nnodes);
        for i1=1:prod(Nelems)
            if  ~any(relts==i1)
                inods=conn(i1,1:elt(i1));
                Un=Ui(inods+0*Nn);
                Vn=Ui(inods+1*Nn);
                Tn=Ti(inods);
                xn=xo(inods);
                yn=yo(inods);
                inde=[];
                if ng==0
                    ipix=max(1,ceil(min(xn)-2)):min(sizeim(1),floor(max(xn)+2));
                    jpix=max(1,ceil(min(yn)-2)):min(sizeim(2),floor(max(yn)+2));
                    im0e=im0(ipix,jpix);
                    [ypix,xpix]=meshgrid(jpix,ipix);
                    
                    if invmap||elt(i1)==3
                        [xg,yg,wg]=GetGaussPointsPixels(elt(i1),xn,yn,xpix(:),ypix(:));
                    else
                        xg=-1+2*(xpix(:)-min(xn))/(max(xn)-min(xn));
                        yg=-1+2*(ypix(:)-min(yn))/(max(yn)-min(yn));
                        switch elt(i1)
                            case 3
                                wg=~(xg<0|yg<0|1-xg-yg<0);
                            case 4
                                wg=~(abs(xg)>1|abs(yg)>1);
                        end
                    end
                    if any(wg>0)
                        N=GetFiniteElementShapeFunctions(elt(i1),xg(:),yg(:));
                        if elt(i1)<4
                            N=N(:,1:elt(i1));
                        end
                        if iscale==1
                            im1e=mexInterpSpline(xpix(:)+(N*Un)/pscale+roip(1)-1,ypix(:)+(N*Vn)/pscale+roip(3)-1,im1);
                        else
                            im1e=mexInterpLinear(xpix(:)+(N*Un)/pscale+roip(1)-1,ypix(:)+(N*Vn)/pscale+roip(3)-1,im1);
                        end
                    maske=im1e<0;
                        phidf=mexPhidf2D(reshape(im1e,size(im0e)),N);
                        wg=diag(sparse(wg(:)));
                        
                        im0e=im0e(:)+(N*Tn).*(tx-im0e(:))/(tx-to);
                        mdisce=sum(wg*abs(im0e(:)-im1e(:)))/sum(diag(wg));
                        mdisc=mdisc+mdisce;
                        im1e=im0e(:)-im1e(:);
                        im1e(maske)=0;
                        indo=[inods,inods+Nn];
                        if do_matrix
                            Mel=(wg*phidf);
                            Mel=phidf'*Mel;
                            [indj,indi]=meshgrid(indo,indo);
                            ne=1:(2*size(N,2))^2;
                            indig(nel+ne)=indi(:);
                            indjg(nel+ne)=indj(:);
                            val(nel+ne)=Mel(:);
                            nel=nel+numel(Mel);
                        end
                        Fel=wg*(im1e(:));
                        Fel=phidf'*Fel;
                        Fc(indo(:))=Fc(indo(:))+Fel(:);
                        disc(i1)=mdisce;
                    end
                end
            end
            mdisc=mdisc/length(elt);
            if do_matrix
                if nel<length(indig)
                    indig(nel+1:length(indig))=[];
                    indjg(nel+1:length(indjg))=[];
                    val(nel+1:length(val))=[];
                end
                M=sparse(indig,indjg,val,2*Nn,2*Nn);
            end
        end
    end
    function AssembleT(do_matrix)
        Ft=0*Ft;
        roip=(roi-1)/pscale+1;
        if do_matrix
            nind=(3)^2*sum(elt==3)+(4)^2*sum(elt==4);
            indig=zeros(nind,1);
            indjg=zeros(nind,1);
            val=zeros(nind,1);
            disc=zeros(length(elt),1);
            nel=0;
        end
        if ii==1
            if ~isfield(param0,'deformed_image')
                im1=double(readim(reader,frames(ijm)));
            else
                if (nim==1)&&(ncamd==1)
                    fildef=param0.deformed_image;
                else
                    fildef=param0.deformed_image{indcam(icamr,icamd),ijm};
                end
                im1=double(readim(fildef));
            end
            if length(size(im1))==3
                im1=mean(im1,3);
            end
            if iscale>1
                for ip=1:(iscale-1)
                    scale=2;
                    imsiz0=size(im1);
                    imsiz1=floor(imsiz0/2);
                    nn=2*imsiz1;
                    im1=im1(1:nn(1),1:nn(2));
                    
                    im1=reshape(im1,scale,prod(nn)/scale);
                    im1=mean(im1,1);
                    nn(1)=nn(1)/scale;
                    im1=reshape(im1,nn);
                    
                    im1=im1';
                    im1=reshape(im1,scale,prod(nn)/scale);
                    im1=mean(im1,1);
                    nn(2)=nn(2)/scale;
                    im1=reshape(im1,nn([2,1]));
                    im1=im1';
                    
                end
            end
        end
        Nn=prod(Nnodes);
        for i1=1:prod(Nelems)
            if  ~any(relts==i1)
                inods=conn(i1,1:elt(i1));
                Un=Ui(inods+0*Nn);
                Vn=Ui(inods+1*Nn);
                Tn=Ti(inods);
                xn=xo(inods);
                yn=yo(inods);
                inde=[];
                ipix=max(1,ceil(min(xn)-2)):min(sizeim(1),floor(max(xn)+2));
                jpix=max(1,ceil(min(yn)-2)):min(sizeim(2),floor(max(yn)+2));
                im0e=im0(ipix,jpix);
                [ypix,xpix]=meshgrid(jpix,ipix);
                
                if invmap||elt(i1)==3
                    [xg,yg,wg]=GetGaussPointsPixels(elt(i1),xn,yn,xpix(:),ypix(:));
                else
                    xg=-1+2*(xpix(:)-min(xn))/(max(xn)-min(xn));
                    yg=-1+2*(ypix(:)-min(yn))/(max(yn)-min(yn));
                    switch elt(i1)
                        case 3
                            wg=~(xg<0|yg<0|1-xg-yg<0);
                        case 4
                            wg=~(abs(xg)>1|abs(yg)>1);
                    end
                end
                N=GetFiniteElementShapeFunctions(elt(i1),xg(:),yg(:));
                if elt(i1)<4
                    N=N(:,1:elt(i1));
                end
                im01=(tx-im0e(:))/(tx-to);
                phidf=diag(sparse(im01))*N;
                if iscale==1
                    im1e=mexInterpSpline(xpix(:)+(N*Un)/pscale+roip(1)-1,ypix(:)+(N*Vn)/pscale+roip(3)-1,im1);
                else
                    im1e=mexInterpLinear(xpix(:)+(N*Un)/pscale+roip(1)-1,ypix(:)+(N*Vn)/pscale+roip(3)-1,im1);
                end
                maske=im1e<0;
                if any(wg>0)
                    wg=diag(sparse(wg(:)));
                    im1e=(im1e(:)-im0e(:));
                    im1e(maske)=0;
                    if do_matrix
                        Mel=(wg*phidf);
                        Mel=phidf'*Mel;
                        [indj,indi]=meshgrid(inods,inods);
                        ne=1:(size(N,2))^2;
                        indig(nel+ne)=indi(:);
                        indjg(nel+ne)=indj(:);
                        val(nel+ne)=Mel(:);
                        nel=nel+numel(Mel);
                    end
                    Fel=wg*(im1e(:));
                    Fel=phidf'*Fel;
                    Ft(inods(:))=Ft(inods(:))+Fel(:);
                end
            end
        end
        if do_matrix
            if nel<length(indig)
                indig(nel+1:length(indig))=[];
                indjg(nel+1:length(indjg))=[];
                val(nel+1:length(val))=[];
            end
            T=sparse(indig,indjg,val,Nn,Nn);
        end
    end
end
