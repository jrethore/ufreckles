function [U,PP]=Solve2D(Uini,iscale,nmod,LL,relts)
clear MedianFilter
if nargout>1, Fu=0*Uini;end
if nargin<4
    cond=0;
else
    if ~isempty(LL)
        cond=1;
    else
        cond=0;
    end
end
if nargin<5,relts=[];end
pscale=2^(iscale-1);
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if strcmp(param.basis,'nurbs')&&iscale==1
    [U]=SolveNURBS2D(Uini,iscale,nmod);
    return
end
pstep=1;
if isfield(param0,'psample')
    pstep=param0.psample;
end
pgd=0;
if isfield(param0,'do_pgd_prediction')
    pgd=param0.do_pgd_prediction;
end
inorm=false;
if isfield(param0,'normalize_grey_level')
    inorm=param0.normalize_grey_level;
end
mfilter=0;
if isfield(param0,'regularization_type')
    if strcmp(param0.regularization_type,'median')
        param0.regularization_type='none';
        lmed=param0.regularization_parameter;
        mfilter=1;
    elseif strcmp(param0.regularization_type,'tiko')
        mfilter=1;lmed=1;
    end
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

mesh_file=fullfile('TMP',sprintf('%d_mesh_%d.mat',nmod,iscale-1));
load(mesh_file,'rflag','xo','yo','Nnodes','Nelems','Smesh','conn','elt','ng','ns','selected');
invmap=~rflag;
%invmap=isfield(param,'mesh_file')&&(iscale==1);
Fc=zeros(2*prod(Nnodes),1);


load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim');
dynamic=max(im0(:))-min(im0(:));
[mean0,std0]=mexImScalling2D(im0);
mean1=0;std1=1;
M=0;
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
    V=repmat(cos(2*pi*(xo)*ki).*sin(2*pi*(yo)*ki),2,1);
    hh=min(abs(diff(xo(conn(:,1:2)),[],2)+1i*diff(yo(conn(:,1:2)),[],2)));
    if isempty(hh)
        hh=mean(param.mesh_size);
    end
    if param0.regularization_parameter>5*hh
        mfilter=0;
    end
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
                    if pgd
                        Ui=U(:,ijm-1,iz);
                        
                        display('Automatic prediction...');
                        L=U(:,max(1,ijm-4):(ijm-1),iz);
                        dL=Inf;
                        ii=1;
                        while norm(dL)>0.001&&ii<100
                            [merror,disc]=Assemble(false);
                            MM=L'*(M-R-P)*L;
                            FF=L'*Fc;
                            dL=MM\FF;
                            Ui=Ui+L*dL;
                            ii=ii+1;
                        end
                        ii=1;
                    else
                        Ui=U(:,ijm-1,iz)+DU;
                        
                    end
                end
                
            end
            
        end
        merrorp=Inf;
        while ( res>conv && ii< maxiter)
            [merror,disc]=Assemble((ii==1)&&(ijm==1)&&(icamd==1)&&(maxiter>2));
            if (ii==1)&&(ijm==1)&&(icamd==1)
                switch reg_type
                    case 'tiko'
                        R=AssembleRegularizationOperator([],iscale,nmod,2);
                        a=(V'*M*V)/(V'*R*V);
                        a=a*(2*lm/lc)^2;
                        R=a*R;
                        P=sparse(size(M,1),size(M,2));
                    case 'equilibrium_gap'
                        CreateGradBasisFunction(iscale,nmod);
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
                    M1=[];
                    M2=[];
                    if maxiter>2
                        try
                            M1 = ichol(M);
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
                end
            end
            
            F=Fc-R*Ui;
            if cond
                F=LL'*(F+M*Ui);
                dP=MM\F;
                dU=LL*dP-Ui;
                PP(:,ijm,iz)=dP;
            else
                if maxiter>2
                    if issparse(M)
                        dU=pcg(M,F+M*Ui,1e-6,1000,M1,M2,Ui);
                        dU=dU-Ui;
                    else
                        dU=F./M;
                    end
                else
                    dU=0*Ui;
                end
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
                Ui=MedianFilter(Ui,mesh_file,lmed);
            end
            merrorp=merror;
        end
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
                im1=(readim(fildef));
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
            im1=double(im1);
            [mean1,std1]=mexImScalling2D(im1(roip(1):roip(2),roip(3):roip(4)));
            
            
        end
        if ng>0
            if numel(ns)==1
                nso=[1,1];
            else
                nso=ns;
            end
            if any(elt==3)
                [xgt,ygt,wgt]=GetGaussPointsTriangle(ng,nso);
                Nt=[0.5*(1-xgt-ygt),0.5*(xgt),0.5*(ygt)];
                Nt_r=[-0.5+0*xgt,0.5+0*xgt,0*ygt];
                Nt_s=[-0.5+0*ygt,0*xgt,0.5+0*ygt];
            end
            if any(elt==4)
                [xgq,ygq,wgq]=GetGaussPointsQuadrangle(ng,nso);
                Nq=[0.25*(1-xgq).*(1-ygq),0.25*(1+xgq).*(1-ygq),0.25*(1+xgq).*(1+ygq),0.25*(1-xgq).*(1+ygq)];
                Nq_r=[-0.25*(1-ygq),0.25*(1-ygq),0.25*(1+ygq),-0.25*(1+ygq)];
                Nq_s=[-0.25*(1-xgq),-0.25*(1+xgq),0.25*(1+xgq),0.25*(1-xgq)];
            end
        end
        Nn=prod(Nnodes);
        for i1=1:prod(Nelems)
            if  ~any(relts==i1)
                inods=conn(i1,1:elt(i1));
                Un=Ui(inods+0*Nn);
                Vn=Ui(inods+1*Nn);
                xn=xo(inods);
                yn=yo(inods);
                inde=[];
                if ng==0
                    ipix=max(1,ceil(min(xn)-2)):min(sizeim(1),floor(max(xn)+2));
                    jpix=max(1,ceil(min(yn)-2)):min(sizeim(2),floor(max(yn)+2));
                    im0e=double(im0(ipix,jpix));
                    if pstep>1
                        [ypix,xpix]=meshgrid(jpix(1:pstep:end),ipix(1:pstep:end));
                    else
                        [ypix,xpix]=meshgrid(jpix,ipix);
                    end
                    
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
                    phidf=mexPhidf2D(im0e,N,pstep);
                    if pstep>1
                        im0e=im0e(1:pstep:end,1:pstep:end);
                    end
                    
                else
                    if elt(i1)==3
                        if numel(ns)==1
                            
                            nsx=max(2,round(ns*max(abs(mean(Nt_s*xn)),abs(mean(Nt_r*xn)))));
                            nsy=max(2,round(ns*max(abs(mean(Nt_s*yn)),abs(mean(Nt_r*yn)))));
                            [xgs,ygs,wg]=GetGaussPointsTriangle(ng,[nsx,nsy],xn,yn);
                            
                            N=[0.5*(1-xgs-ygs),0.5*(xgs),0.5*(ygs)];
                            N_r=[-0.5+0*xgs,0.5+0*xgs,0*ygs];
                            N_s=[-0.5+0*ygs,0*xgs,0.5+0*ygs];
                        else
                            N=Nt;wg=wgt;
                            N_r=Nt_r;N_s=Nt_s;
                        end
                    elseif elt(i1)==4
                        if numel(ns)==1
                            nsx=max(2,round(ns*max(abs(mean(Nq_r*xn)),abs(mean(Nq_s*xn)))));
                            nsy=max(2,round(ns*max(abs(mean(Nq_r*yn)),abs(mean(Nq_s*yn)))));
                            [xgs,ygs,wg]=GetGaussPointsQuadrangle(ng,[nsx,nsy],xn,yn);
                            N=[0.25*(1-xgs).*(1-ygs),0.25*(1+xgs).*(1-ygs),0.25*(1+xgs).*(1+ygs),0.25*(1-xgs).*(1+ygs)];
                            N_r=[-0.25*(1-ygs),0.25*(1-ygs),0.25*(1+ygs),-0.25*(1+ygs)];
                            N_s=[-0.25*(1-xgs),-0.25*(1+xgs),0.25*(1+xgs),0.25*(1-xgs)];
                        else
                            N=Nq;wg=wgq;
                            N_r=Nq_r;N_s=Nq_s;
                        end
                    end
                    xpix=N*xn;
                    ypix=N*yn;
                    dxdr=N_r*xn;
                    dydr=N_r*yn;
                    dxds=N_s*xn;
                    dyds=N_s*yn;
                    detJ=(dxdr.*dyds-dydr.*dxds);
                    wg=wg.*detJ;
                    ipix=max(1,ceil(min(xn)-2)):min(sizeim(1),floor(max(xn)+2));
                    jpix=max(1,ceil(min(yn)-2)):min(sizeim(2),floor(max(yn)+2));
                    im0e=im0(ipix,jpix);
                    gradx=mexFDGradient(im0e);
                    grady=mexFDGradient(im0e');
                    grady=grady';
                    gx=mexInterpLinear(xpix(:)-ipix(1)+1,ypix(:)-jpix(1)+1,gradx);
                    gy=mexInterpLinear(xpix(:)-ipix(1)+1,ypix(:)-jpix(1)+1,grady);
                    gx=diag(sparse(gx(:)));
                    gy=diag(sparse(gy(:)));
                    phidf=[gx*N,gy*N];
                    im0ei=mexInterpSpline(xpix(:)-ipix(1)+1,ypix(:)-jpix(1)+1,im0e);
                    im0e=im0ei;
                end
                
                xpix=xpix(:)+(N*Un)/pscale+roip(1)-1;
                ypix=ypix(:)+(N*Vn)/pscale+roip(3)-1;
                xmin=max(floor(min(xpix))-2,1);xmax=min(ceil(max(xpix))+2,size(im1,1));
                ymin=max(floor(min(ypix))-2,1);ymax=min(ceil(max(ypix))+2,size(im1,2));
                im1el=double(im1(xmin:xmax,ymin:ymax));
                
                if numel(im1el)==0
                    im1e=-ones(size(im0e));
                else
                    if iscale==1
                        im1e=mexInterpSpline(xpix-xmin+1,ypix-ymin+1,im1el);
                    else
                        im1e=mexInterpLinear(xpix-xmin+1,ypix-ymin+1,im1el);
                    end
                end
                %                 if iscale==1
                %                     im1e=mexInterpSpline(xpix(:)+(N*Un)/pscale+roip(1)-1,ypix(:)+(N*Vn)/pscale+roip(3)-1,im1);
                %                 else
                %                     im1e=mexInterpLinear(xpix(:)+(N*Un)/pscale+roip(1)-1,ypix(:)+(N*Vn)/pscale+roip(3)-1,im1);
                %                 end
                maske=(im1e<0)|abs(im1e)>100*mean0|isnan(im1e);
                im1e(maske)=im0e(maske);
                if any(wg>0)
                    wg=diag(sparse(wg(:)));
                    
                    mdisce=sum(wg*abs(im0e(:)-im1e(:)))/sum(diag(wg));
                    mdisc=mdisc+mdisce;
                    if inorm
                        if any(maske)
                            im1e=im1e-mean(im0e(:));
                            sc=1;
                        else
                            im1e=im1e-mean(im1e(:));
                            sc=max(1,std(im0e(:)))/max(1,std(im1e(:)));
                        end
                        im0e=im0e-mean(im0e(:));
                        im1e=sc*im1e;
                    else
                        
                        im0e=im0e-mean0;
                        im1e=im1e-mean1;
                        sc=max(1,std0)/max(1,std1);
                        im1e=sc*im1e;
                    end
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
