function [U]=SolveNURBS2D(Uini,iscale,nmod)
pscale=2^(iscale-1);
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
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
inorm=false;
if isfield(param0,'normalize_grey_level')
    inorm=param0.normalize_grey_level;
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
if iscell(param0.deformed_image)
    nim=size(param0.deformed_image,2);
    ncamd=size(param0.deformed_image,1);
else
    nim=1;
    ncamd=1;
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
invmap=isfield(param,'mesh_file')&&(iscale==1);
disp(sprintf('Starting resolution for scale %d...',iscale));

mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
load(mesh_file,'xo','yo','uo','vo','Nbsnodes','Nbselems','Px','Py','degree');
p=degree;
F=zeros(2*prod(Nbsnodes),1);
Fe=zeros(prod(Nbsnodes),1);

CreateGradBasisFunction(iscale,nmod);
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim');
dynamic=max(im0(:))-min(im0(:));
M=0;
Me=0;
im1=0;
U=Uini;
if ~strcmp(reg_type,'none')
    lc=(min(sizeim)/10);
    lm=lc;
    if isfield(param0,'regularization_parameter')
        lm=param0.regularization_parameter/pscale;
    end
    ki=1/lc;
    load(fullfile('TMP',sprintf('%d_phio_%d',nmod,iscale-1)),'phio');
    L=phio'*phio;
    b=phio'*(cos(2*pi*(xo)*ki).*cos(2*pi*(yo)*ki));
    V=L\b;
    V=repmat(V,2,1);
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
                disc=zeros(prod(Nbselems),1);
            end
        else
            if ijm==1
                Ui=Uini(:,ijm,iz);
                if dotopo&&(indcam(icamr,icamd)==icamr)
                    Ui=0*Ui;
                    res=0;
                    disc=zeros(prod(Nbselems),1);
                end
            elseif ijm<3+dotopo
                Ui=Uini(:,ijm,iz);
            else
                DU=U(:,ijm-1,iz)-U(:,ijm-2,iz);
                if incremental
                    Ui=DU;
                else
                    Ui=U(:,ijm-1,iz)+DU;
                end
                
            end
            
        end
        
        while ( res>conv && ii< maxiter)
            [merror,disc]=Assemble((ii==1)&&(ijm==1)&&(icamd==1));
            if (ii==1)&&(ijm==1)&&(icamd==1)
                switch reg_type
                    case 'tiko'
                        load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                        load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
                        load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uxy','Uyx');
                        R= epsxx'*wdetJ*epsxx+Uxy'*wdetJ*Uxy...
                            +epsyy'*wdetJ*epsyy+Uyx'*wdetJ*Uyx;
                        clear epsxx epsyy Uxy Uyx
                        a=(V'*M*V)/(V'*R*V);
                    a=a*(2*lm/lc)^2;
                        R=a*R;
                        P=sparse(size(M,1),size(M,2));
                    otherwise
                        R=sparse(size(M,1),size(M,2));
                        P=sparse(size(M,1),size(M,2));
                        
                end
                M=M+R+P;
                try
                    M1 = ichol(M);
                    M2 = M1';
                catch err
                    disp([err.identifier ' : ' err.message ' => using diagonal compensation'])
                    diagcomp=(max(sum(abs(M),2)./diag(M))-2)/100;
                    M1 = ichol(M,struct('type','ict','droptol',1e-3,'diagcomp',diagcomp));
                    M2 = M1';
                end
            end
            F=F-R*Ui+M*Ui;
            dU=pcg(M,F,1e-6,1000,M1,M2,Ui);
            dU=dU-Ui;
            
            disp(sprintf('At iteration # %d',ii));
            disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100/dynamic));
            if ii==1
                nU0=norm(Ui+dU(:));
            else
                res=norm(dU(:))/nU0;
                disp(sprintf('|dU|=%f',res));
                
            end
            Ui=Ui+dU;
            ii=ii+1;
        end
        if incremental
            UpdateMesh(Ui,iscale,nmod,ijm);
            if ijm>1
                U(:,ijm,iz)=U(:,ijm-1,iz)+Ui;
            else
                U(:,ijm,iz)=Ui;
            end
        else
            U(:,ijm,iz)=Ui;
        end
        GetError((ijm==1)&&(icamd==1));
        disc=pcg(Me,Fe,1e-6,1000);
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
        F=0*F;
        nel=0;
        if do_matrix
            nind=(2*4)^2*prod(Nbselems);
            indig=zeros(nind,1);
            indjg=zeros(nind,1);
            val=zeros(nind,1);
            disc=zeros(prod(Nbselems),1);
        end
        if ii==1
            if (nim==1)&&(ncamd==1)
                fildef=param0.deformed_image;
            else
                fildef=param0.deformed_image{indcam(icamr,icamd),ijm};
            end
            im1=double(readim(fildef));
            if length(size(im1))==3
                im1=mean(im1,3);
            end
        end
        Nn=prod(Nbsnodes);
        [indjo,indio]=meshgrid((0:p(2)),(0:p(1)));
        for ix=1:Nbselems(1)
            inods=ix+p(1)+[0,1];
            ui=uo(inods(:))';
            ipix=(min(ui)+0.5):(max(ui)-0.5);
            [fx]=NURBSBasisFunc(ix+p(1),p(1),ipix,uo);
            indi=indio+ix;
            for iy=1:Nbselems(2)
                jnods=iy+p(2)+[0,1];
                vi=vo(jnods(:))';
                jpix=(min(vi)+0.5):(max(vi)-0.5);
                [fy]=NURBSBasisFunc(iy+p(2),p(2),jpix,vo);
                indj=indjo+iy;
                [Jpix,Ipix]=meshgrid(1:length(jpix),1:length(ipix));
                indn=sub2ind(Nbsnodes,indi(:),indj(:));
                N=fx(Ipix(:),1+indio(:)).*fy(Jpix(:),1+indjo(:));
                xpix=N*Px(indn);
                ypix=N*Py(indn);
                
                igx=min(max(1,(floor(min(xpix(:)))-1):(ceil(max(xpix(:)))+1)),sizeim(1));
                igy=min(max(1,(floor(min(ypix(:)))-1):(ceil(max(ypix(:)))+1)),sizeim(2));
                im0e=im0(igx,igy);
                [gyi,gxi]=gradient(im0e);
                gx=mexInterpLinear(xpix(:)-min(igx)+1,ypix(:)-min(igy)+1,gxi);
                gy=mexInterpLinear(xpix(:)-min(igx)+1,ypix(:)-min(igy)+1,gyi);
                
                im0e=mexInterpSpline(xpix(:),ypix(:),im0);
                %                [gx,gy]=mexGradLinear(xpix(:),ypix(:),im0);
                gx=diag(sparse(gx(:)));
                gy=diag(sparse(gy(:)));
                phidf=[gx*N,gy*N];
                Un=N*Ui(indn+0*Nn);
                Vn=N*Ui(indn+1*Nn);
                im1e=mexInterpSpline(xpix(:)+Un+roi(1)-1,...
                    ypix(:)+Vn+roi(3)-1,...
                    im1);
                    if inorm
                        im0e=im0e-mean(im0e(:));
                        im1e=im1e-mean(im1e(:));
                        sc=max(1,std(im0e(:)))/max(1,std(im1e(:)));
                        im1e=sc*im1e;
                    end
                im1e=im0e(:)-im1e(:);
                indno=[indn,indn+Nn];
                if do_matrix
                    Mel=phidf'*phidf;
                    [indjj,indii]=meshgrid(indno,indno);
                    ne=1:(2*size(N,2))^2;
                    indig(nel+ne)=indii(:);
                    indjg(nel+ne)=indjj(:);
                    val(nel+ne)=Mel(:);
                    nel=nel+numel(Mel);
                end
                Fel=(phidf'*im1e);
                F(indno(:))=F(indno(:))+Fel(:);
                mdisce=mean(abs(im1e(:)));
                mdisc=mdisc+mdisce;
                i1=sub2ind(Nbselems,ix,iy);
                disc(i1)=mdisce;
            end
        end
        mdisc=mdisc/prod(Nbselems);
        if do_matrix
            if nel<length(indig)
                indig(nel+1:length(indig))=[];
                indjg(nel+1:length(indjg))=[];
                val(nel+1:length(val))=[];
            end
            M=sparse(indig,indjg,val,2*(Nn),2*(Nn));
            clear val indig indjg
        end
    end
    function [disc]=GetError(do_matrix)
        Fe=0*Fe;
        nel=0;
        if do_matrix
            nind=(4)^2*prod(Nbselems);
            indig=zeros(nind,1);
            indjg=zeros(nind,1);
            val=zeros(nind,1);
            disc=zeros(prod(Nbselems),1);
        end
        if (nim==1)&&(ncamd==1)
            fildef=param0.deformed_image;
        else
            fildef=param0.deformed_image{indcam(icamr,icamd),ijm};
        end
        im1=double(readim(fildef));
        if length(size(im1))==3
            im1=mean(im1,3);
        end
        Nn=prod(Nbsnodes);
        [indjo,indio]=meshgrid((0:p(2)),(0:p(1)));
        for ix=1:Nbselems(1)
            inods=ix+p(1)+[0,1];
            ui=uo(inods(:))';
            ipix=(min(ui)+0.5):(max(ui)-0.5);
            [fx]=NURBSBasisFunc(ix+p(1),p(1),ipix,uo);
            indi=indio+ix;
            for iy=1:Nbselems(2)
                jnods=iy+p(2)+[0,1];
                vi=vo(jnods(:))';
                jpix=(min(vi)+0.5):(max(vi)-0.5);
                [fy]=NURBSBasisFunc(iy+p(2),p(2),jpix,vo);
                indj=indjo+iy;
                [Jpix,Ipix]=meshgrid(1:length(jpix),1:length(ipix));
                indn=sub2ind(Nbsnodes,indi(:),indj(:));
                N=fx(Ipix(:),1+indio(:)).*fy(Jpix(:),1+indjo(:));
                xpix=N*Px(indn);
                ypix=N*Py(indn);
                im0e=mexInterpSpline(xpix(:),ypix(:),im0);
                Un=Ui(indn+0*Nn);mUn=N*Un;
                Vn=Ui(indn+1*Nn);mVn=N*Vn;
                im1e=mexInterpSpline(xpix(:)+mUn+roi(1)-1,...
                    ypix(:)+mVn+roi(3)-1,...
                    im1);
                im1e=im0e(:)-im1e(:);
                indno=[indn];
                if do_matrix
                    Mel=N'*N;
                    [indjj,indii]=meshgrid(indno,indno);
                    ne=1:(size(N,2))^2;
                    indig(nel+ne)=indii(:);
                    indjg(nel+ne)=indjj(:);
                    val(nel+ne)=Mel(:);
                    nel=nel+numel(Mel);
                end
                Fel=(N'*im1e);
                Fe(indno(:))=Fe(indno(:))+Fel(:);
            end
        end
        if do_matrix
            if nel<length(indig)
                indig(nel+1:length(indig))=[];
                indjg(nel+1:length(indjg))=[];
                val(nel+1:length(val))=[];
            end
            Me=sparse(indig,indjg,val,(Nn),(Nn));
            clear val indig indjg
        end
    end
end
