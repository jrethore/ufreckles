function [U,disct]=Solve(Uini,iscale,nmod,Fext)
clear MedianFilter
pk=[];
if nargin<4,dokk=1;
else dokk=0;end
ttic=cputime;
pscale=2^(iscale-1);
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
onflight=0;
if isfield(param0,'onflight')
    onflight=param0.onflight;
end
preview=0;
if isfield(param0,'preview')
    preview=param0.preview;
end
if onflight
    global phiy phix Xi Yi wdetJ inde on phidf
end
mfilter=0;
if isfield(param0,'regularization_type')
    if strcmp(param0.regularization_type,'median')
        mfilter=1;
    elseif strcmp(param0.regularization_type,'tiko')&&strcmp(param.basis,'fem')
        mfilter=1;
        param0.regularization_parameter=pscale;
    end
end
pgd=0;
if isfield(param0,'do_pgd_prediction')
    pgd=param0.do_pgd_prediction;
end
if isfield(param0,'time_step')
    if iscale==1
        disct=0;
        [U]=SolveST(Uini,iscale,nmod);
        return
    else
        param0.restart=0;
    end
end
computefint=true;
if isfield(param0,'compute_fint')
    computefint=param0.compute_fint;
end
dotopo=(~isfield(param,'topography'))&&isfield(param0,'calibration_data');
incremental=0;
if isfield(param0,'incremental_resolution')
    incremental=param0.incremental_resolution;
    if iscale>1
        incremental=0;
    end
end
psample=1;
if isfield(param0,'sampling_factor')
    psample=param0.sampling_factor;
end
if iscell(param0.reference_image)
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
    mreader=VideoReader(param0.reference_image);
    nbf=mreader.NumberOfFrames-1;
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
restart=1;
if isfield(param0,'restart')
    restart=param0.restart;
end
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
disct=0;
maxiter=param0.iter_max;
conv=param0.convergance_limit;
disp(sprintf('Starting resolution for scale %d...',iscale));


im1=0;
pixint=1;
U=Uini;
fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
if restart
    nmax=nim;
else
    if iscale==1||preview
        nmax=nim;
    else
        if nim<2+dotopo
            nmax=nim;
        else
            nmax=2+dotopo;
        end
        if isfield(param0,'time_step')
            nmax=param0.time_step;
        end
    end
end
for ijm=1:nmax
    display(sprintf('Image %d/%d',ijm,nmax));
    for icamr=1:ncamr
        
        if ijm==1||ncamr>1||incremental
            load(fullfile('TMP',sprintf('sample%d_%d',(icamr-1),iscale-1)),'im0','sizeim','roi');
            dynamic=max(im0(:))-min(im0(:));
            ndim=length(sizeim);
            roip=(roi-1)/pscale+1;
            
            if ~onflight
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),(iscale-1))),'on');
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),(iscale-1))),'phix','Xi','Yi','wdetJ');
                load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icamr-1),(iscale-1))),'phiy');
                load(fullfile('TMP',sprintf('%d_phidf_%d',nmod*10^(icamr-1),iscale-1)),'phidf');
                
            end
            load(fullfile('TMP',sprintf('%d_operator_%d',nmod*10^(icamr-1),iscale-1)),'M','R','P');
            Mo=M;
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
            if ijm==1&&icamr==1, fwrite(fid,0);fwrite(fid,min(255,dynamic));end
            if iscale==1
                if strcmp(param.basis,'fem')||strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri')
                    mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1));
                    load(mesh_file,'ng');
                    if (ng>0)||isfield(param,'gp_file')||isfield(param,'extrusion_parameters')||strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri')
                        load(fullfile('TMP',sprintf('sample%d',(icamr-1))),'im0','sizeim');
                        [im00]=mexInterpSpline((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
                        im0=im00;
                        pixint=0;
                    end
                end
            end
            mean0=mean(im0(on));
            std0=std(im0(on));
            im0=im0-mean0;
        end
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
                            %                            load('TMP/0_mesh_0','xo','yo')
                            L=U(:,max(1,ijm-4):(ijm-1),iz);
                            
                            %                             Lp=[1+0*xo,xo,yo,xo.^2,yo.^2,xo.*yo,zeros(length(xo),6);...
                            %                                 zeros(length(xo),6),1+0*xo,xo,yo,xo.^2,yo.^2,xo.*yo];
                            %                             L=[L,Lp];
                            dL=Inf;
                            ii=1;
                            MM=L'*(Mo)*L;
                            while norm(dL)>0.001&&ii<100
                                Ux=phix*Ui;
                                Uy=phiy*Ui;
                                NestedResidual(ijm);
                                disc=mask*disc;
                                FF=phidf'*(wdetJ*disc);
                                FF=L'*FF;
                                dL=MM\FF;
                                Ui=Ui+L*dL;
                                ii=ii+1;
                            end
                            ii=1;
                        end
                    end
                    
                end
                
            end
            if isfield(param,'coupling_parameter')&&(iscale==1)
                if ijm==1
                    LoadMat(nmod);
                    Uprev=0*Ui;
                    Fprev=0*Ui;
                else
                    Uprev=U(:,ijm-1,iz);
                    Fprev=Fi;
                end
                AssembleMechanicalOperator(iscale,nmod,ijm);
                load(fullfile('TMP',sprintf('%d_k_operator_%d',nmod,iscale-1)),'K');
                if dokk
                    if ijm==1&&ii==1
                        if exist(fullfile('TMP',sprintf('%d_imodel.mat',nmod)),'file')
                            load(fullfile('TMP',sprintf('%d_imodel',nmod)),'select');
                        else
                            AssembleEquilibriumGapOperator(iscale,nmod)
                            load(fullfile('TMP',sprintf('%d_kk_operator_%d',nmod,iscale-1)),'select');
                        end
                    end
                    K=select*K;
                    if isfield(param0,'penalty_parameter')
                        alpha=param0.penalty_parameter;
                    else
                        alpha=0e6;
                    end
                    Pb=alpha*diag(~diag(select));
                end
                l=param.coupling_parameter;
            end
            merrorp=Inf;
            while ( res>conv && ii< maxiter)
                dalpha=1+0*exp(-1*(ii-1)/maxiter);
                Ux=phix*Ui;
                Uy=phiy*Ui;
                NestedResidual(ijm);
                disc=mask*disc;
                merror=std(disc(on))/dynamic;
                
                F=phidf'*(wdetJ*disc);
                
                if isfield(param,'coupling_parameter')&&(iscale==1)
                    if dokk
                        if ii==1
                            if ijm==1
                                mo=disc'*disc;
                                ko=Ui'*K'*K*Ui;
                                %ko=1
                            end
                            M=l*Mo/mo+(1-l)*K'*K/ko+Pb;
                        end
                        if (~computefint)
                            Fi=Fprev+K*(Ui-Uprev);
                        else
                            ComputeNewInternalState(nmod,Ui);
                            Fi=ComputeInternalForceVector(nmod,Ui,ijm);
                        end
                        F=l*F/mo-(1-l)*K'*(select*Fi)/ko;
                        if ii==1
                            Fio=select*Fi;
                            normFio=Fio'*Fio;
                        else
                            display(sprintf('F/Fo %f',full(sqrt(((select*Fi)'*(select*Fi))/normFio))));
                        end
                    else
                        if (~computefint)
                            Fi=Fprev+K*(Ui-Uprev);
                        else
                            ComputeNewInternalState(nmod,Ui);
                            Fi=ComputeInternalForceVector(nmod,Ui,ijm);
                        end
                        if ii==1
                            Fe=Fext(:,ijm);
                            if ijm==1
                                mo=disc'*disc;
                                ko=norm(Fe-Fi);
                            end
                            M=l*Mo/mo+(1-l)*K/ko;
                        end
                        F=l*F/mo+(1-l)*(Fe-Fi)/ko;
                        if ii==1
                            normFio=norm(Fe-Fi);
                        else
                            display(sprintf('R/Ro %f',norm(Fe-Fi)/normFio));
                        end
                        
                        
                        
                    end
                else
                    M=Mo+dalpha*(R+P);
                    F=F-dalpha*R*Ui;
                    
                end
                if isfield(param0,'solver')
                    solv=param0.solver;
                    
                    if strcmp(solv,'lu')
                        if ii==1
                            [LT,UT]=lu(sparse(M));
                        end
                        dU=UT\(LT\F);
                    elseif strcmp(solv,'chol')
                        if ii==1
                            Rc=chol(sparse(M));
                        end
                        dU=Rc\(Rc'\F);
                    end
                else
                    try
                    if isempty(pk)
                        pk=symamd(M);
                        %                            [LM,UM]=lu(M(pk,pk));
                        [LM]=chol(M(pk,pk));
                        dU=zeros(size(M,1),1);
                    end
                    %                        dU(pk)=UM\(LM\F(pk));
                    dU(pk)=LM\(LM'\F(pk));
                    catch
                        dU=M\F;
                    end
                end
                %                  if merror>merrorp
                %                      Ui=Ui-dU;
                %                      if mfilter
                %                          mesh_file=fullfile('TMP',sprintf('%d_mesh_%d.mat',nmod,iscale-1));
                %                          lc=ceil(param0.regularization_parameter/pscale);
                %                          Ui=MedianFilter(Ui,mesh_file,lc);
                %
                %
                %                      end
                %                      break
                %                  end
                
                disp(sprintf('At iteration # %d',ii));
                disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
                if ii==1
                    nU0=max(sqrt(length(Ui)),norm(Ui+dU(:)));
                else
                    if length(dU)>10
                        res=norm(dU(:))/nU0;
                    else
                        res=max(abs(dU)./abs(Ui));
                    end
                    disp(sprintf('|dU|=%f',res));
                    if mfilter
                        res=min(res,abs(merrorp-merror)/abs(merror));
                    end
                    
                end
                Ui=Ui+dU;
                if mfilter
                    mesh_file=fullfile('TMP',sprintf('%d_mesh_%d.mat',nmod,iscale-1));
                    %                 if iscale==1
                    %                     load(mesh_file,'xo','yo')
                    %                     figure
                    %                     scatter(xo,yo,10+0*xo,Ui(length(xo)+(1:length(xo))))
                    %                 end
                    lc=ceil(param0.regularization_parameter/pscale);
                    Ui=MedianFilter(Ui,mesh_file,lc);
                    %                 Ui=Ui-dU;
                    %                     dU=MedianFilter(dU,mesh_file,lc);
                    %                 Ui=Ui+dU;
                    
                    %                 if iscale==1
                    %                     figure
                    %                     scatter(xo,yo,10+0*xo,Ui(length(xo)+(1:length(xo))))
                    % pause
                    % close all
                    %                 end
                    
                end
                %display(sprintf('uncertainty at iteration %d : %f\n',ii,sqrt(0.5*(var(Ui((1:size(U,1)/2)))+var(Ui(size(U,1)/2+(1:size(U,1)/2)))))))
                ii=ii+1;
                merrorp=merror;
            end
            if isfield(param,'coupling_parameter')&&(iscale==1)
                if ~computefint
                    ComputeNewInternalState(nmod,Ui,0);
                    Fi=ComputeInternalForceVector(nmod,Ui,ijm);
                end
                UpdateInternalState(nmod,Ui,ijm);
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
            disct=disct+disc'*mask*disc;
            if dynamic>255
                disc=255*disc/dynamic;
            end
            disc(~on)=0;
            fwrite(fid,(abs(disc)));
            if pixint&&iscale==1
                VTKExportIntMap(sprintf('camr-%d-camd-%d-scale-%d-%d',icamr,icamd,iscale,ijm),'error',roi(1)-1+(1:sizeim(1)),roi(3)-1+(1:sizeim(2)),reshape(uint8(abs(disc)),sizeim),1);
            end
        end
    end
end
disp(sprintf('Enlapsed time for resolution = %6.2f s',cputime -ttic));
fclose(fid);
    function NestedResidual(kkk)
        if ii==1
            if ~isfield(param0,'deformed_image')
                im1=double(readim(mreader,frames(kkk)));
            else
                if ~iscell(param0.deformed_image)
                    fildef=param0.deformed_image;
                    im1=double(readim(fildef));
                else
                    fildef=param0.deformed_image{indcam(icamr,icamd),kkk};
                    im1=double(readim(fildef));
                end
            end
            if length(size(im1))==3
                im1=mean(im1,3);
            end
            if reverse
                im1=im1';
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
        %                                   if iscale==1
        %                                            figure
        %                                            scatter(Xi,Yi,10+0*Uy,Uy)
        %         %                                  imagesc(reshape(Ux,sizeim*pscale))
        %                                            pause
        %                                   end
        sizeim1=size(im1);
        %        disc=(mexInterpSpline((Xi-1)*psample+Ux+roi(1),(Yi-1)*psample+Uy+roi(3),im1));
        if iscale==1
            disc=(mexInterpSpline((Xi-1)*psample+Ux/pscale+roip(1),(Yi-1)*psample+Uy/pscale+roip(3),im1));
        else
            disc=(mexInterpLinear((Xi-1)*psample+Ux/pscale+roip(1),(Yi-1)*psample+Uy/pscale+roip(3),im1));
        end
        %        disc=mexInterpSpline((Xi+Ux+roi(1)-1),(Yi+Uy+roi(3)-1),im1);
        %         if iscale>1
        %             disc=reshape(disc,sizeim*pscale);
        %             NestedCoarseImage(pscale);
        %         end
        
        maske=disc<0;
        maske(~on)=1;
        mean1=mean(disc(~maske));
        std1=std(disc(~maske));
        disc=disc-mean1;
        st=std0/std1;
        disc=(im0(:)-st*disc(:));
        disc(maske)=0;
        
        
        function NestedCoarseImage(scale)
            
            imsiz0=size(disc);
            imsiz1=floor(imsiz0/2);
            nn=2*imsiz1;
            ndim=length(nn);
            if ndim==2
                disc=disc(1:nn(1),1:nn(2));
            elseif ndim==3
                disc=disc(1:nn(1),1:nn(2),1:nn(3));
            end
            
            disc=reshape(disc,scale,prod(nn)/scale);
            disc=mean(disc,1);
            nn(1)=nn(1)/scale;
            disc=reshape(disc,nn);
            
            if ndim==2
                disc=disc';
                disc=reshape(disc,scale,prod(nn)/scale);
                disc=mean(disc,1);
                nn(2)=nn(2)/scale;
                disc=reshape(disc,nn([2,1]));
                disc=disc';
                %toc();
            elseif ndim==3
                disc=permute(disc,[2,3,1]);
                disc=reshape(disc,scale,prod(nn)/scale);
                disc=mean(disc,1);
                nn(2)=nn(2)/scale;
                disc=reshape(disc,nn([2,3,1]));
                disc=permute(disc,[3,1,2]);
                
                disc=permute(disc,[3,1,2]);
                disc=reshape(disc,scale,prod(nn)/scale);
                disc=mean(disc,1);
                nn(3)=nn(3)/scale;
                disc=reshape(disc,nn([3,1,2]));
                disc=permute(disc,[2,3,1]);
            end
            
        end
    end
end
