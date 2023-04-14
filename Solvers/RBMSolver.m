function [U,A]=RBMSolver(Uini,iscale,nmod,L)
if nargin<4,dokk=1;
else dokk=0;end
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
onflight=0;
if isfield(param0,'onflight')
    onflight=param0.onflight;
end
if onflight
    global phiy phix Xi Yi wdetJ inde on phidf
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
psample=1;
if isfield(param0,'sampling_factor')
    psample=param0.sampling_factor;
end
ncamr=1;
dotopo=(~isfield(param,'topography'))&&isfield(param0,'calibration_data');
docrosscorrelation=0;
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

pscale=2^(iscale-1);

im1=0;
pixint=1;
U=Uini;
A=zeros(size(L,2),size(U,2));
fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
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
        if isfield(param0,'time_step')
            nmax=param0.time_step;
        end
    end
end
for ijm=1:nmax
    display(sprintf('Image %d/%d',ijm,nmax));
    for icamr=1:ncamr
        
        if ijm==1||ncamr>1
            load(fullfile('TMP',sprintf('sample%d_%d',(icamr-1),iscale-1)),'im0','sizeim','roi');
            dynamic=max(im0(:))-min(im0(:));
            ndim=length(sizeim);
            if ~onflight
                
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),(iscale-1))),'on');
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),10*(iscale-1))),'phix','Xi','Yi','wdetJ');
                load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icamr-1),10*(iscale-1))),'phiy');
                load(fullfile('TMP',sprintf('%d_phidf_%d',nmod*10^(icamr-1),iscale-1)),'phidf');
            end
            load(fullfile('TMP',sprintf('%d_operator_%d',nmod*10^(icamr-1),iscale-1)),'M','R','P');
            M=L'*((M+R+P)*L);
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
            mean0=mean(im0(:));
            std0=std(im0(:));
            im0=im0-mean0;
        end
        for icamd=1:size(indcam,2)
            iz=indcam(icamr,icamd)+(icamr-1)*size(indcam,2)*docrosscorrelation;
            res=1;
            ii=1;
            if restart
                Ui=Uini(:,ijm,iz);
            Pi=zeros(size(L,2),1);
            else
                if ijm==1
                    Ui=Uini(:,ijm,iz);
            Pi=zeros(size(L,2),1);
                    if dotopo&&(indcam(icamr,icamd)==icamr)
                        Ui=0*Ui;
                        res=0;
                        disc=zeros(numel(im0),1);
                    end
                elseif ijm<3+dotopo
                    Ui=Uini(:,ijm,iz);
                else
%                    DU=U(:,ijm-1,iz)-U(:,ijm-2,iz);
                    DP=A(:,ijm-1)-A(:,ijm-2);
%                    if incremental
%                        Ui=DU;
%                    else
%                        Ui=U(:,ijm-1,iz)+DU;
                        Pi=A(:,ijm-1)+DP;
 Ui=Uini(:,ijm,iz)+L*Pi;
 %                    end
                    
                end
                
            end
            
            while ( res>conv && ii< maxiter)
                                
                Ux=phix*Ui;
                Uy=phiy*Ui;
                NestedResidual(ijm);
                disc=mask*disc;
                merror=std(disc(on))/dynamic;
                
                F=phidf'*(wdetJ*disc);
                F=L'*(F-R*Ui);
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
                    dU=M\F;
                end
                
                Pi=Pi+dU;
                dU=L*dU;
                
                disp(sprintf('At iteration # %d',ii));
                disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
                if ii==1
                    nU0=max(sqrt(length(Ui)),norm(Ui+dU(:)));
                else
                    if length(dU)>10
                        res=norm(dU(:))/nU0;
                    else
                        res=max(abs(dU)./Ui);
                    end
                    disp(sprintf('|dU|=%f',res));
                    
                end
                Ui=Ui+dU;
                %display(sprintf('uncertainty at iteration %d : %f\n',ii,sqrt(0.5*(var(Ui((1:size(U,1)/2)))+var(Ui(size(U,1)/2+(1:size(U,1)/2)))))))
                ii=ii+1;
            end
            U(:,ijm,iz)=Ui;
            A(:,ijm)=Pi;
            disct=disct+disc'*mask*disc;
            if dynamic>255
                disc=255*disc/dynamic;
            end
            fwrite(fid,(abs(disc)));
            if pixint&&iscale==1
                disc(~on)=0;
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
            
        end
        %                                   if iscale==1
        %                                            figure
        %                                            scatter(Xi,Yi,10+0*Uy,Uy)
        %         %                                  imagesc(reshape(Ux,sizeim*pscale))
        %                                            pause
        %                                   end
        sizeim1=size(im1);
        disc=(mexInterpSpline((Xi-1)*psample+Ux+roi(1),(Yi-1)*psample+Uy+roi(3),im1));
        %        disc=mexInterpSpline((Xi+Ux+roi(1)-1),(Yi+Uy+roi(3)-1),im1);
        if iscale>1
            disc=reshape(disc,sizeim*pscale);
            NestedCoarseImage(pscale);
        end
        
        maske=disc<0;
        mean1=mean(disc(:));
        std1=std(disc(:));
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
