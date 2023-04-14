function [U]=SolveST(Uini,iscale,nmod)
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
dotopo=(~isfield(param,'topography'))&&isfield(param0,'calibration_data');
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
indcam=ncamd:-1:1;
disct=0;
maxiter=param0.iter_max;
conv=param0.convergance_limit;
dt=param0.time_step;
disp(sprintf('Starting resolution for scale %d...',iscale));
pscale=2^(iscale-1);
im1=0;
U=Uini;
fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
if iscale==1
    nmax=nim;
else
    nmax=dt+dotopo;
end
ti=0:dt:nmax;
ti(1)=1;
if ti(end)<nmax
    if nmax-ti(end)>0.5*dt
        ti=[ti,nmax];
    else
        ti(end)=nmax;
    end
end
nbstep=length(ti)-1;
icamr=1;
load(fullfile('TMP',sprintf('sample%d_%d',(icamr-1),iscale-1)),'im0','sizeim','roi');
dynamic=max(im0(:))-min(im0(:));
ndim=length(sizeim);

load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),(iscale-1))),'on');
load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),10*(iscale-1))),'phix','Xi','Yi','wdetJ');
load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icamr-1),10*(iscale-1))),'phiy');
load(fullfile('TMP',sprintf('%d_operator_%d',nmod*10^(icamr-1),iscale-1)),'M','R','P');
M=M+R;
load(fullfile('TMP',sprintf('%d_phidf_%d',nmod*10^(icamr-1),iscale-1)),'phidf');
load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
fwrite(fid,0);fwrite(fid,min(255,dynamic));
if iscale==1
    if strcmp(param.basis,'fem')||strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri')
        mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1));
        load(mesh_file,'ng');
        if (ng>0)||isfield(param,'gp_file')||isfield(param,'extrusion_parameters')||strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri')
            load(fullfile('TMP',sprintf('sample%d',(icamr-1))),'im0','sizeim');
            [im00]=mexInterpSpline((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
            im0=im00;
        end
    end
end
mean0=mean(im0(:));
std0=std(im0(:));
im0=im0-mean0;

for icamd=1:size(indcam,2)
    display(sprintf('Camera %d/%d',indcam(icamd),ncamd));
    iz=indcam(icamr,icamd);
    Ut=Uini(:,1,iz);
    
    res=Inf;
    ii=1;
    display(sprintf('Initialization'));
    
    while ( res>conv && ii< maxiter)
        
        Ux=phix*Ut;
        Uy=phiy*Ut;
        NestedResidual(1);
        disc=mask*disc;
        merror=std(disc(on))/dynamic;
        
        F=(phidf'*(wdetJ*disc));
        
        F=F-R*Ut;
        dU=M\F;
        disp(sprintf('At iteration # %d',ii));
        disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
        if ii==1
            nU0=max(sqrt(length(Ut)),norm(Ut+dU(:)));
        else
            if length(dU)>10
                res=norm(dU(:))/nU0;
            else
                res=max(abs(dU)./Ut);
            end
            disp(sprintf('|dU|=%f',res));
            
        end
        Ut=Ut+dU;
        ii=ii+1;
    end
    if dynamic>255
        disc=255*disc/dynamic;
    end
    fwrite(fid,(abs(disc)));
    U(:,1,iz)=Ut;
    DU=U(:,ti(2))-Ut;
    for istep=1:nbstep
        im_step=(ti(istep)+1):ti(istep+1);
        Nt=(im_step-ti(istep))/(ti(istep+1)-ti(istep));
        display(sprintf('Step %d/%d',istep,nbstep));
        Ui=Ut*ones(1,length(im_step))+DU*Nt;
        Ut=Ui(:,end);
%         figure
%         plot(1:im_step(end),Uini(:,1:im_step(end),iz)','b')
%         hold on
%         plot(im_step,Ui','r')
%         keyboard
        res=Inf;
        ii=1;
        nNt=Nt*Nt';
        nNt=sum(Nt);
        while ( res>conv && ii< maxiter)
            
            F=0;merror=0;
            for iim=1:length(im_step)
                Ux=phix*Ui(:,iim);
                Uy=phiy*Ui(:,iim);
                NestedResidual(im_step(iim));
                disc=mask*disc;
                merror=merror+std(disc(on))/dynamic;
                
                F=F+(phidf'*(wdetJ*disc))*(Nt(iim)/nNt);
                
            end
            F=F-R*Ut;
            merror=merror/length(im_step);
            dU=M\F;
            disp(sprintf('At iteration # %d',ii));
            disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
            if ii==1
                nU0=max(sqrt(length(Ut)),norm(Ut+dU(:)));
            else
                if length(dU)>10
                    res=norm(dU(:))/nU0;
                else
                    res=max(abs(dU)./Ut);
                end
                disp(sprintf('|dU|=%f',res));
                
            end
            Ut=Ut+dU;
            Ui=Ui+dU*Nt;
            ii=ii+1;
        end
        U(:,im_step,iz)=Ui;
        DU=Ut-U(:,ti(istep),iz);
        for iim=1:length(im_step)
            Ux=phix*Ui(:,iim);
            Uy=phiy*Ui(:,iim);
            NestedResidual(im_step(iim));
            disc=mask*disc;
            if dynamic>255
                disc=255*disc/dynamic;
            end
            fwrite(fid,(abs(disc)));
        end
    end
end


disp(sprintf('Enlapsed time for resolution = %6.2f s',cputime -ttic));
fclose(fid);
    function NestedResidual(kkk)
        if ~isfield(param0,'deformed_image')
            im1=double(read(reader,frames(kkk)));
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
        
        %                                   if iscale==1
        %                                            figure
        %                                            scatter(Xi,Yi,10+0*Uy,Uy)
        %         %                                  imagesc(reshape(Ux,sizeim*pscale))
        %                                            pause
        %                                   end
        sizeim1=size(im1);
        disc=(mexInterpSpline((Xi-1)+Ux+roi(1),(Yi-1)+Uy+roi(3),im1));
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
