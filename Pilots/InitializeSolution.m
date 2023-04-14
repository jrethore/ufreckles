function [Uini]=InitializeSolution(U0,iscale,nmod,imin)
load(fullfile('TMP','params'));
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
dotopo=(~(isfield(param,'topography'))&&isfield(param0,'calibration_data'));
restart=1;
if isfield(param0,'restart')
    restart=param0.restart;
end
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
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
    reader=VideoReader(param0.reference_image);
    if isfield(param0,'number_of_frames')
        nbf=param0.number_of_frames;
    else
        nbf=reader.NumberOfFrames-1;
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
nz=numel(indcam);
if isempty(U0)
    do_rbt=1;Urbt=zeros([2,nim,nz]);U0=[0;0];
    if ~isfield(param0,'deformed_image'),do_rbt=0;end
    %    indcam=1;
elseif size(U0,1)==2
    do_rbt=1;Urbt=U0;
    %    indcam=1;
else
    do_rbt=0;
end
Uini=[];
im1=[];
if restart
    nmax=nim;
else
    if nim<2+dotopo
        nmax=nim;
    else
        nmax=2+dotopo;
    end
end

if isfield(param0,'time_step')
    dt=param0.time_step;
    assert(nim>dt);
    nmax=dt+dotopo;
end
if nargin<4
    imin=1;
else
nmax=imin;
end

for iim=imin:nmax
    for icamr=1:ncamr
        load(fullfile('TMP',sprintf('sample%d_0',icamr-1)),'sizeim');
        if do_rbt
        load(fullfile('TMP',sprintf('sample%d_%d',icamr-1,param.nscale-1)),'im0');
        im0n=im0;
        load(fullfile('TMP',sprintf('sample%d_0',icamr-1)),'im0','roi');
        end
        if icamr==1
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,1-1)),'mask');
            if numel(mask)>1
                mask=reshape(full(diag(mask)),sizeim);
            end
        else
            mask=1;
        end
        for icamd=1:size(indcam,2)
            iz=indcam(icamr,icamd)+(icamr-1)*size(indcam,2)*docrosscorrelation;
            if do_rbt
                if nim==1&&ncamd==1
                    fildef=param0.deformed_image;
                else
                    fildef=param0.deformed_image{indcam(icamr,icamd),iim};
                end
                im1=double(readim(fildef));
                if length(size(im1))==3
                    im1=mean(im1,3);
                end
                if reverse
                    im1=im1';
                end
                 if ~any(Urbt(:,iim,iz))
                     Uni=0;Vni=0;
%                      if iim>1
%                          Vni=round(Urbt(2,iim-1,iz)/2^(param.nscale-1));
%                          Uni=round(Urbt(2,iim-1,iz)/2^(param.nscale-1));
%                      end
                im1n=MCoarseImage(im1(Uni+(roi(1):roi(2)),Vni+(roi(3):roi(4))),2^(param.nscale-1));

                
                    [Un,Vn]=rbt(im0n,im1n);
                   
                    Urbt(2,iim,iz)=round((Vn+Vni)*2^(param.nscale-1));
                    Urbt(1,iim,iz)=round((Un+Uni)*2^(param.nscale-1));
                end
    %           if numel(im0)<1e6
    %            [U,V]=rbt(im0,im1((roi(1):roi(2))+Urbt(1,iim,iz),(roi(3):roi(4))+Urbt(2,iim,iz)));
    %            else
                    U=0;V=0;
    %            end
                U0=[U+Urbt(1,iim,iz);V+Urbt(2,iim,iz)];
                
            end
            sbasis=param.basis;
            
            switch sbasis
                
                case 'rbt'
                    if isempty(Uini)
                        Uini=zeros([2,nim,nz]);
                    end
                    Uini(1,iim,iz)=U0(1);
                    Uini(2,iim,iz)=U0(2);
                case 'rbm'
                    if isempty(Uini)
                        Uini=zeros([3,nim,nz]);
                    end
                    Uini(1,iim,iz)=U0(1);
                    Uini(2,iim,iz)=U0(2);
                case 'disto'
                    if isempty(Uini)
                        Uini=zeros([5,nim,nz]);
                    end
                    Uini(1,iim,iz)=U0(1);
                    Uini(2,iim,iz)=U0(2);
                case 'cordin'
                    if size(U0,1)==2
                        if isempty(Uini)
                            Uini=zeros([16,nim,nz]);
                        end
                        Uini(1,iim,iz)=U0(1);
                        Uini(2,iim,iz)=U0(2);
                    else
                        Uini=U0;
                    end
                case 'uni'
                    if size(U0,1)==2
                        if isempty(Uini)
                            Uini=zeros([6,nim,nz]);
                        end
                        Uini(1,iim,iz)=U0(1);
                        Uini(2,iim,iz)=U0(2);
                    else
                        Uini=U0;
                    end
                    
                case 'affine'
                    if size(U0,1)==2
                    if isempty(Uini)
                        Uini=zeros([4,nim,nz]);
                    end
                    Uini(1,iim,iz)=U0(1);
                    Uini(2,iim,iz)=U0(2);
                    else
                        Uini=U0;
                    end
                        
                case 'beam'
                    if isempty(Uini)
                        Uini=zeros([6,nim,nz]);
                    end
                    Uini(1,iim,iz)=U0(2);
                    Uini(2,iim,iz)=U0(2);
                    Uini(3,iim,iz)=U0(1);
                    Uini(4,iim,iz)=0;
                    Uini(5,iim,iz)=U0(1);
                    Uini(6,iim,iz)=0;
                case 'beam-nurbs'
                    assert(iscale==1);
                    if iim==1
                        load(fullfile('TMP',sprintf('%d_mphi_%d',nmod,iscale-1)),'mphi');
                    end
                    if isempty(Uini)
                        Uini=zeros([length(mphi),nim,nz]);
                    end
                    Uini(:,iim,iz)=mphi*U0(1)+(1-mphi)*U0(2);
                case 'nurbs-beam'
                    Uinii=ProjectDisplacement(U0,iscale,nmod,im1,icamr);
                    Uinii=LSProjectDisplacement(Uinii,iscale,nmod);
                    if isempty(Uini)
                        Uini=zeros([length(Uinii),nim,nz]);
                    end
                    Uini(:,iim,iz)=Uinii;
                    
                    
                case 'KM'
                    if isfield(param,'crack_id')
                        ic=param.crack_id;
                    else
                        ic=1;
                    end
                    ind=param.km_indices;
                    modes=param.modes;
                    kappa=param.kolosov;
                    load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'theta');
                    if isempty(Uini)
                        Uini=zeros([length(ind)*length(modes),nim,nz]);
                    end
                    found=find(ind==0);
                    U0=U0(1)+i*U0(2);
                    U0=U0*exp(-i*theta);
                    for m=1:length(modes)
                        if modes(m)==1
                            Uini(found+(modes(m)-1)*length(ind),iim,iz)=real(U0)/(1+kappa);
                        elseif modes(m)==2
                            Uini(found+(modes(m)-1)*length(ind),iim,iz)=-imag(U0)/(1+kappa);
                        end
                    end
                case {'KM+','CZ'}
                    load(fullfile('TMP',sprintf('%d_phix_0',nmod)),'Nddl_tot');
                    ind=[0,2];
                    modes=param.modes;
                    kappa=param.kolosov;
                    if isfield(param,'crack_id')
                        ic=param.crack_id;
                    else
                        ic=1;
                    end
                    load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'theta');
                    
                    if isempty(Uini)
                        Uini=zeros([Nddl_tot,nim,nz]);
                    end
                    found=find(ind==0);
                    U0=U0(1)+i*U0(2);
                    U0=U0*exp(-i*theta);
                    for m=1:length(modes)
                        if modes(m)==1
                            Uini(found+(modes(m)-1)*length(ind),iim,iz)=real(U0)/(1+kappa);
                        elseif modes(m)==2
                            Uini(found+(modes(m)-1)*length(ind),iim,iz)=-imag(U0)/(1+kappa);
                        end
                    end
                    
                case 'fem'
                    if size(U0,1)>2
                        U00=U0(:,iim,iz);
                    else
                        U00=U0;
                    end
                    
                    Uinii=ProjectDisplacement(U00,iscale,nmod,im1,icamr);
                    if (iscale==1)
                        if (isfield(param,'enrichment'))
                            load(fullfile('TMP',sprintf('%d_phix_%d',nmod,iscale-1)),'Nddl_tot');
                            Uinii(length(Uinii)+1:Nddl_tot)=0;
                            %                             if ~iscell(param0.levelset_file)
                            %                                 nbfis=1;
                            %                             else
                            %                                 nbfis=numel(param0.levelset_file);
                            %                             end
                            %                             mesh_file=fullfile('TMP',sprintf('%d_meshx_%d',nmod,iscale-1));
                            %                                 load(mesh_file,'xo','yo','Nnodes');
                            %                             figure
                            %                             imagesc(reshape(Uinii(1:prod(Nnodes)),Nnodes))
                            %
                            %                             dec=0;
                            %                             for ic=1:nbfis
                            %                                 load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'face_nodes');
                            %                                 Uface=Uinii([face_nodes,face_nodes+prod(Nnodes)]);
                            %                                 Uinii([face_nodes,face_nodes+prod(Nnodes)])=0;
                            %                                 load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack');
                            %                                 hn=interp2(crack,yo,xo,'*linear');
                            %                                 hn=2*double(hn(face_nodes)>=0)-1;
                            %                                 Uinii(2*prod(Nnodes)+dec+(1:2*length(face_nodes)))=[hn;hn].*Uface;
                            %                                 dec=dec+2*length(face_nodes);
                            %                             figure
                            %                             imagesc(reshape(Uinii(1:prod(Nnodes)),Nnodes))
                            %                             Uo=zeros(2*prod(Nnodes),1);
                            %                             Uo([face_nodes,face_nodes+prod(Nnodes)])=Uface;
                            %                             figure
                            %                             imagesc(reshape(Uo(1:prod(Nnodes)),Nnodes))
                            %
                            %                             end
                        elseif isfield(param,'gp_file')
                            
                            Uinii=LSProjectDisplacement(Uinii,iscale,nmod);
                        end
                    end
                    if isempty(Uini)
                        Uini=zeros([length(Uinii),nim,nz]);
                    end
                    Uini(:,iim,iz)=Uinii;
                case {'nurbs','btri'}
                    if length(U0)>2
                        U00=U0(:,iim,iz);
                    else
                        U00=U0;
                    end
                    Uinii=ProjectDisplacement(U00,iscale,nmod,im1,icamr);
                    
                    if iscale==1
                        Uinii=LSProjectDisplacement(Uinii,iscale,nmod);
                    end
                    if isempty(Uini)
                        Uini=zeros([length(Uinii),nim,nz]);
                    end
                    Uini(:,iim,iz)=Uinii;
                    
                otherwise
                    error ('INVALID FUNCTIONAL BASIS');
                    
            end
        end
    end
    
end

end