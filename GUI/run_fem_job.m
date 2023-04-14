function run_fem_job(param,model,Up)
nmod=0;
if nargin<3,Up=[];end
clear functions
pshift=0;
if isfield(param,'deformed_image')
if iscell(param.deformed_image)
    if size(param.deformed_image,1)>2
        pshift=size(param.deformed_image,1);
        defs=param.deformed_image;
        param.deformed_image=defs(1,:);
        if size(defs,2)==1,param.deformed_image=defs{1,1};end
    end
end
end
param.onflight=1;
LoadParameters(param);
LoadParameters(model,nmod);
ReferenceImage(nmod);
LoadMask(nmod);
LoadMeshes(nmod);
LoadMat(nmod);

nscale=model.nscale;
filres=param.result_file;
%%
load(fullfile('TMP','sample0_0.mat'),'sizeim')
if prod(sizeim)>4e6
    param.initialization='bilin+';
    %    param.initialization='rbt_elem';
    LoadParameters(param);
end
load(fullfile('TMP','sample0_0.mat'),'sizeim')
% if isfield(param,'detect')
%     if param.detect
%      param.initialization='rbt';
%     LoadParameters(param);
%     end
% end
load(fullfile('TMP','0_mesh_0.mat'),'selected','Nnodes')
if isfield(model,'mesh_type')
    if model.mesh_type==1
        rflag=true;
        save(fullfile('TMP','0_mesh_0.mat'),'rflag','-append')
    end
end
extract=0;
if strcmp(param.regularization_type,'equilibrium_gap')
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        switch zone{4}
            %            case 6
            %                selected(zone{6})=0;
            case 5
                if (zone{8}>0)
                    extract=1;
                end
        end
    end
end
tips=zeros(size(model.zone,2),2);
cracks=[];
inds=[];
if extract
    load(fullfile('TMP','0_mesh_0.mat'),'xo','yo','conn','elt')
    if param.detect
        modes=1:2;
        harms=-3:7;
        
    else
        modes=1:2;
        harms=0:7;
    end
    roi=param.roi;
    freenodes=ones(length(xo),1);
    xo=xo+roi(1)-1;
    yo=yo+roi(3)-1;
    L=zeros(prod(Nnodes),1);
    nu=model.material_parameters.nu;
    kappa=(3-4*nu);
    %    kappa=(3-nu)/(1+nu)
    
    zo=xo+1i*yo;
    nw=0;
    ntip=0;
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        if zone{4}==5
            if (zone{8}>0)
                tips(iz,1)=ntip+1;
                indc=zone{10};
                cnodes=indc{1};
                ztip=zo(cnodes(1,1));
                tt=-diff(zo(cnodes(1:2,1)));
                tt=tt/abs(tt);
                rtip=mean(abs(zo(indc{2})-ztip));
                
                %                rtip=Inf
                
                in=abs(zo-ztip)<=rtip;
                in(indc{2})=1;
                freenodes(in)=0;
                
                
                selected(in)=0;
                selected(indc{2})=1;
                
                
                zzone=(zo(in)-ztip)*tt';
                idin=1:prod(Nnodes);
                idin=idin(in);
                for ip=1:length(zzone)
                    if any(cnodes(:,1)==idin(ip))
                        zzone(ip)=abs(zzone(ip))*exp(1i*pi);
                    end
                    if any(cnodes(:,2)==idin(ip))
                        zzone(ip)=abs(zzone(ip))*exp(-1i*pi);
                    end
                end
                phi=Williams(zzone,modes,harms,kappa);
                phi=phi*tt;
                
                L(in,nw+(1:size(phi,2)))=phi;
                inds(nw+(1:length(harms)))=harms;
                ntip=ntip+1;
                nw=nw+size(phi,2);
                if  (zone{9}>0)
                    tips(iz,2)=ntip+1;
                    ztip=zo(cnodes(end,1));
                    tt=diff(zo(cnodes(end+(-1:0),1)));
                    tt=tt/abs(tt);
                    rtip=mean(abs(zo(indc{3})-ztip));
                    in=abs(zo-ztip)<=rtip;
                    in(indc{3})=1;
                    freenodes(in)=0;
                    
                    selected(in)=0;
                    selected(indc{3})=1;
                    
                    zzone=(zo(in)-ztip)*tt';
                    idin=1:prod(Nnodes);
                    idin=idin(in);
                    for ip=1:length(zzone)
                        if any(cnodes(:,1)==idin(ip))
                            zzone(ip)=abs(zzone(ip))*exp(-1i*pi);
                        end
                        if any(cnodes(:,2)==idin(ip))
                            zzone(ip)=abs(zzone(ip))*exp(1i*pi);
                        end
                    end
                    phi=Williams(zzone,modes,harms,kappa);
                    phi=phi*tt;
                    L(in,nw+(1:size(phi,2)))=phi;
                    inds(nw+(1:length(harms)))=harms;
                    nw=nw+size(phi,2);
                    ntip=ntip+1;
                end
                
            end
        end
    end
    nfree=sum(freenodes);
    indK=2*nfree+(1:size(L,2));
    Ln=sparse(find(freenodes),1:nfree,1,prod(Nnodes),nfree);
    L=[real(L);imag(L)];
    L=[blkdiag(Ln,Ln),L];
    %    selected(~freenodes)=0;
    
    if param.detect
        tnodes=[];
        for iz=1:size(model.zone,2)
            zone=model.zone(:,iz);
            if zone{4}==5
                if (zone{8}>0)
                    indc=zone{10};
                    cnodes=indc{1};
                    %                     htip=zone{7};
                    %                     rtip=zone{8};
                    
                    ztip=zo(cnodes(1,1));
                    rtip=min(abs(zo(indc{4})-ztip));
                    tnodes=abs(zo-ztip)<rtip;
                    
                    
                    %                     ntip=ceil(min(20,0.2*rtip)/htip);
                    %                     ttnodes=cnodes(1);
                    %                     for il=1:ntip
                    %                         ttnodes=AddOneNodeLayer(conn,ttnodes);
                    %                     end
                    %                         tnodes=[tnodes;ttnodes];
                    if  (zone{9}>0)
                        ztip=zo(cnodes(end,1));
                        rtip=min(abs(zo(indc{5})-ztip));
                        tnodes=tnodes|(abs(zo-ztip)<rtip);
                        %                                     rtip=mean(abs(zo(indc{3})-ztip));
                        %                 ttnodes=abs(zo-ztip)<=0.4*rtip;
                        
                        
                        %                         ttnodes=cnodes(end,1);
                        %                         for il=1:ntip
                        %                             ttnodes=AddOneNodeLayer(conn,ttnodes);
                        %                         end
                        %                         tnodes=[tnodes;ttnodes];
                    end
                end
            end
        end
        
        ielt=GetEltsFromNodes(conn,elt,tnodes,0);
        elt(ielt)=[];
        conn=conn(:,1:3);
        conn(ielt,:)=[];
        keep=unique(conn(:));
        xo=xo(keep)-roi(1)+1;
        yo=yo(keep)-roi(3)+1;
        selected=selected(keep);
        L=L([keep(:);keep+prod(Nnodes)],:);
        
        newids=zeros(prod(Nnodes),1);
        newids(keep)=1:length(keep);
        conn=newids(conn(:,1:3));
        
        elt=3*ones(size(conn,1),1);
        Nnodes=[length(xo),1,1];
        Nelems=[length(elt),1,1];
        conn=[conn,zeros(length(elt),1)];
        save(fullfile('TMP','0_mesh_0.mat'),'xo','yo','conn','elt','Nnodes','Nelems','-append')
    end
end
% param.regularization_type='none';
% param0.regularization_parameter=64;
% LoadParameters(param)
save(fullfile('TMP','0_mesh_0.mat'),'selected','-append')
%%
U1=[];Ks=[];
if ~isempty(Up),nscale=1;end
if param.thermo==1
    if pshift, error('NOT CODED YET');end
    for iscale=nscale:-1:1
        if isempty(Up)
            Uini=InitializeSolution(U1,iscale,nmod);
        else
            Uini=Up;
        end
        
        if (iscale==1)&&extract
            Pini=(L'*L)\(L'*Uini);
            %            Pini(inds<0,:)=0;
            Uini=L*Pini;
            [U1,dKs]=SolveIR2D(Uini,iscale,nmod,L);
            Ks=dKs(indK,:);
        else
            [U1]=SolveIR2D(Uini,iscale,nmod);
        end
    end
else
    if prod(sizeim)<4e6&&~param.normalize_grey_level&&(param.psample==1)
        
        for iscale=nscale:-1:1
            disp(sprintf('Pre-processing scale %d...',iscale));
            CreateBasisFunction(iscale,nmod);
            ComputeGradFPhi(iscale,nmod);
            CreateGradBasisFunction(iscale,nmod);
            AssembleCorrelationOperator(iscale,nmod);
            
            %%
            
            
            if isempty(Up)
                Uini=InitializeSolution(U1,iscale,nmod);
            else
                Uini=Up;
            end
            if (iscale==1)&&extract
                if pshift, error('NOT CODED YET');end
                Pini=(L'*L)\(L'*Uini);
                %            Pini(indK(inds<0),:)=0;
                %Pini(indK)'
                %keyboard
                %pause
                Uini=L*Pini;
                [U1,dKs]=RBMSolver(Uini,iscale,nmod,L);
                Ks=Pini(indK,:)+dKs(indK,:);
            else
                [U1]=Solve(Uini,iscale,nmod);
                
                if pshift&&iscale==1
                    Uini=U1;
                    for icam=2:pshift
                        if size(defs,2)==1
                            param.deformed_image=defs{icam,1};
                        else
                            param.deformed_image=defs(icam,:);
                        end
                        param.restart=1;
                        LoadParameters(param);
                        [Ui]=Solve(Uini,iscale,nmod);
                        Ux=mean(Uini((1:prod(Nnodes)))-Ui((1:prod(Nnodes))));
                        Uy=mean(Uini(prod(Nnodes)+(1:prod(Nnodes)))-Ui(prod(Nnodes)+(1:prod(Nnodes))));
                        Ui((1:prod(Nnodes)))=Ui((1:prod(Nnodes)))+Ux;
                        Ui(prod(Nnodes)+(1:prod(Nnodes)))=Ui(prod(Nnodes)+(1:prod(Nnodes)))+Uy;
                        U1=U1+Ui;
                    end
                    U1=U1/pshift;
                    param.iter_max=2;
                    if size(defs,2)==1
                        param.deformed_image=defs{1,1};
                    else
                        param.deformed_image=defs(1,:);
                    end
                    LoadParameters(param)
                    [~]=Solve(U1,1,nmod);
                    
                end
            end
            
        end
        
    else
        for iscale=nscale:-1:1
            if isempty(Up)
                Uini=InitializeSolution(U1,iscale,nmod);
            else
                Uini=Up;
            end
            
            if (iscale==1)&&extract
                if pshift, error('NOT CODED YET');end
                Pini=(L'*L)\(L'*Uini);
                %            Pini(inds<0,:)=0;
                Uini=L*Pini;
                [U1,dKs]=Solve2D(Uini,iscale,nmod,L);
                Ks=dKs(indK,:);
            else
                [U1]=Solve2D(Uini,iscale,nmod);
                if pshift&&iscale==1
                    Uini=U1;
                    for icam=2:pshift
                        if size(defs,2)==1
                            param.deformed_image=defs{icam,1};
                        else
                            param.deformed_image=defs(icam,:);
                        end
                        param.restart=1;
                        LoadParameters(param);
                        [Ui]=Solve2D(Uini,iscale,nmod);
                        Ux=mean(Uini((1:prod(Nnodes)))-Ui((1:prod(Nnodes))));
                        Uy=mean(Uini(prod(Nnodes)+(1:prod(Nnodes)))-Ui(prod(Nnodes)+(1:prod(Nnodes))));
                        Ui((1:prod(Nnodes)))=Ui((1:prod(Nnodes)))+Ux;
                        Ui(prod(Nnodes)+(1:prod(Nnodes)))=Ui(prod(Nnodes)+(1:prod(Nnodes)))+Uy;
                        U1=U1+Ui;
                    end
                    U1=U1/pshift;
                    param.iter_max=2;
                    if size(defs,2)==1
                        param.deformed_image=defs{1,1};
                    else
                        param.deformed_image=defs(1,:);
                    end
                    LoadParameters(param)
                    [~]=Solve2D(U1,1,nmod);
                    
                end
                
            end
        end
    end
end
%%
filreso=strrep(filres,'.res','');
U=U1;
copyfile(fullfile('TMP','0_error_0.mat'), [filreso,'-error.res']);
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
model=param;
load(fullfile('TMP','params'),'param');
load(fullfile('TMP','0_mesh_0'),'Nnodes','Nelems','xo','yo','conn','elt','ng','rflag','rint');
try
    tips(~(cell2mat(model.zone(4,:))==5),:)=[];
    model.zone(:,~(cell2mat(model.zone(4,:))==5))=[];
    cracks=find(cell2mat(model.zone(4,:))==5);
catch, end
%ExportImageToVTK(filres)
delete([fullfile('VTK',[filreso,'-error']),'-0*.vtk']);
if prod(sizeim)<4e6&~param.normalize_grey_level&&(param.thermo==0)&&(param.psample==1)
    if isfield(param,'deformed_image')
        if isfield(param,'image_number')
            images=param.image_number;
        else
            if iscell(param.deformed_image)
                images=1:size(param.deformed_image,2);
            else
                images=1;
            end
        end
    else
        images=(2:param.video_sampling:param.number_of_frames)-1;
    end
    for iim=1:size(U,2)
        movefile(fullfile('VTK',sprintf('camr-1-camd-1-scale-1-%d-error.vtk',iim)),fullfile('VTK',sprintf('%s-error-%04d.vtk',filreso,images(iim))));
    end
    VTKExportScalarMap([filreso,sprintf('-%d',0)],'error',param.roi(1)-1+(1:sizeim(1)),param.roi(3)-1+(1:sizeim(2)),zeros(sizeim),1);
    movefile(fullfile('VTK',[filreso,sprintf('-%d',0),'-error.vtk']),fullfile('VTK',sprintf('%s-error-%04d.vtk',filreso,0)));
end
delete(fullfile('VTK',['camr*','-error.vtk']));
save(filres,'U','Ks','tips','cracks','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
if isfield(param,'calibration_data')
    [U,Xo,Yo,Zo]=Extract3DFields(U1,nmod);
    CreateGradBasisFunction25D(iscale,nmod);
    
    save(filres,'U','U1','Xo','Yo','Zo','-append');
    writeVTKmesh25D(filres);
    postproVTK25D(filres,1);
else
    postproVTK([filres,''],0);
end
end