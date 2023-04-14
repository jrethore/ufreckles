function [Uf]=RefineMesh(nmod,fac,Uo)
Uf=[];


load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
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
if isfield(param0,'deformed_image')
    if iscell(param0.deformed_image)
        nim=size(param0.deformed_image,2);
        ncamd=size(param0.deformed_image,1);
    else
        nim=1;
        ncamd=1;
    end
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
for icamr=1:ncamr
    iscale=1;
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'elt','conn','xo','yo','Nnodes','selected');
    selected=selected(:);
    if isfield(param,'topography')
     load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'Xo','Yo','Zo');
          [elt2,conn2,Xo,Yo,Nnodes2,Nelems2,selected2,Zo]=BuiltRefinedConnectivity(elt,conn,Xo,Yo,selected,fac);
     
     save(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'Xo','Yo','Zo','-append');
        
    end

    if ~strcmp(param.basis,'nurbs')&&~isfield(param,'enrichment')&&nargin>2
        [elt,conn,xo,yo,Nnodes,Nelems,selected,zo,Uf]=BuiltRefinedConnectivity(elt,conn,xo,yo,selected,fac,0*xo,Uo);
    else
        [elt,conn,xo,yo,Nnodes,Nelems,selected]=BuiltRefinedConnectivity(elt,conn,xo,yo,selected,fac);
    end
    save(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'elt','conn','xo','yo','Nnodes','Nelems','selected','-append');
    if isfield(param,'enrichment')
        save(fullfile('TMP',sprintf('%d_meshx_%d',nmod*10^(icamr-1),iscale-1)),'elt','conn','xo','yo','Nnodes','Nelems','selected','-append');
    end
    if (strcmp(param.basis,'nurbs'))&&nargin>2
        if isempty(Uf)
            Uf=repmat(0,[2*length(xo),size(Uo,2),size(Uo,3)]);
        end

        load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),10*(1-1))),'phix','sizeim');
        load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icamr-1),10*(1-1))),'phiy');
        phixo=phix;phiyo=phiy;
        if ~strcmp(param.basis,'nurbs')
            CreateBasisFunction(1,nmod);
            load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),10*(1-1))),'phix','sizeim');
            load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icamr-1),10*(1-1))),'phiy');
            M=phix'*phix+phiy'*phiy;
        end
        for ijm=1:nim
            for icamd=1:size(indcam,2)
                iz=indcam(icamr,icamd)+(icamr-1)*size(indcam,2);
                if strcmp(param.basis,'nurbs')
                    Ux=interp2(reshape(phixo*Uo(:,ijm,iz),sizeim),min(max(ceil(min(yo+0.)),yo),floor(max(yo-0.))),min(max(ceil(min(xo+0.)),xo),floor(max(xo-0.))),'*linear');
                    Uy=interp2(reshape(phiyo*Uo(:,ijm,iz),sizeim),min(max(ceil(min(yo+0.)),yo),floor(max(yo-0.))),min(max(ceil(min(xo+0.)),xo),floor(max(xo-0.))),'*linear');
                    Uf(:,ijm,iz)=[Ux;Uy];
                else
                    Ui=M\(phix'*phixo*Uo(:,ijm,iz)+phiy'*phiyo*Uo(:,ijm,iz));
                    Uf(:,ijm,iz)=Ui;
                end
            end
        end
    end
end