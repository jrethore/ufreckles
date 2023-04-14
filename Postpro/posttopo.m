function posttopo(filreso,comp)
filres=filreso;
if nargin<2, comp=false;end
disp(['Result file:',filreso]);
lim=[];
epsflag=false;
nmod=1;
load([filreso,'.mat'])
filexp='';
if ~exist('model1')
model1=model;
end
roi=param.roi;
if comp
    filreso=['comp-',filreso];
end
if isfield(param,'pixel_size')
    pix2m=param.pixel_size;
else 
    pix2m=1;
end
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
if isfield(model1,'crack_id')
    ncrack=model1.crack_id;
else
   ncrack=1; 
end
if strcmp(param.analysis,'mechanics')
    fac=pix2m;
else
    fac=1;
end
submean=0;

unstruc=0;
if isfield(model1,'mesh_file')
    unstruc=1;
end

%% mask

load(fullfile('TMP',[num2str(nmod),'_mask_0']),'mask');
mask=double(diag(mask)>0);
found=find(mask==0);
mask=full(mask);
mask(found)=NaN;
if unstruc
    xn=xo;
    yn=yo;
xo=[ceil(min(xo));floor(max(xo))];
yo=[ceil(min(yo));floor(max(yo))];
else
xo=reshape(xo,Nnodes)+0.5;
yo=reshape(yo,Nnodes)+0.5;
xn=xo(1:Nnodes(1),1)';
yn=yo(1,1:Nnodes(2));
% xn=(xo);
% xn(length(xn))=xo(length(xo))-1;
% yn=(yo);
% yn(length(yn))=yo(length(yo))-1;
end



    load(fullfile('TMP','sample0'),'sizeim');
    sizeim0=sizeim;
load(fullfile('TMP',[num2str(nmod),'_phi_0']),'phi','sizeim','Xi','Yi');
if unstruc
   hasval=(sum(phix,2)>0)|(sum(phiy,2)>0);
   found=find(~hasval);
   if numel(mask)==1
       mask=repmat(1,[prod(sizeim),1]);
   end
   mask(found)=NaN;
end
%% displacement map
umap=fac*reshape((phi*U).*mask,sizeim);
umap=umap(xn(1):xn(length(xn))-1,yn(1):yn(length(yn))-1);
%[limx]=FindMinMaxDisp(min(umap(:)),max(umap(:)));
if comp
    limx=[-10,10];
else
    limx=[];
end
PlotMap(umap,'Z',limx,[filres,'-u'],epsflag);

[ys,xs]=meshgrid((yn-1)*psample+roi(3)-0.5*sizeim0(2),(xn-1)*psample+roi(1)-0.5*sizeim0(1));
C=param.calibration_data{1};
%if isempty(unmasked_nodes)
    Zw=reshape(U,Nnodes);
%
    [Xw,Yw]=GetXYFromabZ(C,xs(:),ys(:),Zw(:));
    Xw=reshape(Xw,size(Zw));
    Yw=reshape(Yw,size(Zw));
    
DeformedMesh25D(Xw,Yw,Zw,Zw,lim,1,[filres,'-topo'],epsflag,1);

%% error map
nim=length(param.reference_image)-1;
    if exist([filreso,'-error.mat'],'file')
        fide=fopen([filreso,'-error.mat']);
    else
        fide=fopen(fullfile('TMP',[num2str(nmod),'_error_0.mat']));
    end
    dynamic=fread(fide,1);

for iim=1:nim
        if nim>1
        filexp=sprintf('-%03d',iim+1);
        else
            filexp='';
    end
    filres=[filreso,filexp];

    disc=fread(fide,prod(sizeim));
    disc=reshape(disc,sizeim);
    disc=100*abs(disc(xn(1):xn(length(xn))-1,yn(1):yn(length(yn))-1))/dynamic;
    if comp
        limd=[0,10];
    else
    limd=[0,0.5*100*max(disc(:))/dynamic];
    end
   PlotMap(disc,'Discrepancy wrt dyn (%)',limd,[filres,'-error'],epsflag,0);
    disp(sprintf('Mean error wrt image dyn.: %6.2f %%',mean(disc(:))));
end
    fclose(fide);

%% ROI
    load(fullfile('TMP','sample0'),'im0')
    PlotZOI(roi,(xo-1)*psample+roi(1),(yo-1)*psample+roi(3),im0,'ZOI',[filres,'-zoi'],epsflag);
    clear im0

end


