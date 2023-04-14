function postproVTKW(filres,submean)
if nargin<2,submean=0;end
display(sprintf('In postproVTKW\nResult file : %s',filres));
nmod=1;
tic
load([filres,'.mat'])
if ~exist('U','var')
    U=U2;
    nmod=2;
else
    nmod=1;
end
if ~exist('model2')
    model2=model;
end
if ~exist('unmasked_nodes','var')
    unmasked_nodes=[];
end
if ~isempty(unmasked_nodes)
    Nnodes=[length(unmasked_nodes),1,1];
end
if isfield(param,'image_number')
    images=param.image_number;
else
    images=1:size(U,2);
end
if isfield(param,'pixel_size')
    pix2m=param.pixel_size;
else
    pix2m=1;
end
if strcmp(param.analysis,'mechanics')
    dfac=1;
    fac=1+0*pix2m;
else
    fac=1;dfac=1;
end
roi=param.roi;
%%
set.ascii=0;
set.remark=' computed by MIC';
load(fullfile('TMP','sample0_0'),'sizeim');
load(fullfile('TMP',[num2str(nmod),'_mask_0']),'mask');
mask=double(diag(mask));
found=find((mask(:)==1));
[yi,xi]=meshgrid(1:sizeim(2),1:sizeim(1));
yi=yi(found);
xi=xi(found);
images=[0,images];
U=[zeros(size(U,1),1),U];

load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix','sizeim');
load(fullfile('TMP',[num2str(nmod),'_phiy_0']),'phiy','sizeim');
Ux=phix*U;
Uy=phiy*U;
Ux=Ux(found,:);
Uy=Uy(found,:);

ind=model.km_indices;
Uo=U;
foundi=find(ind==0);
Uo(foundi,:)=0;
Uo(foundi+length(ind),:)=0;
foundi=find(ind==2);
Uo(foundi+length(ind),:)=0;
Uxo=phix*Uo;
Uyo=phiy*Uo;
Uxo=Uxo(found,:);
Uyo=Uyo(found,:);

nx=length(found);
nn=prod(Nnodes);

for iim=1:size(U,2)
    set.vtkname=[filres , sprintf('-w-%06d.vtk',images(iim))];
    %%
    unix(['rm ',fullfile('VTK',set.vtkname)]);
    if set.ascii
    fwid = fopen(fullfile('VTK',set.vtkname),'w'); % IMPORTANT: big endian
    else
        fwid = fopen(fullfile('VTK',set.vtkname),'w','b'); % IMPORTANT: big endian
    end
    count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
    count = fprintf(fwid,[set.vtkname,set.remark,'\n']);
    if set.ascii
        count = fprintf(fwid,'ASCII\n');
    else
        count = fprintf(fwid,'BINARY\n');
    end
    count = fprintf(fwid,'DATASET UNSTRUCTURED_GRID\n');
    count = fprintf(fwid,'POINTS %u float\n',nx);
    data=[(xi-1)+roi(1),(yi-1)+roi(3),0*yi]';


    if set.ascii
        fprintf(fwid, '%f %f %f \n', dfac*data);
    else
        fwrite(fwid, data,'float');
    end
    count = fprintf(fwid,'CELLS %u %u\n',nx,2*nx);
    data=[repmat(1,1,nx);(1:nx)-1];


    if set.ascii
        fprintf(fwid, '%d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end


    count = fprintf(fwid,'CELL_TYPES %u\n',nx);
    data=repmat(2,1,nx);

    if set.ascii
        fprintf(fwid, '%d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    count = fprintf(fwid,'POINT_DATA %u\n',nx);
    count = fprintf(fwid,['VECTORS Displacement float\n']);
    Uxi=Ux(:,iim);
    Uyi=Uy(:,iim);
    UU=fac*[Uxi';Uyi';repmat(0,1,nx)];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
    count = fprintf(fwid,['VECTORS U float\n']);
    Uxi=Uxo(:,iim);
    Uyi=Uyo(:,iim);
    UU=fac*[Uxi';Uyi';repmat(0,1,nx)];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end

    fclose(fwid);
end
end