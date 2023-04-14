function postproVTKX3D(filres,submean)
display(sprintf('In postproVTKX3D\nResult file : %s',filres));
nmod=1;

tic
load([filres,'.mat'])
if isfield(param,'image_number')
    images=param.image_number;
else
images=1:size(U,2);
end
roi=param.roi;
if isfield(param,'calibration_data')
    roi=0*roi+1;
end
set.ascii=0;
set.remark=' computed by MIC';
nn=prod(Nnodes);
ne=length(elt);
if ~exist('ns'),ns=8*ones(1,2);end
images=[0,images];
U=[zeros(size(U,1),1),U];


nnew=sum(elt(face_elts)==8)*prod(ns)+sum(elt(face_elts)==6)*2*prod(ns);    
Ux=zeros(nnew*8,size(U,2));
Uy=zeros(nnew*8,size(U,2));
Uz=zeros(nnew*8,size(U,2));
xe=zeros(nnew*8,1);
ye=zeros(nnew*8,1);
ze=zeros(nnew*8,1);
xg=zeros(nnew,1);
yg=zeros(nnew,1);
zg=zeros(nnew,1);
nconn=reshape((1:8*nnew),8,nnew)';
inn=0;
ine=0;

lx=2/ns(1);ly=2/ns(2);lz=2/ns(3);

        [ygsub,xgsub,zgsub]=meshgrid((-1+ly/2):ly:(1-ly/2),(-1+lx/2):lx:(1-lx/2),(-1+lz/2):lz:(1-lz/2));
        xgsub=xgsub(:);ygsub=ygsub(:);zgsub=zgsub(:);
        xisub=lx/2*repmat([-1;1;1;-1],2,1);
        yisub=ly/2*repmat([-1;-1;1;1],2,1);
        zisub=lz/2*[-1;-1;-1;-1;1;1;1;1];
                xsub=repmat(xgsub',8,1);xsub=xsub(:)+repmat(xisub,prod(ns),1);
                ysub=repmat(ygsub',8,1);ysub=ysub(:)+repmat(yisub,prod(ns),1);
                zsub=repmat(zgsub',8,1);zsub=zsub(:)+repmat(zisub,prod(ns),1);

        nsub=length(xsub);
        nesub=length(xgsub);
          Ngq=[0.125*(1-xgsub).*(1-ygsub).*(1-zgsub),0.125*(1+xgsub).*(1-ygsub).*(1-zgsub),0.125*(1+xgsub).*(1+ygsub).*(1-zgsub),0.125*(1-xgsub).*(1+ygsub).*(1-zgsub),...
            0.125*(1-xgsub).*(1-ygsub).*(1+zgsub),0.125*(1+xgsub).*(1-ygsub).*(1+zgsub),0.125*(1+xgsub).*(1+ygsub).*(1+zgsub),0.125*(1-xgsub).*(1+ygsub).*(1+zgsub)];
          Nnq=[0.125*(1-xsub).*(1-ysub).*(1-zsub),0.125*(1+xsub).*(1-ysub).*(1-zsub),0.125*(1+xsub).*(1+ysub).*(1-zsub),0.125*(1-xsub).*(1+ysub).*(1-zsub),...
            0.125*(1-xsub).*(1-ysub).*(1+zsub),0.125*(1+xsub).*(1-ysub).*(1+zsub),0.125*(1+xsub).*(1+ysub).*(1+zsub),0.125*(1-xsub).*(1+ysub).*(1+zsub)];
      
%          figure
%         plot3(xsub,ysub,zsub,'bx')
%         hold on
%         plot3(xgsub,ygsub,zgsub,'r+')
for ie=1:length(face_elts)
    ielt=face_elts(ie);
    inods=conn(ielt,1:elt(ielt));
    xn=xo(inods);yn=yo(inods);zn=zo(inods);
    if elt(ielt)==8
        xe(inn+(1:nsub))=Nnq*xn;
        ye(inn+(1:nsub))=Nnq*yn;
        ze(inn+(1:nsub))=Nnq*zn;
        inn=inn+nsub;
        xg(ine+(1:nesub))=Ngq*xn;
        yg(ine+(1:nesub))=Ngq*yn;
        zg(ine+(1:nesub))=Ngq*zn;
        ine=ine+nesub;
    else
        error('not coded yet');

    end
end
load(model.levelset_file,'crack','front','zone');
xn=xo(face_nodes)+roi(1)-zone(1);yn=yo(face_nodes)+roi(3)-zone(3);zn=zo(face_nodes)+roi(5)-zone(5);
 hn=mexInterpLinear3D(xn,yn,zn,crack);
hn=double(hn>=0);
hs=mexInterpLinear3D(xg+roi(1)-zone(1),yg+roi(3)-zone(3),zg+roi(5)-zone(5),crack);
hs=double(hs>=0);
%hs=repmat(hs,1,3);
inn=0;
ine=0;
non=sum((crack(:)<=0)&(crack(:)>-1)&(front(:)<=0));
xf=zeros(non,1);
yf=zeros(non,1);
zf=zeros(non,1);
ff=zeros(non,1);
Dx=zeros(non,size(U,2));
Dy=zeros(non,size(U,2));
Dz=zeros(non,size(U,2));
ion=0;
for ie=1:length(face_elts)
    ielt=face_elts(ie);
    inods=conn(ielt,1:elt(ielt));
    if elt(ielt)==8
    hG=repmat(hs(ine+(1:nesub))',8,1);
    Uxc=Nnq*U(inods,:);
    Uyc=Nnq*U(inods+nn,:);
    Uzc=Nnq*U(inods+2*nn,:);
Uxx=0;Uyx=0;Uzx=0;
for in=1:8
   ienr=find(face_nodes==inods(in));
   if ~isempty(ienr)
       H=diag(sparse(hG(:)-hn(ienr)));
   Uxx=Uxx+H*Nnq(:,in)*U(3*nn+ienr,:);
   Uyx=Uyx+H*Nnq(:,in)*U(3*nn+ienr+length(face_nodes),:);
   Uzx=Uzx+H*Nnq(:,in)*U(3*nn+ienr+2*length(face_nodes),:);
   end
end
    Ux(inn+(1:nsub),:)=Uxc+Uxx;
    Uy(inn+(1:nsub),:)=Uyc+Uyx;
    Uz(inn+(1:nsub),:)=Uzc+Uzx;
            inn=inn+nsub;
        ine=ine+nesub;
    else
        error('not coded yet');
    end
    xn=xo(inods);yn=yo(inods);zn=zo(inods);
    xpix=ceil(min(xn)):floor(max(xn));
    ypix=ceil(min(yn)):floor(max(yn));
    zpix=ceil(min(zn)):floor(max(zn));
    cracke=crack(xpix+roi(1)-zone(1),ypix+roi(3)-zone(3),zpix+roi(5)-zone(5));
    fronte=front(xpix+roi(1)-zone(1),ypix+roi(3)-zone(3),zpix+roi(5)-zone(5));
    
    oncrack=find((cracke<=0)&(cracke>-1)&(fronte<=0));
    [ypix,xpix,zpix]=meshgrid(ypix,xpix,zpix);
    xpix=xpix(oncrack);ypix=ypix(oncrack);zpix=zpix(oncrack);
            if elt(ielt)==6
                [xge,yge,zge,wge]=GetGaussPointsWedge(0,[1,1,1],xn,yn,zn,xpix(:),ypix(:),zpix(:));
                          N=[0.5*(1-xge-yge).*(1-zge),0.5*(xge).*(1-zge),0.5*(yge).*(1-zge),...
                    0.5*(1-xge-yge).*(1+zge),0.5*(xge).*(1+zge),0.5*(yge).*(1+zge)];
  elseif elt(ielt)==8
                [xge,yge,zge,wge]=GetGaussPointsHexaedron(0,[1,1,1],xn,yn,zn,xpix(:),ypix(:),zpix(:));
                N=[0.125*(1-xge).*(1-yge).*(1-zge),0.125*(1+xge).*(1-yge).*(1-zge),0.125*(1+xge).*(1+yge).*(1-zge),0.125*(1-xge).*(1+yge).*(1-zge),...
                    0.125*(1-xge).*(1-yge).*(1+zge),0.125*(1+xge).*(1-yge).*(1+zge),0.125*(1+xge).*(1+yge).*(1+zge),0.125*(1-xge).*(1+yge).*(1+zge)];
            end
            indx=ion+(1:length(oncrack));
    xf(indx)=xpix;
    yf(indx)=ypix;
    zf(indx)=zpix;
    ff(indx)=-fronte(oncrack);
Uxx=0;Uyx=0;Uzx=0;
for in=1:8
   ienr=find(face_nodes==inods(in));
   if ~isempty(ienr)
   Uxx=Uxx+N(:,in)*U(3*nn+ienr,:);
   Uyx=Uyx+N(:,in)*U(3*nn+ienr+length(face_nodes),:);
   Uzx=Uzx+N(:,in)*U(3*nn+ienr+2*length(face_nodes),:);
   end
end
    Dx(indx,:)=Uxx;
    Dy(indx,:)=Uyx;
    Dz(indx,:)=Uzx;
    ion=ion+length(oncrack);
    
    
end
if ion<non
   xf(ion+1:non)=[];
   yf(ion+1:non)=[];
   zf(ion+1:non)=[];
   ff(ion+1:non)=[];
   Dx(ion+1:non,:)=[];
   Dy(ion+1:non,:)=[];
   Dz(ion+1:non,:)=[];
   non=ion;
end

unix(['rm ',fullfile('VTK',filres),'-x-0*.vtk']);
L=[1+0*xo,0*xo,0*xo,0*xo,-zo,yo;...
    0*yo,1+0*yo,0*yo,zo,0*yo,-xo;...
    0*zo,0*zo,1+0*zo,-yo,xo,0*zo];
Lxs=[1+0*xe,0*xe,0*xe,0*xe,-ze,ye];
Lys=[0*ye,1+0*ye,0*ye,ze,0*ye,-xe];
Lzs=[0*ze,0*ze,1+0*ze,-ye,xe,0*ze];

for iim=1:size(U,2)
    set.vtkname=[filres , sprintf('-x-%06d.vtk',images(iim))];
    %%
    %fwid = fopen(set.vtkname,'w','b'); % IMPORTANT: big endian
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
    count = fprintf(fwid,'POINTS %u float\n',8*nnew);
    data=[(xe-1)+roi(1),(ye-1)+roi(3),(ze-1)+roi(5)]';


    if set.ascii
        fprintf(fwid, '%f %f %f \n', data);
    else
        fwrite(fwid, data,'float');
    end
    count = fprintf(fwid,'CELLS %u %u\n',nnew,9*nnew);
    data=[repmat(8,1,nnew);nconn'-1];


    if set.ascii
        fprintf(fwid, '%d %d %d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end


    count = fprintf(fwid,'CELL_TYPES %u\n',nnew);
    data=repmat(12,1,nnew);

    if set.ascii
        fprintf(fwid, '%d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    count = fprintf(fwid,'POINT_DATA %u\n',8*nnew);
    count = fprintf(fwid,['VECTORS Displacement float\n']);
    Ui=U(1:(3*nn),iim);
    Uxi=Ux(:,iim);
    Uyi=Uy(:,iim);
    Uzi=Uz(:,iim);
    if submean
    A=L\Ui;
        Uxi=Uxi-Lxs*A;
        Uyi=Uyi-Lys*A;
        Uzi=Uzi-Lzs*A;
    end
%     figure
%     trimesh(nconn,xe,ye,0*ye,Uxi)
%     axis image
    UU=[Uxi';Uyi';Uzi'];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end

    fclose(fwid);
end

unix(['rm ',fullfile('VTK',filres),'-xon-0*.vtk']);
for iim=1:size(U,2)
    set.vtkname=[filres , sprintf('-xon-%06d.vtk',images(iim))];
    %%
    %fwid = fopen(set.vtkname,'w','b'); % IMPORTANT: big endian
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
    count = fprintf(fwid,'POINTS %u float\n',non);
    data=[(xf-1)+roi(1),(yf-1)+roi(3),(zf-1)+roi(5)]';


    if set.ascii
        fprintf(fwid, '%f %f %f \n', data);
    else
        fwrite(fwid, data,'float');
    end
    count = fprintf(fwid,'CELLS %u %u\n',non,2*non);
    data=[repmat(1,1,non);(1:non)-1];


    if set.ascii
        fprintf(fwid, '%d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end


    count = fprintf(fwid,'CELL_TYPES %u\n',non);
    data=repmat(2,1,non);

    if set.ascii
        fprintf(fwid, '%d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    count = fprintf(fwid,'POINT_DATA %u\n',non);
    count = fprintf(fwid,['VECTORS DU float\n']);
    Uxi=Dx(:,iim);
    Uyi=Dy(:,iim);
    Uzi=Dz(:,iim);
    UU=[Uxi';Uyi';Uzi'];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
            count = fprintf(fwid,['SCALARS Front float, 1\n']);
            count = fprintf(fwid,'LOOKUP_TABLE default\n');
            if set.ascii
                fprintf(fwid, '%f\n', ff);
            else
                fwrite(fwid,ff,'float');
            end

    fclose(fwid);
end






end