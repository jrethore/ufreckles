function postproVTKX(filres,submean)
display(sprintf('In postproVTKX\nResult file : %s',filres));
nmod=1;
tic
load([filres,'.mat'])
if ~exist('U','var')
    U=U1;
end
if ~exist('model1')
    model1=model;
end
if ~exist('unmasked_nodes','var')
    unmasked_nodes=[];
end
if ~isempty(unmasked_nodes)
    Nnodes=[length(unmasked_nodes),1,1];
    idmask=zeros(prod(Nnodes),1);
    idmask(unmasked_nodes)=1:length(unmasked_nodes);
else
    idmask=1:length(xo);
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
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
roi=param.roi;
images=[0,images];
U=[zeros(size(U,1),1),U];
%%
set.ascii=0;
set.remark=' computed by MIC';
nn=prod(Nnodes);
ne=length(elt);
if ~exist('ns'),ns=8*ones(1,2);end
nnew4=sum(elt(face_elts)==4)*prod(ns);
nnew3=sum(elt(face_elts)==3)*2*prod(ns);
Ux4=zeros(nnew4*4,size(U,2));
Uy4=zeros(nnew4*4,size(U,2));
Ux3=zeros(nnew3*3,size(U,2));
Uy3=zeros(nnew3*3,size(U,2));
xe3=zeros(nnew3*3,1);
ye3=zeros(nnew3*3,1);
xg3=zeros(nnew3,1);
yg3=zeros(nnew3,1);
xe4=zeros(nnew4*4,1);
ye4=zeros(nnew4*4,1);
xg4=zeros(nnew4,1);
yg4=zeros(nnew4,1);
nconn3=reshape((1:(3*nnew3)),3,nnew3)';
nconn4=reshape((1:(4*nnew4)),4,nnew4)';

lx=2/ns(1);ly=2/ns(2);

[ycsub,xcsub]=meshgrid((-1+ly/2):ly:(1-ly/2),(-1+lx/2):lx:(1-lx/2));
xisub3=[-lx/2;0;-lx/2;...
    -lx/2;0;lx/2;
    lx/2;0;lx/2;
    lx/2;0;-lx/2];
yisub3=[ ly/2;0;-ly/2;...
    -ly/2;0;-ly/2;
    -ly/2;0;ly/2;
    ly/2;0;ly/2];
xcsub3=repmat(xcsub(:)',4*3,1);
ycsub3=repmat(ycsub(:)',4*3,1);
xsub3=repmat(xisub3,prod(ns),1)+xcsub3(:);
ysub3=repmat(yisub3,prod(ns),1)+ycsub3(:);
xsub3=reshape(xsub3,3,4*prod(ns))';
ysub3=reshape(ysub3,3,4*prod(ns))';
xgsub3=mean(xsub3,2);
ygsub3=mean(ysub3,2);
%           figure
%         plot(xsub3,ysub3,'bx')
%         hold on
%         plot(xcsub3,ycsub3,'ko')
%         plot(xgsub3,ygsub3,'r+')

found=find(ygsub3<=-xgsub3);
ygsub3=(ygsub3(found)+1)/2;
xgsub3=(xgsub3(found)+1)/2;
nesub3=length(xgsub3);
ysub3=(ysub3(found,:)'+1)/2;
xsub3=(xsub3(found,:)'+1)/2;
xsub3=xsub3(:);ysub3=ysub3(:);
nsub3=length(xsub3);
%           figure
%         plot(xsub3,ysub3,'bx')
%         hold on
%         plot(xgsub3,ygsub3,'r+')

xisub4=[-lx/2;lx/2;...
    lx/2;-lx/2];
yisub4=[ -ly/2;-ly/2;...
    ly/2;ly/2];
xcsub4=repmat(xcsub(:)',4,1);
ycsub4=repmat(ycsub(:)',4,1);
xsub4=repmat(xisub4,prod(ns),1)+xcsub4(:);
ysub4=repmat(yisub4,prod(ns),1)+ycsub4(:);
nsub4=length(xsub4);
xgsub4=mean(reshape(xsub4,4,prod(ns))',2);
ygsub4=mean(reshape(ysub4,4,prod(ns))',2);
nesub4=length(xgsub4);
Nnq=[0.25*(1-xsub4).*(1-ysub4),0.25*(1+xsub4).*(1-ysub4),0.25*(1+xsub4).*(1+ysub4),0.25*(1-xsub4).*(1+ysub4)];
Ngq=[0.25*(1-xgsub4).*(1-ygsub4),0.25*(1+xgsub4).*(1-ygsub4),0.25*(1+xgsub4).*(1+ygsub4),0.25*(1-xgsub4).*(1+ygsub4)];
Nnt=[1-xsub3-ysub3,xsub3,ysub3];
Ngt=[1-xgsub3-ygsub3,xgsub3,ygsub3];
%         figure
%         plot(xsub4,ysub4,'bx')
%         hold on
%         plot(xcsub4,ycsub4,'ko')
%         plot(xgsub4,ygsub4,'r+')
inn3=0;ine3=0;inn4=0;ine4=0;
for ie=1:length(face_elts)
    ielt=face_elts(ie);
    inods=conn(ielt,1:elt(ielt));
    xn=xo(inods);yn=yo(inods);
    if elt(ielt)==4
        xe4(inn4+(1:nsub4))=Nnq*xn;
        ye4(inn4+(1:nsub4))=Nnq*yn;
        inn4=inn4+nsub4;
        xg4(ine4+(1:nesub4))=Ngq*xn;
        yg4(ine4+(1:nesub4))=Ngq*yn;
        ine4=ine4+nesub4;
    else
        xe3(inn3+(1:nsub3))=Nnt*xn;
        ye3(inn3+(1:nsub3))=Nnt*yn;
        inn3=inn3+nsub3;
        xg3(ine3+(1:nesub3))=Ngt*xn;
        yg3(ine3+(1:nesub3))=Ngt*yn;
        ine3=ine3+nesub3;
    end
end
% figure
% triplot(nconn,xe,ye)
load(fullfile('TMP',sprintf('%d_levelsets',1)),'crack','front');
hn=interp2(crack,yo,xo,'*linear',0);
hn=double(hn>=0);
hs3=interp2(crack,yg3,xg3,'*linear',0);
hs3=double(hs3>=0);
hs4=interp2(crack,yg4,xg4,'*linear',0);
hs4=double(hs4>=0);
%hs=repmat(hs,1,3);
inn3=0;ine3=0;inn4=0;ine4=0;
lcz=0;
if isfield(param,'cz_length')
    lcz=param.cz_length;
end
non=sum((crack(:)<=0)&(crack(:)>-1)&(front(:)<=lcz));
xf=zeros(non,1);
yf=zeros(non,1);
ff=zeros(non,1);
Dx=zeros(non,size(U,2));
Dy=zeros(non,size(U,2));
ion=0;
for ie=1:length(face_elts)
    ielt=face_elts(ie);
    inods=conn(ielt,1:elt(ielt));
    hN=hn(inods);
    if elt(ielt)==4
        hG=repmat(hs4(ine4+(1:nesub4))',4,1);
        Uxc=Nnq*U(idmask(inods),:);
        Uyc=Nnq*U(idmask(inods)+prod(Nnodes),:);
        Uxx=0;Uyx=0;
        for in=1:4
            ienr=find(face_nodes==inods(in));
            if ~isempty(ienr)
                H=diag(sparse(hG(:)-hN(in)));
                Uxx=Uxx+H*Nnq(:,in)*U(2*prod(Nnodes)+ienr,:);
                Uyx=Uyx+H*Nnq(:,in)*U(2*prod(Nnodes)+ienr+length(face_nodes),:);
            end
        end
        Ux4(inn4+(1:nsub4),:)=Uxc+Uxx;
        Uy4(inn4+(1:nsub4),:)=Uyc+Uyx;
        inn4=inn4+nsub4;
        ine4=ine4+nesub4;
    else
        hG=repmat(hs3(ine3+(1:nesub3))',3,1);
        Uxc=Nnt*U(idmask(inods),:);
        Uyc=Nnt*U(idmask(inods)+prod(Nnodes),:);
        Uxx=0;Uyx=0;
        for in=1:3
            ienr=find(face_nodes==inods(in));
            if ~isempty(ienr)
                H=diag(sparse(hG(:)-hN(in)));
                Uxx=Uxx+H*Nnt(:,in)*U(2*prod(Nnodes)+ienr,:);
                Uyx=Uyx+H*Nnt(:,in)*U(2*prod(Nnodes)+ienr+length(face_nodes),:);
            end
        end
        Ux3(inn3+(1:nsub3),:)=Uxc+Uxx;
        Uy3(inn3+(1:nsub3),:)=Uyc+Uyx;
        inn3=inn3+nsub3;
        ine3=ine3+nesub3;
    end
    xn=xo(inods);yn=yo(inods);
    xpix=ceil(min(xn)):floor(max(xn));
    ypix=ceil(min(yn)):floor(max(yn));
    cracke=crack(xpix,ypix);
    fronte=front(xpix,ypix);

    oncrack=find((cracke<=0)&(cracke>-1)&(fronte<=lcz));
    [ypix,xpix]=meshgrid(ypix,xpix);
    xpix=xpix(oncrack);ypix=ypix(oncrack);
    Uxx=0;Uyx=0;
    if elt(ielt)==3
        [xge,yge,wge]=GetGaussPointsTriangle(0,[1,1,1],xn,yn,xpix(:),ypix(:));
        ind=find(wge);
        xge=xge(ind);yge=yge(ind);
        oncrack=oncrack(ind);
        N=[1-xge-yge,xge,yge];
    elseif elt(ielt)==4
        [xge,yge,wge]=GetGaussPointsQuadrangle(0,[1,1,1],xn,yn,xpix(:),ypix(:));
        N=[0.25*(1-xge).*(1-yge),0.25*(1+xge).*(1-yge),0.25*(1+xge).*(1+yge),0.25*(1-xge).*(1+yge)];
    end
    indx=ion+(1:length(oncrack));
    xf(indx)=N*xn;
    yf(indx)=N*yn;
    ff(indx)=-fronte(oncrack);
    for in=1:elt(ielt)
        ienr=find(face_nodes==inods(in));
        if ~isempty(ienr)
            Uxx=Uxx+N(:,in)*U(2*prod(Nnodes)+ienr,:);
            Uyx=Uyx+N(:,in)*U(2*prod(Nnodes)+ienr+length(face_nodes),:);
        end

    end

    Dx(indx,:)=Uxx;
    Dy(indx,:)=Uyx;
    ion=ion+length(oncrack);


end
if ion<non
    xf=xf(1:ion);
    yf=yf(1:ion);
    ff=ff(1:ion);
    Dx=Dx(1:ion,:);
    Dy=Dy(1:ion,:);
    non=ion;
end
non=length(xf);
unix(['rm ',fullfile('VTK',filres),'-x-0*.vtk']);
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
    count = fprintf(fwid,'POINTS %u float\n',3*nnew3+4*nnew4);
    data=[(xe3-1)*psample+roi(1),(ye3-1)*psample+roi(3),0*ye3]';


    if set.ascii
        fprintf(fwid, '%f %f %f \n', dfac*data);
    else
        fwrite(fwid, data,'float');
    end
    data=[(xe4-1)*psample+roi(1),(ye4-1)*psample+roi(3),0*ye4]';


    if set.ascii
        fprintf(fwid, '%f %f %f \n', dfac*data);
    else
        fwrite(fwid, data,'float');
    end
    count = fprintf(fwid,'CELLS %u %u\n',nnew3+nnew4,4*nnew3+5*nnew4);
    data=[repmat(3,1,nnew3);nconn3'-1];


    if set.ascii
        fprintf(fwid, '%d %d %d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    data=[repmat(4,1,nnew4);nconn4'-1];


    if set.ascii
        fprintf(fwid, '%d %d %d %d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end


    count = fprintf(fwid,'CELL_TYPES %u\n',nnew3+nnew4);
    data=repmat(5,1,nnew3);

    if set.ascii
        fprintf(fwid, '%d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    data=repmat(9,1,nnew4);

    if set.ascii
        fprintf(fwid, '%d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    count = fprintf(fwid,'POINT_DATA %u\n',3*nnew3+4*nnew4);
    count = fprintf(fwid,['VECTORS Displacement float\n']);
    Ui=U(1:(2*length(xo)),iim);
    Uxi=Ux3(:,iim);
    Uyi=Uy3(:,iim);
    if submean
        Uxi=Uxi-mean(Ui(1:nn));
        Uyi=Uyi-mean(Ui(nn+(1:nn)));
    end
    %     figure
    %     trimesh(nconn,xe,ye,0*ye,Uxi)
    %     axis image
    UU=fac*[Uxi';Uyi';repmat(0,1,3*nnew3)];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
    Uxi=Ux4(:,iim);
    Uyi=Uy4(:,iim);
    if submean
        Uxi=Uxi-mean(Ui(1:nn));
        Uyi=Uyi-mean(Ui(nn+(1:nn)));
    end
    UU=fac*[Uxi';Uyi';repmat(0,1,4*nnew4)];
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
    fwid = fopen(fullfile('VTK',set.vtkname),'w'); % IMPORTANT: big endian
    count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
    count = fprintf(fwid,[set.vtkname,set.remark,'\n']);
    if set.ascii
        count = fprintf(fwid,'ASCII\n');
    else
        count = fprintf(fwid,'BINARY\n');
    end
    count = fprintf(fwid,'DATASET UNSTRUCTURED_GRID\n');
    count = fprintf(fwid,'POINTS %u float\n',non);
    data=[(xf-1)*psample+roi(1),(yf-1)*psample+roi(3),0*yf]';


    if set.ascii
        fprintf(fwid, '%f %f %f \n', dfac*data);
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
    count = fprintf(fwid,['VECTORS D float\n']);
    Uxi=Dx(:,iim);
    Uyi=Dy(:,iim);
    UU=fac*[Uxi';Uyi';repmat(0,1,non)];
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