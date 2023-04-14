function ExportUnwarped3DImage(filres)
display(sprintf('In ExportUnwarped3DImage\nResult file : %s',filres));
nmod=1;
iscale=1;
tic
load([filres,'.mat'])
if ~exist('zone','var')
    zone=ones(1,6);
end
if isfield(param,'image_number')
    images=param.image_number;
else
images=1:size(U,2);
end
roi=param.roi;
set.ascii=0;
set.remark=' computed by MIC';

unix(['rm ',fullfile('VTK',filres),'-images-0*.vtk']);
p=model.degree;
Nny=length(yo)+p(1)-1;
Nnx=length(xo)+p(2)-1;
pscale=2^(1-1);
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','ls1','ls2','nx','ny','nz','sizeim','nband');
nim=size(U,2);
non=length(nband);





for iim=1:nim
    if nim==1
        fildef=param.deformed_image;
    else
        fildef=param.deformed_image{kkk};
    end
    load(fildef,'jm3');
    jm3=double(jm3);
    set.vtkname=[filres , sprintf('-images-%05d.vtk',images(iim))];
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
xi=zeros(non,1);
yi=zeros(non,1);
zi=zeros(non,1);
Ux=zeros(non,1);
Uy=zeros(non,1);
Uz=zeros(non,1);
im0on=zeros(non,1);
im1on=zeros(non,1);
npix=0;
        for iy=1:Nelems(1)
            if iy==Nelems(1)
                condy=(ls1(:)>=yo(iy))&(ls1(:)<=yo(iy+1));
            else
                condy=(ls1(:)>=yo(iy))&(ls1(:)<yo(iy+1));
            end
            for ix=1:Nelems(2)
                if ix==Nelems(2)
                    condx=(ls2(:)>=xo(ix))&(ls2(:)<=xo(ix+1));
                else
                    condx=(ls2(:)>=xo(ix))&(ls2(:)<xo(ix+1));
                end

                cond=condx&condy;
                yp=ls1(cond);
                xp=ls2(cond);
                nxp=-nx(cond);
                nyp=-ny(cond);
                nzp=-nz(cond);
                [N1y]=NURBSBasisFunc(iy+p(1),p(1),yp',vo,0);
                [N1x]=NURBSBasisFunc(ix+p(2),p(2),xp',uo,0);
                [xpix,ypix,zpix]=ind2sub(sizeim,nband(cond));
                [indx,indy]=meshgrid(1:size(N1x,2),1:size(N1y,2));
                N=N1x(:,indx(:),1).*N1y(:,indy(:),1);
                indn=sub2ind([Nny,Nnx],indy(:)+iy-1,indx(:)+ix-1);
                Un=N*U(indn,iim);
                im1e=mexInterpLinear3D(xpix(:)+nxp.*Un+roi(1)-1,...
                    ypix(:)+nyp.*Un+roi(3)-1,...
                    zpix(:)+nzp.*Un+roi(5)-1,...
                    jm3);
indp=npix+(1:numel(xpix));
Ux(indp)=nxp.*Un;
Uy(indp)=nyp.*Un;
Uz(indp)=nzp.*Un;
                im0on(indp)=im0(nband(cond));
im1on(indp)=im1e;
xi(indp)=xpix(:);
yi(indp)=ypix(:);
zi(indp)=zpix(:);
npix=npix+numel(xpix);
            end
        end
        im0on=uint8(floor(255*(im0on-min(im0on(:)))/(max(im0on(:))-min(im0on(:)))));
        im1on=uint8(floor(255*(im1on-min(im1on(:)))/(max(im1on(:))-min(im1on(:)))));
    count = fprintf(fwid,'DATASET UNSTRUCTURED_GRID\n');
    count = fprintf(fwid,'POINTS %u float\n',non);
                    data=[(xi-1)+roi(1)-1+zone(1),(yi-1)+roi(3)-1+zone(3),(zi-1)+roi(5)-1+zone(5)]';
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
    count = fprintf(fwid,['VECTORS D float\n']);
    UU=[Ux';Uy';Uz'];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
                count = fprintf(fwid,['SCALARS IM0 int, 1\n']);
 count = fprintf(fwid,'LOOKUP_TABLE default\n');
 if set.ascii
   fprintf(fwid, '%e\n', im0on);
else
 fwrite(fwid,im0on,'uint8');     
 end
                count = fprintf(fwid,['SCALARS IM1 int, 1\n']);
 count = fprintf(fwid,'LOOKUP_TABLE default\n');
 if set.ascii
   fprintf(fwid, '%e\n', im1on);
else
 fwrite(fwid,im1on,'uint8');     
 end

    fclose(fwid);
    xi=xi-min(xi)+1;
    yi=yi-min(yi)+1;
    zi=zi-min(zi)+1;
    sizeim=[max(xi)-min(xi)+1,max(yi)-min(yi)+1,max(zi)-min(zi)+1];
    ind=sub2ind(sizeim,xi,yi,zi);
    im1on(ind)=im1on;
        fwid = fopen([strrep(fildef,'.mat',''),sprintf('-unwarped-%d-%d-%d.raw',sizeim)],'w'); % IMPORTANT: big endian
fwrite(fwid,im1on);
fclose(fwid);

    
    
    
end
end