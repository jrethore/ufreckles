function postproVIC3D(filres)
display(sprintf('In postproVIC3D\nResult file : %s',filres));
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
nn=prod(Nnodes);
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'ls1','ls2','sizeim','nband','on');
load(fullfile('TMP',sprintf('%d_3d_phix_%d',nmod,(iscale-1))),'phix','Xi','Yi','Zi');
load(fullfile('TMP',sprintf('%d_3d_phiy_%d',nmod,(iscale-1))),'phiy');
load(fullfile('TMP',sprintf('%d_3d_phiz_%d',nmod,(iscale-1))),'phiz');
non=length(on);
unix(['rm ',fullfile('VTK',filres),'-0*.vtk']);
unix(['rm ',filres,'-xyzon-0*.dat']);
if size(U,2)==1
    U=[zeros(size(U,1),1),U];
    images=[0,images];
end
for iim=1:size(U,2)
    set.vtkname=[filres , sprintf('-%06d.vtk',images(iim))];
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
    data=[(Xi-1)+roi(1)-1+zone(1),(Yi-1)+roi(3)-1+zone(3),(Zi-1)+roi(5)-1+zone(5)]';


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
    Uxi=phix*U(:,iim);
    Uyi=phiy*U(:,iim);
    Uzi=phiz*U(:,iim);
    UU=[Uxi';Uyi';Uzi'];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
    
    dlmwrite([filres , sprintf('-xyzon-%06d.dat',images(iim))],[(Xi-1)+roi(1)-1+zone(1)+Uxi,(Yi-1)+roi(3)-1+zone(3)+Uyi,(Zi-1)+roi(5)-1+zone(5)+Uzi]);

    
                    count = fprintf(fwid,['SCALARS LS1 double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%f\n', ls1(on));
                    else
                        fwrite(fwid, ls1(on),'double');
                    end
                    count = fprintf(fwid,['SCALARS LS2 double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%f\n', ls2(on));
                    else
                        fwrite(fwid, ls2(on),'double');
                    end

    fclose(fwid);
end






end