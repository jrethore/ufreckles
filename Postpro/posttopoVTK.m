function posttopoVTK(filres)
nmod=1;

tic
load([filres,'.mat'])

%%
set.ascii=1;
set.remark=' computed by MIC';
nn=prod(Nnodes);
ne=length(elt);
nt3=sum(elt==3);
nq4=sum(elt==4);
 foundt3=find(elt==3);
 foundq4=find(elt==4);
 foundt=[foundt3;foundq4];
C=param.calibration_data{1};
set.vtkname=[filres,'-topo.vtk'];
%%
fwid = fopen(fullfile('VTK',set.vtkname),'w'); % IMPORTANT: big endian
count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
count = fprintf(fwid,[set.vtkname,set.remark,'\n']);
if set.ascii
    count = fprintf(fwid,'ASCII\n');
else
    count = fprintf(fwid,'BINARY\n');
end
count = fprintf(fwid,'DATASET UNSTRUCTURED_GRID\n');
count = fprintf(fwid,'POINTS %u float\n',nn);
data=[Xo,Yo,Zo]';


if set.ascii
   fprintf(fwid, '%f %f %f \n', data);
else  
    fwrite(fwid, data,'float');
end
 count = fprintf(fwid,'CELLS %u %u\n',nt3+nq4,4*nt3+5*nq4);
data=[repmat(3,1,length(foundt3));conn(foundt3,1:3)'-1];


if set.ascii
   fprintf(fwid, '%d %d %d %d\n', data);
else
      fwrite(fwid, data,'uint');
end
data=[repmat(4,1,length(foundq4));conn(foundq4,1:4)'-1];


if set.ascii
   fprintf(fwid, '%d %d %d %d %d\n', data);
else
      fwrite(fwid, data,'uint');
end


      count = fprintf(fwid,'CELL_TYPES %u\n',ne);
data=[repmat(5,1,nt3),repmat(9,1,nq4)];

if set.ascii
   fprintf(fwid, '%d\n', data);
else
fwrite(fwid, data,'uint');
end
% 
 count = fprintf(fwid,'POINT_DATA %u\n',nn);
 count = fprintf(fwid,['SCALARS Elevation float, 1\n']);
count = fprintf(fwid,'LOOKUP_TABLE default\n');
 UU=[Zo'];
if set.ascii
   fprintf(fwid, '%f\n', UU);
else
 fwrite(fwid,UU,'float');     
end

fclose(fwid);
fprintf(1,'vtkexport done in %5.3f s\n',toc);
