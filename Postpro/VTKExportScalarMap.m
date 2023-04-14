function VTKExportScalarMap(filres,var,xo,yo,map,dx)

if nargin<6, dx=1;end

% xo=xo-xo(1);
% yo=yo-yo(1);

%%
set.vtkname=[ filres '-'  var  '.vtk'];
set.ascii=0;
set.remark=set.vtkname;
set.dim=[(xo(length(xo))-xo(1)+1) (yo(length(yo))-yo(1)+1) 1];
set.origin=[xo(1) yo(1) 0];
set.spacing=[dx dx 0];
set.varname=var;
np=size(map,1)*size(map,2);
%%
fwid = fopen(fullfile('VTK',set.vtkname),'w','b'); % IMPORTANT: big endian
count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
count = fprintf(fwid,[set.remark,'\n']);
if set.ascii
    count = fprintf(fwid,'ASCII\n');
else
    count = fprintf(fwid,'BINARY\n');
end
count = fprintf(fwid,'DATASET STRUCTURED_POINTS\n');
count = fprintf(fwid,'DIMENSIONS %u %u %u\n',set.dim);
%count = fprintf(fwid,'DIMENSIONS %5.3e %5.3e %5.3e\n',set.dim);
count = fprintf(fwid,'ORIGIN %u %u %u\n',set.origin);
%count = fprintf(fwid,'SPACING %3.2f %3.2f %3.2f\n',set.spacing);
count = fprintf(fwid,'SPACING %5.3e %5.3e %5.3e\n',set.spacing);
count = fprintf(fwid,'POINT_DATA %u\n',np);
count = fprintf(fwid,['SCALARS ',set.varname,' float',' 1\n']);
%count = fprintf(fwid,['SCALARS ',set.varname,' ',vtkprecision,' 1\n']);
count = fprintf(fwid,'LOOKUP_TABLE default\n');

% write data to vtk file
tic
display(sprintf('writing %s...',set.vtkname));
%if ~set.ascii
%    count = fwrite(fwid, M, set.precision);
%else
%    fprintf(fwid, '%g \n', map);
    fwrite(fwid, map,'float');
%end
fclose(fwid);
display(sprintf('vtkexport done in %5.3f s',toc));
