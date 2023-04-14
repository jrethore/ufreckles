function VTKExportScalarVol(filres,var,xo,yo,zo,map,dx)

if nargin<7, dx=1;end

% xo=xo-xo(1);
% yo=yo-yo(1);
% zo=zo-zo(1);

%%
set.vtkname=[filres   '.vtk'];
set.ascii=0;
set.remark=set.vtkname;
%set.dim=[(xo(length(xo))-xo(1)+1) (yo(length(yo))-yo(1)+1) (zo(length(zo))-zo(1)+1)];
set.dim=size(map);
set.origin=[xo(1) yo(1) zo(1)];
if length(dx)==1
set.spacing=[dx dx dx];
else
 set.spacing=[dx(1) dx(2) dx(3)];   
end
set.varname=var;
np=numel(map);
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
fprintf(1,['writing ',set.vtkname,' ... ']);
%if ~set.ascii
%    count = fwrite(fwid, M, set.precision);
%else
%    fprintf(fwid, '%g \n', map);
    fwrite(fwid, map(:),'float');
%end
fclose(fwid);
fprintf(1,'vtkexport done in %5.3f s\n',toc);
