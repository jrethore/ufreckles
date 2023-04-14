function writeVTKmesh3D(filreso,selected2)
nmod=1;

tic
[pp,filres,ext]=fileparts(filreso);
if isempty(ext) filreso=[filreso,'.mat'];end
load(filreso,'-mat')

if (nargin<2)
    if ~exist('selected','var')
        selected=ones(length(xo),1);
    end
else
    selected=selected2;
end
    conn=[conn,zeros(size(conn,1),8-size(conn,2))];

%%
set.ascii=1;
set.remark=' computed by MIC';
nn=numel(xo);
ne=length(elt);
 foundt3=find(elt==6);
 foundt4=find(elt==4);
 foundq4=find(elt==8);
 foundt=[foundt3;foundt4;foundq4];
nt4=sum(elt==4);
np6=sum(elt==6);
nh8=sum(elt==8);
belt=[];
for ie=1:length(elt)
    inods=conn(ie,1:elt(ie));
    if sum(~selected(inods))>1
        found=find(~selected(inods));
        for ib=1:(length(found)-1)
            belt=[belt;inods(found(ib+(0:1)))];
        end
    end

end
nb2=size(belt,1);
set.vtkname=[filres,'.vtk'];
%%
fwid = fopen(set.vtkname,'w'); % IMPORTANT: big endian
count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
count = fprintf(fwid,[set.vtkname,set.remark,'\n']);
if set.ascii
    count = fprintf(fwid,'ASCII\n');
else
    count = fprintf(fwid,'BINARY\n');
end
count = fprintf(fwid,'DATASET UNSTRUCTURED_GRID\n');
count = fprintf(fwid,'POINTS %u double\n',nn);

data=[xo,yo,zo]';

if set.ascii
    fprintf(fwid, '%f %f %f \n', data);
else
    fwrite(fwid, data,'double');
end

 count = fprintf(fwid,'CELLS %u %u\n',nt4+np6+nh8,5*nt4+7*np6+9*nh8);
data=[repmat(4,1,length(foundt4));conn(foundt4,1:4)'-1];


if set.ascii
   fprintf(fwid, '%d %d %d %d %d\n', data);
else
      fwrite(fwid, data,'uint');
end
data=[repmat(6,1,length(foundt3));conn(foundt3,1:6)'-1];


if set.ascii
   fprintf(fwid, '%d %d %d %d %d %d %d\n', data);
else
      fwrite(fwid, data,'uint');
end

data=[repmat(8,1,length(foundq4));conn(foundq4,1:8)'-1];


if set.ascii
   fprintf(fwid, '%d %d %d %d %d %d %d %d %d\n', data);
else
      fwrite(fwid, data,'uint');
end


      count = fprintf(fwid,'CELL_TYPES %u\n',ne);
data=[repmat(10,1,nt4),repmat(13,1,np6),repmat(12,1,nh8)];

if set.ascii
   fprintf(fwid, '%d\n', data);
else
fwrite(fwid, data,'uint');
end



%

fclose(fwid);
fprintf(1,'vtkexport done in %5.3f s\n',toc);
end