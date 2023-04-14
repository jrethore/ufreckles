function writeVTKmesh25D(filres,selected2)
nmod=1;

tic
load([filres,'.mat'])
if (nargin<2)
    if ~exist('selected')
        selected=ones(length(Xo),1);
    end
else
        selected=selected2;
end
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
set.vtkname=[filres,'-topo.vtk'];
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

data=[Xo,Yo,Zo]';

if set.ascii
    fprintf(fwid, '%f %f %f \n', data);
else
    fwrite(fwid, data,'double');
end
fprintf(fwid, '\n', data);
count = fprintf(fwid,'CELLS %u %u\n',nt3+nq4+nb2,4*nt3+5*nq4+3*nb2);
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
data=[repmat(2,1,nb2);belt'-1];


if set.ascii
    fprintf(fwid, '%d %d %d\n', data);
else
    fwrite(fwid, data,'uint');
end

count = fprintf(fwid,'CELL_TYPES %u\n',ne+nb2);
data=[repmat(5,1,nt3),repmat(9,1,nq4),repmat(2,1,nb2)];

if set.ascii
    fprintf(fwid, '%d\n', data);
else
    fwrite(fwid, data,'uint');
end
%

fclose(fwid);
fprintf(1,'vtkexport done in %5.3f s\n',toc);
end