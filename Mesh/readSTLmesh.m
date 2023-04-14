function [xo,yo,zo,conn,elt]=readSTLmesh(filename)
fwid = fopen(filename,'r');
FN = [];
while isempty(FN)
    FN=sscanf(strtrim(fgetl(fwid)),'facet normal %f %f %f');
end
xyzo=[];
invc=[];
ie=1;
while ~isempty(FN)
    count = fscanf(fwid,'  outer loop\n');
        xyz = fscanf(fwid,'    vertex %f %f %f\n');
xyz=reshape(xyz,3,3);
nn=cross(diff(xyz(:,1:2),[],2),diff(xyz(:,2:3),[],2));
    sc=dot(nn,FN);
    if sc>0
        invc(ie)=0;
    else
        invc(ie)=1;
    end
    count = fscanf(fwid,'  endloop\n');
xyzo=[xyzo;xyz'];

count = fscanf(fwid,'endfacet\n');
FN = fscanf(fwid,'facet normal %f %f %f\n');
ie=ie+1;
end
fclose(fwid);

conn=reshape((1:size(xyzo,1))',3,length(invc))';
conn(invc>0,:)=conn(invc>0,[3,2,1]);
[xo,yo,conn,zo]=CleanTriMesh(xyzo(:,1),xyzo(:,2),conn,xyzo(:,3));
elt=3*ones(size(conn,1),1);


end