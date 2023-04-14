function writeSTLmesh(filreso)

tic
[pp,filres,ext]=fileparts(filreso);
if isempty(ext) filreso=[filreso,'.mat'];end
load(filreso,'-mat')
if ~exist('zo','var')
    zo=zeros(size(xo));
end
TR=triangulation(conn(:,1:3),[xo,yo,zo]);
FN=faceNormal(TR);
fwid = fopen([filres,'.stl'],'w');
count = fprintf(fwid,'solid Created by UFRECKLES\n');

for ie=1:numel(elt)
inods=conn(ie,1:3);
xn=xo(inods);
yn=yo(inods);
zn=zo(inods);

count = fprintf(fwid,'facet normal %f %f %f\n',FN(ie,:));
count = fprintf(fwid,'  outer loop\n');
for in=1:3
count = fprintf(fwid,'    vertex %f %f %f\n',xn(in),yn(in),zn(in));
    
end


count = fprintf(fwid,'  endloop\n');


count = fprintf(fwid,'endfacet\n');

end
count = fprintf(fwid,'endsolid Created by UFRECKLES\n');
fclose(fwid);





end