function WritePVDFile(filres,images,time)
fwid = fopen(fullfile('VTK',sprintf('%s.pvd',filres)),'w');
count = fprintf(fwid,'<?xml version="1.0"?>\n');
count = fprintf(fwid,'<VTKFile type="Collection" version="0.1">\n');
count = fprintf(fwid,'  <Collection>\n');
for iim=1:numel(images)
filname=sprintf('%s-%05d.vtk',filres,images(iim));
count = fprintf(fwid,'    <DataSet timestep="%f" group="" part="0"file="%s"/>\n',time(iim),filname);
end
count = fprintf(fwid,'  </Collection>\n');
count = fprintf(fwid,'</VTKFile>\n');
end