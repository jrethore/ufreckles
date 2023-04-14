function genPVSMFile(filres)
load(filres,'param');
nim=1;
if iscell(param.deformed_image)
    nim=size(param.deformed_image,2);
end
%copyfile('pvsm-part-1.txt',[filres,'.pvsm']);
if isfield(param,'image_number')
    images=param.image_number;
else
    images=1:nim;
end
images=[0,images];
nim=nim+1;
fwid = fopen([filres,'.pvsm'],'w');
count = fprintf(fwid,'<ParaView>\n');
count = fprintf(fwid,'  <ServerManagerState version="3.14.1">\n');
count = fprintf(fwid,'    <Proxy group="sources" type="LegacyVTKFileReader" id="3487" servers="1">\n');
count = fprintf(fwid,'      <Property name="FileNameInfo" id="3487.FileNameInfo" number_of_elements="1">\n');
vtkname=fullfile('VTK',[filres , sprintf('-%05d.vtk',0)]);
count = fprintf(fwid,'        <Element index="0" value="%s"/>\n',vtkname);
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'      <Property name="FileNames" id="3487.FileNames" number_of_elements="%d">\n',nim);
for iim=1:nim
vtkname=fullfile('VTK',[filres , sprintf('-%05d.vtk',images(iim))]);
count = fprintf(fwid,'        <Element index="%d" value="%s"/>\n',iim-1,vtkname);
end
count = fprintf(fwid,'        <Domain name="files" id="3487.FileNames.files"/>\n');
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'      <Property name="TimestepValues" id="3487.TimestepValues" number_of_elements="%d">\n',nim);
for iim=1:nim
count = fprintf(fwid,'        <Element index="%d" value="%d"/>\n',iim-1,iim-1);
end
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'    </Proxy>\n');

count = fprintf(fwid,'    <Proxy group="sources" type="LegacyVTKFileReader" id="3786" servers="1">\n');
count = fprintf(fwid,'      <Property name="FileNameInfo" id="3786.FileNameInfo" number_of_elements="1">\n');
vtkname=fullfile('VTK',[filres , sprintf('-im-%04d.vtk',0)]);
count = fprintf(fwid,'        <Element index="0" value="%s"/>\n',vtkname);
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'      <Property name="FileNames" id="3786.FileNames" number_of_elements="%d">\n',nim);
for iim=1:nim
vtkname=fullfile('VTK',[filres , sprintf('-im-%04d.vtk',images(iim))]);
count = fprintf(fwid,'        <Element index="%d" value="%s"/>\n',iim-1,vtkname);
end
count = fprintf(fwid,'        <Domain name="files" id="3486.FileNames.files"/>\n');
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'      <Property name="TimestepValues" id="3486.TimestepValues" number_of_elements="%d">\n',nim);
for iim=1:nim
count = fprintf(fwid,'        <Element index="%d" value="%d"/>\n',iim-1,iim-1);
end
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'    </Proxy>\n');

count = fprintf(fwid,'    <Proxy group="misc" type="TimeKeeper" id="256" servers="16">\n');
count = fprintf(fwid,'      <Property name="Time" id="256.Time" number_of_elements="1">\n');
count = fprintf(fwid,'        <Element index="0" value="0"/>\n');
count = fprintf(fwid,'        <Domain name="range" id="256.Time.range"/>\n');
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'      <Property name="TimeRange" id="256.TimeRange" number_of_elements="2">\n');
count = fprintf(fwid,'        <Element index="0" value="0"/>\n');
count = fprintf(fwid,'        <Element index="1" value="%d"/>\n',nim-1);
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'      <Property name="TimeSources" id="256.TimeSources" number_of_elements="3">\n');
count = fprintf(fwid,'        <Proxy value="3487"/>\n');
count = fprintf(fwid,'        <Proxy value="3786"/>\n');
count = fprintf(fwid,'        <Proxy value="4084"/>\n');
count = fprintf(fwid,'      </Property>\n');
count = fprintf(fwid,'      <Property name="TimestepValues" id="256.TimestepValues" number_of_elements="%d">\n',nim);
for iim=1:nim
count = fprintf(fwid,'        <Element index="%d" value="%d"/>\n',iim-1,iim-1);
end
count = fprintf(fwid,'      </Property>\n');


%count = fprintf(fwid,'  </ServerManagerState>\n');
%count = fprintf(fwid,'</ParaView>\n');

fclose(fwid);
end