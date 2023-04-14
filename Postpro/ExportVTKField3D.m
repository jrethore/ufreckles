function ExportVTKField3D(filres,images,xo,yo,zo,elt,conn,smap,sname)

set.ascii=0;
set.remark=' computed by UFreckles';
nn=length(xo);
ne=length(elt);
 foundt3=find(elt==6);
 foundt4=find(elt==4);
 foundq4=find(elt==8);
 foundt=[foundt3;foundt4;foundq4];
nt4=sum(elt==4);
np6=sum(elt==6);
nh8=sum(elt==8);
    conn=[conn,zeros(size(conn,1),8-size(conn,2))];

display(sprintf('Result file : %s',filres));
delete([fullfile('VTK',filres),'-0*.vtk']);
images=[0,images];
tic;
for iim=1:length(images)
    
    set.vtkname=[filres , sprintf('-%06d.vtk',images(iim))];
    if set.ascii
        fwid = fopen(fullfile('VTK',set.vtkname),'w');
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
    count = fprintf(fwid,'POINTS %u float\n',nn);
    data=[xo,yo,zo]';
    if set.ascii
        fprintf(fwid, '%f %f %f \n', data);
    else
        fwrite(fwid, data,'float');
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
    
    count = fprintf(fwid,'POINT_DATA %u\n',nn);
    for kk=1:length(smap)
        FF=smap{kk};
        FF=FF(:,max(1,iim-1));
        if iim==1,FF=0*FF;end
        switch size(FF,1)/nn
            case 1
                count = fprintf(fwid,['SCALARS %s float, 1\n'],sname{kk});
                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                if set.ascii
                    fprintf(fwid, '%f\n', FF);
                else
                    fwrite(fwid,FF,'float');
                end
                
            case 3
                count = fprintf(fwid,['VECTORS %s float\n'],sname{kk});
                UU=[FF(1:nn)';FF(nn+(1:nn))';FF(2*nn+(1:nn))'];
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', UU);
                else
                    fwrite(fwid,UU,'float');
                end
            case 6
     count = fprintf(fwid,['TENSORS %s float\n'],sname{kk});
     E=[FF((1:nn)),FF(3*nn+(1:nn)),FF(5*nn+(1:nn)),FF(3*nn+(1:nn)),FF(1*nn+(1:nn)),FF(4*nn+(1:nn)),FF(5*nn+(1:nn)),FF(4*nn+(1:nn)),FF(2*nn+(1:nn))]';
    if set.ascii
        fprintf(fwid, '%f %f %f \n', E);
    else
        fwrite(fwid, E,'float');
    end
                
                
        end
        
        
    end
    
    count = fprintf(fwid,'CELL_DATA %u\n',ne);
    for kk=1:length(smap)
        FF=smap{kk};
        FF=FF(:,max(1,iim-1));
        if iim==1,FF=0*FF;end
        
        switch size(FF,1)/ne
            case 1
                count = fprintf(fwid,['SCALARS %s float, 1\n'],sname{kk});
                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                if set.ascii
                    fprintf(fwid, '%f\n', FF);
                else
                    fwrite(fwid,FF,'float');
                end
                
            case 3
                count = fprintf(fwid,['VECTORS %s float\n'],sname{kk});
                UU=[FF(1:ne)';FF(ne+(1:ne))';FF(2*ne+(1:ne))'];
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', UU);
                else
                    fwrite(fwid,UU,'float');
                end
            case 6
     count = fprintf(fwid,['TENSORS %s float\n'],sname{kk});
     E=[FF((1:ne)),FF(3*ne+(1:ne)),FF(5*ne+(1:ne)),FF(3*ne+(1:ne)),FF(1*ne+(1:ne)),FF(4*ne+(1:ne)),FF(5*ne+(1:ne)),FF(4*ne+(1:ne)),FF(2*ne+(1:ne))]';
    if set.ascii
        fprintf(fwid, '%f %f %f \n', E);
    else
        fwrite(fwid, E,'float');
    end
              
                
        end
        
        
    end
    fclose(fwid);  
end

fprintf(1,'vtkexport done in %5.3f s\n',toc);
end