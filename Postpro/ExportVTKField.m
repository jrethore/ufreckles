function ExportVTKField(filres,images,xo,yo,elt,conn,smap,sname)
set.ascii=0;
set.remark=' computed by UFreckles';
nn=length(xo);
ne=length(elt);
foundb2=find(elt==2);
foundt3=find(elt==3);
foundq4=find(elt==4);
foundt=[foundt3;foundq4];
nb2=sum(elt==2);
nt3=sum(elt==3);
nq4=sum(elt==4);
    conn=[conn,zeros(size(conn,1),4-size(conn,2))];
    
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
    data=[xo,yo,0*yo]';
    if set.ascii
        fprintf(fwid, '%f %f %f \n', data);
    else
        fwrite(fwid, data,'float');
    end
    count = fprintf(fwid,'CELLS %u %u\n',nt3+nq4+nb2,4*nt3+5*nq4+3*nb2);
    data=[repmat(2,1,length(foundb2));conn(foundb2,1:2)'-1];
    
    if set.ascii
        fprintf(fwid, '%d %d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
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
    data=[repmat(3,1,nb2),repmat(5,1,nt3),repmat(9,1,nq4)];
    
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
                    fwrite(fwid,full(FF),'float');
                end
                
            case 2
                count = fprintf(fwid,['VECTORS %s float\n'],sname{kk});
                UU=[FF(1:nn)';FF(nn+(1:nn))';repmat(0,1,nn)];
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', UU);
                else
                    fwrite(fwid,full(UU),'float');
                end
              case {3,4}
     count = fprintf(fwid,['TENSORS %s float\n'],sname{kk});
     E0=zeros(nn,1);
     if size(FF,1)<4*nn
         FF=[FF;zeros(nn,1)];
     end
     E=[FF((1:nn)),FF(2*nn+(1:nn)),E0,FF(2*nn+(1:nn)),FF(1*nn+(1:nn)),E0,E0,E0,FF(3*nn+(1:nn))]';
    if set.ascii
        fprintf(fwid, '%f %f %f \n', E);
    else
        fwrite(fwid, full(E),'float');
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
                    fwrite(fwid,full(FF),'float');
                end
                
            case 2
                count = fprintf(fwid,['VECTORS %s float\n'],sname{kk});
                UU=[FF(1:ne)';FF(ne+(1:ne))';repmat(0,1,ne)];
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', UU);
                else
                    fwrite(fwid,full(UU),'float');
                end
            case {3,4}
     count = fprintf(fwid,['TENSORS %s float\n'],sname{kk});
     E0=zeros(ne,1);
     if size(FF,1)<4*ne
         FF=[FF;zeros(ne,1)];
     end
     E=[FF((1:ne)),FF(2*ne+(1:ne)),E0,FF(2*ne+(1:ne)),FF(1*ne+(1:ne)),E0,E0,E0,FF(3*ne+(1:ne))]';
    if set.ascii
        fprintf(fwid, '%f %f %f \n', E);
    else
        fwrite(fwid, full(E),'float');
    end
              
                
        end
        
        
    end
    fclose(fwid);  
end

fprintf(1,'vtkexport done in %5.3f s\n',toc);
end