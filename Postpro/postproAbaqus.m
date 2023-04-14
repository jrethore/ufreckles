function postproAbaqus(filres,submean)
nmod=1;

tic
load([filres,'.mat'])
if nargin<2,submean=0;end
if isfield(param,'image_number')
    images=param.image_number;
else
    images=1:size(U,2);
end
images=[0,images];
U=[zeros(size(U,1),1),U];
roi=param.roi;
dflag=exist('zo');
rint=true;
if isfield(model,'reduced_integration')
    rint=model.reduced_integration;
end
nlflag=false;
if isfield(model,'nlgeom')
    nlflag=model.nlgeom;
end
cpflag=false;
if isfield(model,'plane_stress')
    cpflag=model.plane_stress;
end
shell=false;
if isfield(model,'shell_element')
    shell=model.shell_element;
end
cpus=1;
if isfield(param,'abaqus_cpus');
    cpus=param.abaqus_cpus;
end
abqcmd='abaqus';
if isfield(param,'abaqus_command_line')
    abqcmd=param.abaqus_command_line;
end
filres=[filres,'-abaqus'];
filres0=[filres,'.inp'];

if exist(filres0,'file')
    %            unix(['rm ',strrep(filres0,'.inp','.*')]);
end



%%
fwid = fopen(filres0,'w');
count = fprintf(fwid,'*HEADING\n');
count = fprintf(fwid,'*Preprint,echo=NO,history=NO,model=NO, contact=NO\n');
count = fprintf(fwid,'UFreckles - CNRS\n');
count = fprintf(fwid,'Laboratoire de Mecanique des Contacts et des Structures - INSA de Lyon\n');
count = fprintf(fwid,'Institut de Recherche Génie Civil et Mecanique - Centrale Nantes\n');


count = fprintf(fwid,'*Node, System=R,Nset=allnodes\n');
if rint
    ng=ones(length(elt),1);
else
    ng=4*(elt==4)+1*(elt==3)+2*(elt==6)+8*(elt==8);
end
ngc=[0;cumsum(ng(1:length(ng)-1))];
if dflag
    foundt3=find(elt==6);
    foundq4=find(elt==8);
    data=[(1:prod(Nnodes))',xo,yo,zo]';
    count = fprintf(fwid, '%d, %f, %f , %f \n', data);
    if any(foundt3)
        if shell
            count = fprintf(fwid,'*Element, type=SC6R,Elset=MAIN\n');
        else
            if rint
                count = fprintf(fwid,'*Element, type=C3D6R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=C3D6,Elset=MAIN\n');
            end
        end
        
        data=[(1:length(foundt3))',conn(foundt3,1:6)]';
        count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d\n',data);
    end
    if any(foundq4)
        if shell
            count = fprintf(fwid,'*Element, type=SC8R,Elset=MAIN\n');
        else
            if rint
                count = fprintf(fwid,'*Element, type=C3D8R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=C3D8,Elset=MAIN\n');
            end
        end
        
        data=[(1:length(foundq4))',conn(foundq4,1:8)]';
        count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',data);
    end
else
    data=[(1:prod(Nnodes))',xo,yo]';
    foundt3=find(elt==3);
    foundq4=find(elt==4);
    count = fprintf(fwid, '%d, %f, %f \n', data);
    if any(foundt3)
        if cpflag
            count = fprintf(fwid,'*Element, type=CPS3,Elset=MAIN\n');
        else
            count = fprintf(fwid,'*Element, type=CPE3,Elset=MAIN\n');
        end
        data=[(1:length(foundt3))',conn(foundt3,1:3)]';
        count = fprintf(fwid,'%d, %d, %d, %d\n',data);
    end
    if any(foundq4)
        if cpflag
            if rint
                count = fprintf(fwid,'*Element, type=CPS4R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=CPS4,Elset=MAIN\n');
            end
        else
            if rint
                count = fprintf(fwid,'*Element, type=CPE4R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=CPE4,Elset=MAIN\n');
            end
        end
        data=[(1:length(foundq4))',conn(foundq4,1:4)]';
        count = fprintf(fwid,'%d, %d, %d, %d, %d\n',data);
    end
end
%%
if isfield(model,'material_model')
    matmod=model.material_model;
    if shell
        count = fprintf(fwid,'*SHELL SECTION,ELSET=MAIN,NODAL THICKNESS,MATERIAL=%s\n',matmod);
    else
        count = fprintf(fwid,'*SOLID SECTION,ELSET=MAIN,MATERIAL=%s\n',matmod);
    end
    count = fprintf(fwid,'*MATERIAL,NAME=%s\n',matmod);
    
    buffer=writeAbaqusMat(matmod,model.material_parameters);
    
    count = fprintf(fwid,'%s',buffer);
else
    if shell
        count = fprintf(fwid,'*SHELL SECTION,ELSET=MAIN,NODAL THICKNESS,MATERIAL=elastic_homogeneous_isotropic\n');
    else
        count = fprintf(fwid,'*SOLID SECTION,ELSET=MAIN,MATERIAL=elastic_homogeneous_isotropic\n');
    end
    count = fprintf(fwid,'*MATERIAL,NAME=elastic_homogeneous_isotropic\n');
    
    count = fprintf(fwid,'*ELASTIC,TYPE=ISOTROPIC\n');
    count = fprintf(fwid,'%12.5e, %12.5e\n',1,0.3);
    
    
end
%%
for iim=1:size(U,2)
    if nlflag
        count = fprintf(fwid,'*STEP,NLGEOM=YES\n');
    else
        count = fprintf(fwid,'*STEP,NLGEOM=NO\n');
    end
    count = fprintf(fwid,'*STATIC\n');
    count = fprintf(fwid,'1,1,0.0000001,1\n');
    count = fprintf(fwid,'*BOUNDARY, OP=NEW\n');
    
    Ux=U((1:prod(Nnodes)),iim);
    Uy=U(prod(Nnodes)+(1:prod(Nnodes)),iim);
    data=[(1:prod(Nnodes))',repmat(1,prod(Nnodes),2),Ux-mean(Ux(:))]';
    count = fprintf(fwid, '%d, %d, %d, %12.5e\n', data);
    data=[(1:prod(Nnodes))',repmat(2,prod(Nnodes),2),Uy-mean(Uy(:))]';
    
    count = fprintf(fwid, '%d, %d, %d, %12.5e\n', data);
    if dflag
        Uz=U(2*prod(Nnodes)+(1:prod(Nnodes)),iim);
        
        data=[(1:prod(Nnodes))',repmat(3,prod(Nnodes),2),Uz-mean(Uz(:))]';
        count = fprintf(fwid, '%d, %d, %d, %12.5e\n', data);
        
    end
    for io=1:length(model.abaqus_output)
        count = fprintf(fwid,sprintf('*EL PRINT, FREQUENCY=1000000000,ELSET=MAIN\n%s\n',model.abaqus_output{io}));
    end
    
    count = fprintf(fwid,'*END STEP\n');
    
end


%%
fclose(fwid);
if cpus==1
    unix([abqcmd,' job=',filres,' double  interactive']);
else
    unix([abqcmd,' job=',filres,' cpus=',num2str(cpus),' double  interactive']);
end
fprintf(1,'abaqusexport done in %5.3f s\n',toc);
%%
frid = fopen([filres,'.dat'],'r');

set.ascii=1;
set.remark=' computed by UFreckles';
nn=prod(Nnodes);
ne=length(elt);
np6=sum(elt==6);
nh8=sum(elt==8);
nt3=sum(elt==3);
nq4=sum(elt==4);
if isfield(param,'calibration_data')
    roi=0*roi+1;
end
unix(['rm ',fullfile('VTK',filres),'-0*.vtk']);
for iim=1:size(U,2)
    set.vtkname=[ filres , sprintf('-%05d.vtk',images(iim))];
    %fwid = fopen(set.vtkname,'w','b'); % IMPORTANT: big endian
    fwid = fopen(fullfile('VTK',set.vtkname),'w'); % IMPORTANT: big endian
    count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
    count = fprintf(fwid,[set.vtkname,set.remark,'\n']);
    if set.ascii
        count = fprintf(fwid,'ASCII\n');
    else
        count = fprintf(fwid,'BINARY\n');
    end
    count = fprintf(fwid,'DATASET UNSTRUCTURED_GRID\n');
    count = fprintf(fwid,'POINTS %u double\n',nn);
    if dflag
        data=[xo+roi(1)-1,yo+roi(3)-1,zo]';
    else
        data=[xo+roi(1)-1,yo+roi(3)-1,0*yo]';
    end
    if set.ascii
        fprintf(fwid, '%f %f %f \n', data);
    else
        fwrite(fwid, data,'double');
    end
    if dflag
        count = fprintf(fwid,'CELLS %u %u\n',np6+nh8,7*np6+9*nh8);
        data=[repmat(6,1,length(foundt3));conn(foundt3,1:6)'-1];
        
        if set.ascii
            fprintf(fwid, '%d %d %d %d\n', data);
        else
            fwrite(fwid, data,'uint');
        end
        
        data=[repmat(8,1,length(foundq4));conn(foundq4,1:8)'-1];
        
        
        if set.ascii
            fprintf(fwid, '%d %d %d %d %d\n', data);
        else
            fwrite(fwid, data,'uint');
        end
        
        
        count = fprintf(fwid,'CELL_TYPES %u\n',ne);
        data=[repmat(13,1,np6),repmat(12,1,nh8)];
        
        if set.ascii
            fprintf(fwid, '%d\n', data);
        else
            fwrite(fwid, data,'uint');
        end
    else
        count = fprintf(fwid,'CELLS %u %u\n',nt3+nq4,4*nt3+5*nq4);
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
        data=[repmat(5,1,nt3),repmat(9,1,nq4)];
        
        if set.ascii
            fprintf(fwid, '%d\n', data);
        else
            fwrite(fwid, data,'uint');
        end
    end
    Ui=U(:,iim);
    if submean
        Ui(1:nn)=Ui(1:nn)-mean(Ui(1:nn));
        Ui(nn+(1:nn))=Ui(nn+(1:nn))-mean(Ui(nn+(1:nn)));
        if dflag
            Ui(2*nn+(1:nn))=Ui(2*nn+(1:nn))-mean(Ui(2*nn+(1:nn)));
        end
    end
    if ~dflag
        UU=[Ui(1:nn)';Ui(nn+(1:nn))';zeros(1,nn)];
    else
        UU=[Ui(1:nn)';Ui(nn+(1:nn))';Ui(2*nn+(1:nn))'];
    end
    
    count = fprintf(fwid,'POINT_DATA %u\n',nn);
    count = fprintf(fwid,['VECTORS Displacement double\n']);
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'double');
    end
    count = fprintf(fwid,'CELL_DATA %u\n',ne);
    for io=1:length(model.abaqus_output)
        switch  model.abaqus_output{io}
            case 'E'
                count = fprintf(fwid,['TENSORS Strain double\n']);
                titi=[];
                while isempty(findstr(titi,'MAIN'))%ELEMENT  PT FOOT'))
                    toto=textscan(frid,'%4c%*[^\n]',1);
                    titi=cell2mat(toto);
                end
                toto=textscan(frid,'%1c%*[^\n]',1);
                toto=textscan(frid,'%1c%*[^\n]',1);
                
                if dflag
                    data=fscanf(frid,'%d%d%f%f%f%f%f%f');
                    data=reshape(data',8,length(data)/8)';
                    indj=data(:,1);
                    E=data(:,2+[1,4,5,4,2,6,5,6,3]);
                else
                    if cpflag
                        data=fscanf(frid,'%d%d%f%f%f');
                        data=reshape(data',5,length(data)/5)';
                        indj=data(:,1);
                        Eo=zeros(size(data,1),1);
                        E=[data(:,3),data(:,5),Eo,data(:,5),data(:,4),Eo,Eo,Eo,Eo];
                    else
                        data=fscanf(frid,'%d%d%f%f%f%f');
                        data=reshape(data',6,length(data)/6)';
                        indj=data(:,1);
                        Eo=zeros(size(data,1),1);
                        E=[data(:,3),data(:,6),Eo,data(:,6),data(:,4),Eo,Eo,Eo,data(:,5)];
                        
                    end
                end
                indi=1:size(E,1);
                P=sparse(indj,indi,1./ng(indj),length(elt),size(E,1));
                
                E=(P*E)';
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', E);
                else
                    fwrite(fwid, E,'double');
                end
                 save(sprintf('%s-%s-%05d',filres,model.abaqus_output{io},images(iim)),'E')
           case 'LE'
                count = fprintf(fwid,['TENSORS LE double\n']);
                titi=[];
                while isempty(findstr(titi,'MAIN'))%ELEMENT  PT FOOT'))
                    toto=textscan(frid,'%4c%*[^\n]',1);
                    titi=cell2mat(toto);
                end
                toto=textscan(frid,'%1c%*[^\n]',1);
                toto=textscan(frid,'%1c%*[^\n]',1);
                
                if dflag
                    data=fscanf(frid,'%d%d%f%f%f%f%f%f');
                    data=reshape(data',8,length(data)/8)';
                    indj=data(:,1);
                    E=data(:,2+[1,4,5,4,2,6,5,6,3]);
                else
                    if cpflag
                        data=fscanf(frid,'%d%d%f%f%f');
                        data=reshape(data',5,length(data)/5)';
                        indj=data(:,1);
                        Eo=zeros(size(data,1),1);
                        E=[data(:,3),data(:,5),Eo,data(:,5),data(:,4),Eo,Eo,Eo,Eo];
                    else
                        data=fscanf(frid,'%d%d%f%f%f%f');
                        data=reshape(data',6,length(data)/6)';
                        indj=data(:,1);
                        Eo=zeros(size(data,1),1);
                        E=[data(:,3),data(:,6),Eo,data(:,6),data(:,4),Eo,Eo,Eo,data(:,5)];
                        
                    end
                end
                indi=1:size(E,1);
                P=sparse(indj,indi,1./ng(indj),length(elt),size(E,1));
                
                E=(P*E)';
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', E);
                else
                    fwrite(fwid, E,'double');
                end
                save(sprintf('%s-%s-%05d',filres,model.abaqus_output{io},images(iim)),'E')
            case 'S'
                count = fprintf(fwid,['TENSORS Stress double\n']);
                titi=[];
                while isempty(findstr(titi,'MAIN'))%ELEMENT  PT FOOT'))
                    toto=textscan(frid,'%4c%*[^\n]',1);
                    titi=cell2mat(toto);
                end
                toto=textscan(frid,'%1c%*[^\n]',1);
                toto=textscan(frid,'%1c%*[^\n]',1);
                
                if dflag
                    data=fscanf(frid,'%d%d%f%f%f%f%f%f');
                    data=reshape(data',8,length(data)/8)';
                    indj=data(:,1);
                    E=data(:,2+[1,4,5,4,2,6,5,6,3]);
                else
                    if cpflag
                        data=fscanf(frid,'%d%d%f%f%f');
                        data=reshape(data',5,length(data)/5)';
                        indj=data(:,1);
                        Eo=zeros(size(data,1),1);
                        E=[data(:,3),data(:,5),Eo,data(:,5),data(:,4),Eo,Eo,Eo,Eo];
                    else
                        data=fscanf(frid,'%d%d%f%f%f%f');
                        data=reshape(data',6,length(data)/6)';
                        indj=data(:,1);
                        Eo=zeros(size(data,1),1);
                        E=[data(:,3),data(:,6),Eo,data(:,6),data(:,4),Eo,Eo,Eo,data(:,5)];
                    end
                end
                vm=sqrt(E(:,1).^2+E(:,5).^2+E(:,9).^2-E(:,1).*E(:,9)-E(:,5).*E(:,9)-E(:,1).*E(:,5)...
                    +3*(E(:,2).^2+E(:,3).^2+E(:,6).^2));
                indi=1:size(E,1);
                P=sparse(indj,indi,1./ng(indj),length(elt),size(E,1));
                
                E=(P*E)';
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', E);
                else
                    fwrite(fwid, E,'double');
                end
                
                
                count = fprintf(fwid,['SCALARS VM double, 1\n']);
                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                vm=(P*vm)';
                if set.ascii
                    fprintf(fwid, '%f\n', vm);
                else
                    fwrite(fwid, vm,'double');
                end
                
                
                
                save(sprintf('%s-%s-%05d',filres,model.abaqus_output{io},images(iim)),'E','vm')
                
            case 'PEEQ'
                count = fprintf(fwid,['SCALARS PEEQ double, 1\n']);
                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                titi=[];
                while isempty(findstr(titi,'MAIN'))%ELEMENT  PT FOOT'))
                    toto=textscan(frid,'%4c%*[^\n]',1);
                    titi=cell2mat(toto);
                end
                toto=textscan(frid,'%1c%*[^\n]',1);
                toto=textscan(frid,'%1c%*[^\n]',1);
                
                data=fscanf(frid,'%d%d%f%*s');
                if ~isempty(data)
                    data=reshape(data',3,length(data)/3)';
                    E=zeros(sum(ng),1);
                    indg=ngc(data(:,1))+data(:,2);
                    E(indg)=data(:,3);
                    E=(P*E)';
                else
                    E=repmat(0,size(P,1),1)';
                end
                if set.ascii
                    fprintf(fwid, '%f\n', E);
                else
                    fwrite(fwid, E,'double');
                end
                 save(sprintf('%s-%s-%05d',filres,model.abaqus_output{io},images(iim)),'E')
           case {'MISES','TRESC','PRESS'}
                count = fprintf(fwid,sprintf('SCALARS %s double, 1\n',model.abaqus_output{io}));
                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                titi=[];
                while isempty(findstr(titi,'MAIN'))%ELEMENT  PT FOOT'))
                    toto=textscan(frid,'%4c%*[^\n]',1);
                    titi=cell2mat(toto);
                end
                toto=textscan(frid,'%1c%*[^\n]',1);
                toto=textscan(frid,'%1c%*[^\n]',1);
                
                data=fscanf(frid,'%d%d%f');
                if ~isempty(data)
                    data=reshape(data',3,length(data)/3)';
                    E=zeros(sum(ng),1);
                    indg=ngc(data(:,1))+data(:,2);
                    E(indg)=data(:,3);
                    E=(P*E)';
                else
                    E=repmat(0,size(P,1),1)';
                end
                if set.ascii
                    fprintf(fwid, '%f\n', E);
                else
                    fwrite(fwid, E,'double');
                end
                
                save(sprintf('%s-%s-%05d',filres,model.abaqus_output{io},images(iim)),'E')
            otherwise
                display(sprintf('OUTPUT FIELDS %s NOT AVAILABLE',model.abaqus_output{io}));
        end
    end
    
    fclose(fwid);
end
fclose(frid);
%%

end
