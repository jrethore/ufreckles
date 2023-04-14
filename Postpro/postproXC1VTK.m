function postproXC1VTK(filreso)
if nargin<2,submean=0;end
if nargin<3,do_error=1;end
nmod=1;
[pp,filres,ext]=fileparts(filreso);
if isempty(ext)
    filreso=[filreso,'.mat'];
    ext='.mat';
end
tic
load(filreso,'-mat')
if ~exist('U','var')
    U=U1;
end
if ~exist('model1','var')
    model1=model;
end
if isfield(param,'image_number')
    images=param.image_number;
else
    images=1:size(U,2);
end
if isfield(param,'pixel_size')
    pix2m=param.pixel_size;
else
    pix2m=1;
end
if strcmp(param.analysis,'mechanics')
    dfac=1;
    fac=1+0*pix2m;
else
    fac=1;dfac=1;
end
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
roi=param.roi;
sizeim=[roi(2)-roi(1),roi(4)-roi(3)]+1;
%%
set.ascii=0;
set.remark=' computed by UFRECKLES';

clear wdetJ

save(filreso,'rint','-append');
Us=U;
conns=conn;
[xoo,yoo,~,conno,~,~]=ReadVTK(model.mesh_file);
elto=size(conno,1);
xos=xo;
yos=yo;
nn=numel(xoo);
display(sprintf('Result file : %s',filres));
delete([fullfile('VTK',filres),'-0*.vtk']);
delete([fullfile('VTK',filres),'-x-0*-0*.vtk']);


meshfile='postproXC1VTK';
ns=5*[1,1];

%%iim=0;

set.vtkname=[filres , sprintf('-%06d.vtk',0)];
%%
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
data=[xoo,yoo,0*yoo]';


if set.ascii
    fprintf(fwid, '%f %f %f \n', dfac*data);
else
    fwrite(fwid, data,'float');
end
count = fprintf(fwid,'CELLS %u %u\n',elto,4*elto);
data=[repmat(3,1,size(conno,1));conno(:,1:3)'-1];


if set.ascii
    fprintf(fwid, '%d %d %d %d\n', data);
else
    fwrite(fwid, data,'uint');
end



count = fprintf(fwid,'CELL_TYPES %u\n',elto);
data=[repmat(5,1,elto)];

if set.ascii
    fprintf(fwid, '%d\n', data);
else
    fwrite(fwid, data,'uint');
end
fclose(fwid);

for iz=1:size(model.zone,2)
    zone=model.zone(:,iz);
    switch zone{4}
        case 5 %CRACK
            
            set.vtkname=[filres , sprintf('-x-%03d-%06d.vtk',iz,0)];
            %%
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
            count = fprintf(fwid,'POINTS %u float\n',0);
            count = fprintf(fwid,'CELLS %u %u\n',0,0);
            count = fprintf(fwid,'CELL_TYPES %u\n',0);
            fclose(fwid);
    end
end
%
for iim=1:numel(Us)
    set.vtkname=[filres , sprintf('-%06d.vtk',images(iim))];
    %%
    U=Us{iim};
    conn=conns{iim};
    conn=conn(:,1:3);
    xo=xos{iim};
    yo=yos{iim};
    Nnodes=[numel(xo),1,1];
    lf=max(abs(U(6*prod(Nnodes)+1)));
    U=U/lf;
    Ux=U(1+6*((1:nn)-1),:);
    Uxx=U(2+6*((1:nn)-1),:);
    Uxy=U(3+6*((1:nn)-1),:);
    Uxxx=U(4+6*((1:nn)-1),:);
    Uxxy=U(5+6*((1:nn)-1),:);
    Uxyy=U(6+6*((1:nn)-1),:);
    Uy=U(6*prod(Nnodes)+1+6*((1:nn)-1),:);
    Uyx=U(6*prod(Nnodes)+2+6*((1:nn)-1),:);
    Uyy=U(6*prod(Nnodes)+3+6*((1:nn)-1),:);
    Uyxx=U(6*prod(Nnodes)+4+6*((1:nn)-1),:);
    Uyxy=U(6*prod(Nnodes)+5+6*((1:nn)-1),:);
    Uyyy=U(6*prod(Nnodes)+6+6*((1:nn)-1),:);
    
    face_elt=repmat(((elto+1):size(conn,1))',1,2);
    xg=mean(xo(conn),2);
    yg=mean(yo(conn),2);
    
    for ie=1:size(face_elt,1)
        [~,ide]=min(abs((xg-xg(face_elt(ie,1)))+1i*(yg-yg(face_elt(ie,1)))));
        face_elt(ie,1)=(ide);
    end
    
    conn=conn(1:elto,:);
    conn(face_elt(:,1),:)=[];
    xo=xo(1:nn);
    yo=yo(1:nn);
    
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
    data=[(xo-1)*psample+roi(1),(yo-1)*psample+roi(3),0*yo]';
    
    
    if set.ascii
        fprintf(fwid, '%f %f %f \n', dfac*data);
    else
        fwrite(fwid, data,'float');
    end
    count = fprintf(fwid,'CELLS %u %u\n',size(conn,1),4*size(conn,1));
    data=[repmat(3,1,size(conn,1));conn(:,1:3)'-1];
    
    
    if set.ascii
        fprintf(fwid, '%d %d %d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    
    
    
    count = fprintf(fwid,'CELL_TYPES %u\n',size(conn,1));
    data=[repmat(5,1,size(conn,1))];
    
    if set.ascii
        fprintf(fwid, '%d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    %
    count = fprintf(fwid,'POINT_DATA %u\n',nn);
    count = fprintf(fwid,['VECTORS Displacement float\n']);
    UU=fac*[Ux';Uy';zeros(size(Ux'))];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
    
    
    count = fprintf(fwid,['TENSORS Strain float\n']);
    E0=zeros(nn,1);
    E=[Uxx,0.5*(Uxy+Uyx),E0,0.5*(Uxy+Uyx),Uyy,E0,E0,E0,E0]';
    if set.ascii
        fprintf(fwid, '%f %f %f\n', E);
    else
        fwrite(fwid,E,'float');
    end
    
    count = fprintf(fwid,['TENSORS DStrain float\n']);
    E0=zeros(nn,1);
    E=[Uxxx,Uxxy,0.5*(Uxxy+Uyxx),0.5*(Uxyy+Uyxy),Uyxy,Uyyy,E0,E0,E0]';
    if set.ascii
        fprintf(fwid, '%f %f %f\n', E);
    else
        fwrite(fwid,E,'float');
    end
    
    if isfield(model1,'vtk_export')
        for iex=1:length(model1.vtk_export)
            go=0;
            switch model1.vtk_export{iex}
                case 'S'
                    go=1;
                    H=model1.material_parameters.H;
                    
                    Exx=H(1,1)*Uxx+H(1,2)*Uyy+0.5*H(1,3)*(Uxy+Uyx)*sqrt(2);
                    Eyy=H(2,1)*Uxx+H(2,2)*Uyy+0.5*H(2,3)*(Uxy+Uyx)*sqrt(2);
                    Exy=(H(3,1)*Uxx+H(3,2)*Uyy+0.5*H(3,3)*(Uxy+Uyx)*sqrt(2))/sqrt(2);
                    Eyx=Exy;
                    Exz=sparse(nn,1);
                    Eyz=sparse(nn,1);
                    Ezx=sparse(nn,1);
                    Ezy=sparse(nn,1);
                    Ezz=sparse(nn,1);
                    vm=sqrt(Exx.^2+Eyy.^2+Ezz.^2 ...
                        -Exx.*Ezz-Eyy.*Ezz-Exx.*Eyy ...
                        +3*(Exy.^2+Eyz.^2+Exz.^2));
                    
                    
                    
                    count = fprintf(fwid,['SCALARS VM double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%f\n', full(vm'));
                    else
                        fwrite(fwid, full(vm),'double');
                    end
                    count = fprintf(fwid,['TENSORS Stress double\n']);
                case 'HS'
                    go=1;
                    A=model1.material_parameters.A;
                    Exx=A(1,1)*Uxxx+A(1,2)*Uxxy+(0.5*A(1,3))*(Uxxy+Uyxx)+(0.5*A(1,4))*(Uxyy+Uyxy)+A(1,5)*Uyxy+A(1,6)*Uyyy;
                    Exy=A(2,1)*Uxxx+A(2,2)*Uxxy+(0.5*A(2,3))*(Uxxy+Uyxx)+(0.5*A(2,4))*(Uxyy+Uyxy)+A(2,5)*Uyxy+A(2,6)*Uyyy;
                    Exz=A(3,1)*Uxxx+A(3,2)*Uxxy+(0.5*A(3,3))*(Uxxy+Uyxx)+(0.5*A(3,4))*(Uxyy+Uyxy)+A(3,5)*Uyxy+A(3,6)*Uyyy;
                    Eyx=A(4,1)*Uxxx+A(4,2)*Uxxy+(0.5*A(4,3))*(Uxxy+Uyxx)+(0.5*A(4,4))*(Uxyy+Uyxy)+A(4,5)*Uyxy+A(4,6)*Uyyy;
                    Eyy=A(5,1)*Uxxx+A(5,2)*Uxxy+(0.5*A(5,3))*(Uxxy+Uyxx)+(0.5*A(5,4))*(Uxyy+Uyxy)+A(5,5)*Uyxy+A(5,6)*Uyyy;
                    Eyz=A(6,1)*Uxxx+A(6,2)*Uxxy+(0.5*A(6,3))*(Uxxy+Uyxx)+(0.5*A(6,4))*(Uxyy+Uyxy)+A(6,5)*Uyxy+A(6,6)*Uyyy;
                    Ezx=sparse(nn,1);
                    Ezy=sparse(nn,1);
                    Ezz=sparse(nn,1);
                    
                    
                    count = fprintf(fwid,['TENSORS HStress double\n']);
            end
            if go
                E=[ Exx,Exy,Exz,...
                    Eyx,Eyy,Eyz,...
                    Ezx,Ezy,Ezz]';
                if set.ascii
                    fprintf(fwid, '%e %e %e \n', full(E));
                else
                    fwrite(fwid, full(E),'double');
                end
            end
            
        end
    end
    fclose(fwid);
    Smesh=[1,1,1];
    rflag='false';
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        switch zone{4}
            case 5 %CRACK
                xyc=zone{2};
                
                conni=conns{iim};
                conni=conni(:,1:3);
                xo=xos{iim};
                yo=yos{iim};
                Uxi=U(1:6*numel(xo),:);
                Uyi=U(6*numel(xo)+(1:6*numel(xo)),:);
                Nnodes=[numel(xo),1,1];
                conn=conni(face_elt(:,1),:);
                elt=3*ones(size(conn,1),1);
                ng=0;
                Nelems=[numel(elt),1,1];
                save(meshfile,'xo','yo','Nnodes','Nelems','conn','elt','rint','ng','ns','rflag','Smesh')
                phi=CreateFiniteElementBasis(meshfile,sizeim,1,[],'sub_cells',0,1);
                xg=phi*xo+param.roi(1)-1;
                yg=phi*yo+param.roi(3)-1;
                [maskg,~]=GetSignedDistanceToCrack(xyc,xg+1i*yg);
                
                conn=conni(face_elt(:,1),:);
                save(meshfile,'conn','-append')
                [phi,dphidx,dphidy,dphidxx,dphidyy,dphidxy]=CreateGradC1FiniteElementBasis(meshfile,sizeim,1,[],'sub_cells',[],1);%calcul le gradient des fonctions de formes
                Ux=(maskg>0).*(phi*Uxi);
                Uxx=(maskg>0).*(dphidx*Uxi);
                Uxy=(maskg>0).*(dphidy*Uxi);
                Uxxx=(maskg>0).*(dphidxx*Uxi);
                Uxxy=(maskg>0).*(dphidxy*Uxi);
                Uxyy=(maskg>0).*(dphidyy*Uxi);
                Uy=(maskg>0).*(phi*Uyi);
                Uyx=(maskg>0).*(dphidx*Uyi);
                Uyy=(maskg>0).*(dphidy*Uyi);
                Uyxx=(maskg>0).*(dphidxx*Uyi);
                Uyxy=(maskg>0).*(dphidxy*Uyi);
                Uyyy=(maskg>0).*(dphidyy*Uyi);
                
                
                conn=conni(face_elt(:,2),:);
                save(meshfile,'conn','-append')
                [phi,dphidx,dphidy,dphidxx,dphidyy,dphidxy]=CreateGradC1FiniteElementBasis(meshfile,sizeim,1,[],'sub_cells',[],1);%calcul le gradient des fonctions de formes
                Ux=Ux+(maskg<0).*(phi*Uxi);
                Uxx=Uxx+(maskg<0).*(dphidx*Uxi);
                Uxy=Uxy+(maskg<0).*(dphidy*Uxi);
                Uxxx=Uxxx+(maskg<0).*(dphidxx*Uxi);
                Uxxy=Uxxy+(maskg<0).*(dphidxy*Uxi);
                Uxyy=Uxyy+(maskg<0).*(dphidyy*Uxi);
                Uy=Uy+(maskg<0).*(phi*Uyi);
                Uyx=Uyx+(maskg<0).*(dphidx*Uyi);
                Uyy=Uyy+(maskg<0).*(dphidy*Uyi);
                Uyxx=Uyxx+(maskg<0).*(dphidxx*Uyi);
                Uyxy=Uyxy+(maskg<0).*(dphidxy*Uyi);
                Uyyy=Uyyy+(maskg<0).*(dphidyy*Uyi);
                
                
                
                set.vtkname=[filres , sprintf('-x-%03d-%06d.vtk',iz,iim)];
                %%
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
                
                
                count = fprintf(fwid,'POINTS %u float\n',numel(xg));
                data=[xg,yg,0*yg]';
                
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', dfac*data);
                else
                    fwrite(fwid, data,'float');
                end
                count = fprintf(fwid,'CELLS %u %u\n',numel(xg),3*numel(xg));
                data=[repmat(2,1,numel(xg));repmat((1:numel(xg))-1,2,1)];
                
                if set.ascii
                    fprintf(fwid, '%d %d %d\n', data);
                else
                    fwrite(fwid, data,'uint');
                end
                count = fprintf(fwid,'CELL_TYPES %u\n',numel(xg));
                data=[repmat(2,1,numel(xg))];
                if set.ascii
                    fprintf(fwid, '%d\n', data);
                else
                    fwrite(fwid, data,'uint');
                end
                
                count = fprintf(fwid,'POINT_DATA %u\n',numel(xg));
                count = fprintf(fwid,['VECTORS Displacement float\n']);
                UU=fac*[Ux';Uy';zeros(size(Ux'))];
                if set.ascii
                    fprintf(fwid, '%f %f %f \n', UU);
                else
                    fwrite(fwid,UU,'float');
                end
                count = fprintf(fwid,['TENSORS Strain float\n']);
                E0=zeros(numel(xg),1);
                E=[Uxx,0.5*(Uxy+Uyx),E0,0.5*(Uxy+Uyx),Uyy,E0,E0,E0,E0]';
                if set.ascii
                    fprintf(fwid, '%f %f %f\n', E);
                else
                    fwrite(fwid,E,'float');
                end
                
                count = fprintf(fwid,['TENSORS DStrain float\n']);
                E0=zeros(numel(xg),1);
                E=[Uxxx,Uxxy,0.5*(Uxxy+Uyxx),0.5*(Uxyy+Uyxy),Uyxy,Uyyy,E0,E0,E0]';
                if set.ascii
                    fprintf(fwid, '%f %f %f\n', E);
                else
                    fwrite(fwid,E,'float');
                end
                
                if isfield(model1,'vtk_export')
                    for iex=1:length(model1.vtk_export)
                        go=0;
                        switch model1.vtk_export{iex}
                            case 'S'
                                go=1;
                                H=model1.material_parameters.H;
                                
                                Exx=H(1,1)*Uxx+H(1,2)*Uyy+0.5*H(1,3)*(Uxy+Uyx)*sqrt(2);
                                Eyy=H(2,1)*Uxx+H(2,2)*Uyy+0.5*H(2,3)*(Uxy+Uyx)*sqrt(2);
                                Exy=(H(3,1)*Uxx+H(3,2)*Uyy+0.5*H(3,3)*(Uxy+Uyx)*sqrt(2))/sqrt(2);
                                Eyx=Exy;
                                Exz=sparse(numel(xg),1);
                                Eyz=sparse(numel(xg),1);
                                Ezx=sparse(numel(xg),1);
                                Ezy=sparse(numel(xg),1);
                                Ezz=sparse(numel(xg),1);
                                vm=sqrt(Exx.^2+Eyy.^2+Ezz.^2 ...
                                    -Exx.*Ezz-Eyy.*Ezz-Exx.*Eyy ...
                                    +3*(Exy.^2+Eyz.^2+Exz.^2));
                                
                                
                                
                                count = fprintf(fwid,['SCALARS VM double, 1\n']);
                                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                                if set.ascii
                                    fprintf(fwid, '%f\n', full(vm'));
                                else
                                    fwrite(fwid, full(vm),'double');
                                end
                                count = fprintf(fwid,['TENSORS Stress double\n']);
                            case 'HS'
                                go=1;
                                A=model1.material_parameters.A;
                                Exx=A(1,1)*Uxxx+A(1,2)*Uxxy+(0.5*A(1,3))*(Uxxy+Uyxx)+(0.5*A(1,4))*(Uxyy+Uyxy)+A(1,5)*Uyxy+A(1,6)*Uyyy;
                                Exy=A(2,1)*Uxxx+A(2,2)*Uxxy+(0.5*A(2,3))*(Uxxy+Uyxx)+(0.5*A(2,4))*(Uxyy+Uyxy)+A(2,5)*Uyxy+A(2,6)*Uyyy;
                                Exz=A(3,1)*Uxxx+A(3,2)*Uxxy+(0.5*A(3,3))*(Uxxy+Uyxx)+(0.5*A(3,4))*(Uxyy+Uyxy)+A(3,5)*Uyxy+A(3,6)*Uyyy;
                                Eyx=A(4,1)*Uxxx+A(4,2)*Uxxy+(0.5*A(4,3))*(Uxxy+Uyxx)+(0.5*A(4,4))*(Uxyy+Uyxy)+A(4,5)*Uyxy+A(4,6)*Uyyy;
                                Eyy=A(5,1)*Uxxx+A(5,2)*Uxxy+(0.5*A(5,3))*(Uxxy+Uyxx)+(0.5*A(5,4))*(Uxyy+Uyxy)+A(5,5)*Uyxy+A(5,6)*Uyyy;
                                Eyz=A(6,1)*Uxxx+A(6,2)*Uxxy+(0.5*A(6,3))*(Uxxy+Uyxx)+(0.5*A(6,4))*(Uxyy+Uyxy)+A(6,5)*Uyxy+A(6,6)*Uyyy;
                                Ezx=sparse(numel(xg),1);
                                Ezy=sparse(numel(xg),1);
                                Ezz=sparse(numel(xg),1);
                                
                                
                                count = fprintf(fwid,['TENSORS HStress double\n']);
                        end
                        if go
                            E=[ Exx,Exy,Exz,...
                                Eyx,Eyy,Eyz,...
                                Ezx,Ezy,Ezz]';
                            if set.ascii
                                fprintf(fwid, '%e %e %e \n', full(E));
                            else
                                fwrite(fwid, full(E),'double');
                            end
                        end
                        
                    end
                end
                
                fclose(fwid);
        end
        
    end
    
    
    
end
fprintf(1,'vtkexport done in %5.3f s\n',toc);
