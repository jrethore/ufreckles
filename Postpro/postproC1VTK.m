function postproC1VTK(filreso,submean,do_error)
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
if iscell(conn)&&model.phantom_nodes==1
postproXC1VTK(filreso)
return
end
if ~exist('U','var')
    U=U1;
end
if ~exist('model1','var')
    model1=model;
end
if strcmp(model1.basis,'nurbs')||strcmp(model1.basis,'btri')
    postproNURBSVTK(filreso,submean);
    return
end
if ~exist('unmasked_nodes','var')
    unmasked_nodes=[];
else
    if ~isempty(unmasked_nodes)
        Nnodes=[length(unmasked_nodes),1,1];
    end
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
nn=prod(Nnodes);
ne=length(elt);
neo=length(elt);
%load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix');
if isfield(model1,'nb_gauss_points')
    ng=model1.nb_gauss_points;
else
    ng=0;
end
if do_error&&~strcmp(param.analysis,'mechanics')
    if exist([filres,'-error',ext],'file')
        fide=fopen([filres,'-error',ext]);
    else
        fide=fopen(fullfile('TMP',[num2str(nmod),'_error_0.mat']));
    end
    erroronelt=fread(fide,1);
    dynamic=fread(fide,1);
    
    if ~erroronelt
        phix=CreateFiniteElementBasis(filreso,sizeim,1,[]);
        [wdetJ,inde]=GetWeigthDetJ(filreso,sizeim);
        M=phix'*wdetJ*phix;
        phix=wdetJ*phix;
        if ng>0
            load(fullfile('TMP',[num2str(nmod),'_phix_0']),'inde');
            gpc2cell=sparse(inde,1:length(inde),diag(wdetJ),ne,length(inde));
            acell=sum(gpc2cell,2);
        end
    end
end

clear wdetJ
rint=false;
if isfield(model,'reduced_integration')
    rint=model.reduced_integration;
end
ngq=4;
if rint,ngq=1;end

save(filreso,'rint','-append');
foundt3=find(elt==3);
foundq4=find(elt==4);
foundt=[foundt3;foundq4];

if ~strcmp(param.analysis,'mechanics')&&do_error
    if ~erroronelt
        if ng>0
            gpc2cell=selected_elt*gpc2cell;
            acell=sum(gpc2cell,2);
        end
    end
end
nt3=sum(elt==3);
nq4=sum(elt==4);

display(sprintf('Result file : %s',filres));
delete([fullfile('VTK',filres),'-0*.vtk']);
L=[1+0*xo,0*xo,yo;...
    0*yo,1+0*yo,-xo];
images=[0,images];
U=[zeros(size(U,1),1),U];
Ux=U(1+6*((1:prod(Nnodes))-1),:);
Uxx=U(2+6*((1:prod(Nnodes))-1),:);
Uxy=U(3+6*((1:prod(Nnodes))-1),:);
Uxxx=U(4+6*((1:prod(Nnodes))-1),:);
Uxxy=U(5+6*((1:prod(Nnodes))-1),:);
Uxyy=U(6+6*((1:prod(Nnodes))-1),:);
Uy=U(6*prod(Nnodes)+1+6*((1:prod(Nnodes))-1),:);
Uyx=U(6*prod(Nnodes)+2+6*((1:prod(Nnodes))-1),:);
Uyy=U(6*prod(Nnodes)+3+6*((1:prod(Nnodes))-1),:);
Uyxx=U(6*prod(Nnodes)+4+6*((1:prod(Nnodes))-1),:);
Uyxy=U(6*prod(Nnodes)+5+6*((1:prod(Nnodes))-1),:);
Uyyy=U(6*prod(Nnodes)+6+6*((1:prod(Nnodes))-1),:);
for iim=1:size(U,2)
    set.vtkname=[filres , sprintf('-%06d.vtk',images(iim))];
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
    %data=[xo+roi(1)-1;yo+roi(3)-1;0*yo]';
    data=[(xo-1)*psample+roi(1),(yo-1)*psample+roi(3),0*yo]';
    
    
    if set.ascii
        fprintf(fwid, '%f %f %f \n', dfac*data);
    else
        fwrite(fwid, data,'float');
    end
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
    %
    count = fprintf(fwid,'POINT_DATA %u\n',nn);
    count = fprintf(fwid,['VECTORS Displacement float\n']);
    UU=fac*[Ux(:,iim)';Uy(:,iim)';repmat(0,1,nn)];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
    
    
    count = fprintf(fwid,['TENSORS Strain float\n']);
    E0=zeros(nn,1);
    E=[Uxx(:,iim),0.5*(Uxy(:,iim)+Uyx(:,iim)),E0,0.5*(Uxy(:,iim)+Uyx(:,iim)),Uyy(:,iim),E0,E0,E0,E0]';
    if set.ascii
        fprintf(fwid, '%f %f %f\n', E);
    else
        fwrite(fwid,E,'float');
    end
    
    count = fprintf(fwid,['TENSORS DStrain float\n']);
    E0=zeros(nn,1);
    E=[Uxxx(:,iim),Uxxy(:,iim),0.5*(Uxxy(:,iim)+Uyxx(:,iim)),0.5*(Uxyy(:,iim)+Uyxy(:,iim)),Uyxy(:,iim),Uyyy(:,iim),E0,E0,E0]';
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
                    if iim==1
                        Exx=sparse(nn,1);
                        Exy=sparse(nn,1);
                        Exz=sparse(nn,1);
                        Eyx=sparse(nn,1);
                        Eyy=sparse(nn,1);
                        Eyz=sparse(nn,1);
                        Ezx=sparse(nn,1);
                        Ezy=sparse(nn,1);
                        Ezz=sparse(nn,1);
                        vm=sparse(nn,1);
                    else
                        Exx=H(1,1)*Uxx(:,iim)+H(1,2)*Uyy(:,iim)+0.5*H(1,3)*(Uxy(:,iim)+Uyx(:,iim))*sqrt(2);
                        Eyy=H(2,1)*Uxx(:,iim)+H(2,2)*Uyy(:,iim)+0.5*H(2,3)*(Uxy(:,iim)+Uyx(:,iim))*sqrt(2);
                        Exy=(H(3,1)*Uxx(:,iim)+H(3,2)*Uyy(:,iim)+0.5*H(3,3)*(Uxy(:,iim)+Uyx(:,iim))*sqrt(2))/sqrt(2);
                        Eyx=Exy;
                        Exz=sparse(nn,1);
                        Eyz=sparse(nn,1);
                        Ezx=sparse(nn,1);
                        Ezy=sparse(nn,1);
                        Ezz=sparse(nn,1);
                        vm=sqrt(Exx.^2+Eyy.^2+Ezz.^2 ...
                            -Exx.*Ezz-Eyy.*Ezz-Exx.*Eyy ...
                            +3*(Exy.^2+Eyz.^2+Exz.^2));
                    end
                    
                    
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
                    if iim==1
                        Exx=sparse(nn,1);
                        Exy=sparse(nn,1);
                        Exz=sparse(nn,1);
                        Eyx=sparse(nn,1);
                        Eyy=sparse(nn,1);
                        Eyz=sparse(nn,1);
                        Ezx=sparse(nn,1);
                        Ezy=sparse(nn,1);
                        Ezz=sparse(nn,1);
                    else
                        Exx=A(1,1)*Uxxx(:,iim)+A(1,2)*Uxxy(:,iim)+(0.5*A(1,3))*(Uxxy(:,iim)+Uyxx(:,iim))+(0.5*A(1,4))*(Uxyy(:,iim)+Uyxy(:,iim))+A(1,5)*Uyxy(:,iim)+A(1,6)*Uyyy(:,iim);
                        Exy=A(2,1)*Uxxx(:,iim)+A(2,2)*Uxxy(:,iim)+(0.5*A(2,3))*(Uxxy(:,iim)+Uyxx(:,iim))+(0.5*A(2,4))*(Uxyy(:,iim)+Uyxy(:,iim))+A(2,5)*Uyxy(:,iim)+A(2,6)*Uyyy(:,iim);
                        Exz=A(3,1)*Uxxx(:,iim)+A(3,2)*Uxxy(:,iim)+(0.5*A(3,3))*(Uxxy(:,iim)+Uyxx(:,iim))+(0.5*A(3,4))*(Uxyy(:,iim)+Uyxy(:,iim))+A(3,5)*Uyxy(:,iim)+A(3,6)*Uyyy(:,iim);
                        Eyx=A(4,1)*Uxxx(:,iim)+A(4,2)*Uxxy(:,iim)+(0.5*A(4,3))*(Uxxy(:,iim)+Uyxx(:,iim))+(0.5*A(4,4))*(Uxyy(:,iim)+Uyxy(:,iim))+A(4,5)*Uyxy(:,iim)+A(4,6)*Uyyy(:,iim);
                        Eyy=A(5,1)*Uxxx(:,iim)+A(5,2)*Uxxy(:,iim)+(0.5*A(5,3))*(Uxxy(:,iim)+Uyxx(:,iim))+(0.5*A(5,4))*(Uxyy(:,iim)+Uyxy(:,iim))+A(5,5)*Uyxy(:,iim)+A(5,6)*Uyyy(:,iim);
                        Eyz=A(6,1)*Uxxx(:,iim)+A(6,2)*Uxxy(:,iim)+(0.5*A(6,3))*(Uxxy(:,iim)+Uyxx(:,iim))+(0.5*A(6,4))*(Uxyy(:,iim)+Uyxy(:,iim))+A(6,5)*Uyxy(:,iim)+A(6,6)*Uyyy(:,iim);
                        Ezx=sparse(nn,1);
                        Ezy=sparse(nn,1);
                        Ezz=sparse(nn,1);
                    end
                    
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
    
    
    
    
    if exist('cn','var')
        count = fprintf(fwid,'CELL_DATA %u\n',ne);
        count = fprintf(fwid,'SCALARS Eroded float, 1\n');
        count = fprintf(fwid,'LOOKUP_TABLE default\n');
        if set.ascii
            fprintf(fwid, '%f\n', cn(:,max(1,iim-1)));
        else
            fwrite(fwid,cn(:,max(1,iim-1)),'float');
        end
    end
    fclose(fwid);
end
fprintf(1,'vtkexport done in %5.3f s\n',toc);
