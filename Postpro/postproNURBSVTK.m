function postproNURBSVTK(filreso,submean)
if nargin<2,submean=0;end
nmod=1;
[pp,filres,ext]=fileparts(filreso);
if isempty(ext) filreso=[filreso,'.mat'];end
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
%%
set.ascii=0;
set.remark=' computed by MIC';
nn=prod(Nnodes);
ne=length(elt);
neo=length(elt);
foundt3=find(elt==3);
foundq4=find(elt==4);
foundt=[foundt3;foundq4];
if strcmp(param.analysis,'correlation')
    if exist([filres,'-error.mat'],'file')
        fide=fopen([filres,'-error.mat']);
    else
        fide=fopen(fullfile('TMP',[num2str(nmod),'_error_0.mat']));
    end
    erroronelt=fread(fide,1);
    dynamic=fread(fide,1);
end
if ~strcmp(param.analysis,'mechanics')
    if ~erroronelt
        load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix','wdetJ');
        phix=phix(:,1:size(phix,2)/2);
        M=phix'*(wdetJ*phix);
        phix=wdetJ*phix;
    end
end


switch model1.basis
    case 'nurbs'
        [dphidx,dphidy]=CreateGradNURBSBasis(filres,model1.degree,'nodes',1,'physical');
        [phin]=CreateNURBSBasis(filres,model1.degree,'nodes');
        if isfield(model1,'vtk_export')
            [phig,xg,yg,wg]=CreateNURBSBasis(filres,model1.degree,'Gauss_points');
            [dphidxg,dphidyg,xg,yg,wg]=CreateGradNURBSBasis(filres,model1.degree,'Gauss_points',1,'physical');
            P=phig'*(wg*phig);
            phig=wg*phig;
        end
    case 'btri'
        [phin]=CreateBezierTriangleBasis(filres,'nodes');
        [phig,xg,yg,wg]=CreateBezierTriangleBasis(filres,'Gauss_points');
        [dphidx,dphidy]=CreateGradBezierTriangleBasis(filres,'Gauss_points');
        P=phig'*(wg*phig);
        phig=wg*phig;
end
phi0=sparse(size(dphidx,1),size(dphidx,2));
Uxx=[dphidx,phi0];
Uxy=[dphidy,phi0];
Uyx=[phi0,dphidx];
Uyy=[phi0,dphidy];

nt3=sum(elt==3);
nq4=sum(elt==4);

delete([fullfile('VTK',filres),'-0*.vtk']);
images=[0,images];
U=[zeros(size(U,1),1),U];
UP=[zeros(size(UP,1),1),UP];
L=[1+0*xo,0*xo,yo;...
    0*yo,1+0*yo,-xo];

for iim=1:size(U,2)
    set.vtkname=[filres , sprintf('-%06d.vtk',images(iim))];
    %%
    %ewid = fopen(set.vtkname,'w','b'); % IMPORTANT: big endian
    if set.ascii
    fwid = fopen(fullfile('VTK',set.vtkname),'w'); % IMPORTANT: big endian
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
    count = fprintf(fwid,'POINTS %u double\n',nn);
    %data=[xo+roi(1)-1;yo+roi(3)-1;0*yo]';
    data=[(xo-1)*psample+roi(1),(yo-1)*psample+roi(3),0*yo]';
    
    
    if set.ascii
        fprintf(fwid, '%e %e %e \n', dfac*data);
    else
        fwrite(fwid, data,'double');
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
    Ui=U(1:(2*length(xo)),iim);
    UPi=UP(:,iim);
    if submean
        A=L\Ui;
        Ui=Ui-L*A;
    end
    count = fprintf(fwid,'POINT_DATA %u\n',nn);
    count = fprintf(fwid,['VECTORS Displacement double\n']);
    UU=fac*[Ui(1:nn)';Ui(nn+(1:nn))';repmat(0,1,nn)];
    if set.ascii
        fprintf(fwid, '%e %e %e \n', UU);
    else
        fwrite(fwid,UU,'double');
    end
    if ~strcmp(param.analysis,'mechanics')
        if erroronelt
            if images(iim)==0
                disc=zeros(prod(Nbsnodes),1);
            else
                disc=fread(fide,prod(Nbsnodes));
            end
        else
            if images(iim)==0
                disc=zeros(size(phix,1),1);
            else
                disc=fread(fide,size(phix,1));
            end
            disc=M\(phix'*(disc));
        end
        disc=100*abs(phin*disc)/dynamic;
        count = fprintf(fwid,['SCALARS Error double, 1\n']);
        count = fprintf(fwid,'LOOKUP_TABLE default\n');
        if set.ascii
            fprintf(fwid, '%e\n', full(disc));
        else
            fwrite(fwid,disc,'double');
        end
    end
    count = fprintf(fwid,['TENSORS Strain double\n']);
    Exx=Uxx*UPi;
    Eyy=Uyy*UPi;
    Exy=0.5*(Uxy+Uyx)*UPi;
    if strcmp(model.basis,'btri')
        En=phin*(P\(phig'*[Exx,Eyy,Exy]));
        Exx=En(:,1);Eyy=En(:,2);Exy=En(:,3);
    end
    
    E0=0*Exx;
    E=[Exx,Exy,E0,Exy,Eyy,E0,E0,E0,E0]';
    if set.ascii
        fprintf(fwid, '%e %e %e \n', E);
    else
        fwrite(fwid, E,'double');
    end
    count = fprintf(fwid,['TENSORS Green-Lagrange double\n']);
    Fxx=1+Uxx*UPi;Fxy=Uxy*UPi;
    Fyy=1+Uyy*UPi;Fyx=Uyx*UPi;
    FTFxx=Fxx.*Fxx+Fyx.*Fyx;
    FTFxy=Fxx.*Fxy+Fyx.*Fyy;
    FTFyy=Fxy.*Fxy+Fyy.*Fyy;
    if strcmp(model.basis,'btri')
        En=phin*(P\(phig'*[Fxx,Fyy,FTFxx,FTFxy,FTFyy]));
        Fxx=En(:,1);Fyy=En(:,2);FTFxx=En(:,3);FTFxy=En(:,4);FTFyy=En(:,5);
    end
    E0=0*Exx;
    EGL=0.5*[(FTFxx-1),FTFxy,E0,FTFxy,(FTFyy-1),E0,E0,E0,E0]';
    
    if set.ascii
        fprintf(fwid, '%e %e %e \n', EGL);
    else
        fwrite(fwid, EGL,'double');
    end
%     count = fprintf(fwid,['TENSORS Hencky double\n']);
%     if iim==1
%         H=0*EGL;
%     else
%         H=2*EGL;
%         H(1,:)=H(1,:)+1;
%         H(3,:)=1;
%         H(5,:)=H(5,:)+1;
%         H(6,:)=1;
%         H(7,:)=1;
%         H(8,:)=1;
%         H(9,:)=H(9,:)+1;
%         H=0.5*log(H);
%     end
%     if set.ascii
%         fprintf(fwid, '%e %e %e \n', H);
%     else
%         fwrite(fwid, H,'double');
%     end
    
    if isfield(model1,'vtk_export')
        
        for iex=1:length(model1.vtk_export)
            go=0;
            switch model1.vtk_export{iex}
                case 'S'
                    go=1;
                    if iim==1
                        S=sparse(prod(Nnodes),6);
                        vm=sparse(prod(Nnodes),1);
                    else
                        load(sprintf('%s_%04d',param.result_file,iim-1),'S');
                        Sn=P\(phig'*S);
                        S=phin*Sn;
                        vm=sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2 ...
                            -S(:,1).*S(:,3)-S(:,2).*S(:,3)-S(:,1).*S(:,2) ...
                            +3*(S(:,4).^2+S(:,5).^2+S(:,6).^2));
                    end
                    
                    
                    count = fprintf(fwid,['SCALARS VM double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%e\n', full(vm'));
                    else
                        fwrite(fwid, vm,'double');
                    end
                    
                    
                    count = fprintf(fwid,['TENSORS Stress double\n']);
                case 'EP'
                    go=1;
                    if iim==1
                        S=sparse(prod(Nnodes),6);
                    else
                        load(sprintf('%s_%04d',param.result_file,iim-1),'Ep');
                        Sn=P\(phig'*Ep);
                        S=phin*Sn;
                    end
                    count = fprintf(fwid,['TENSORS EP double\n']);
                    
                    
            end
            if go
                Exx=S(:,1);Eyy=S(:,2);Ezz=S(:,3);
                Exy=S(:,4);Eyz=S(:,5);Exz=S(:,6);
                E=[ Exx,Exy,Exz,...
                    Exy,Eyy,Eyz,...
                    Exz,Eyz,Ezz]';
                if set.ascii
                    fprintf(fwid, '%e %e %e \n', full(E));
                else
                    fwrite(fwid, E,'double');
                end
            end
        end
        
        for iex=1:length(model.vtk_export)
            go=0;
            switch model.vtk_export{iex}
                case 'PEEQ'
                    go=1;
                    load(sprintf('%s_%04d',param.result_file,iim),'Eeqp');
                    Sn=P\(phig'*Eeqp);
                    S=phin*Sn;
                    count = fprintf(fwid,['SCALARS PEEQ double, 1\n']);
            end
            if go
                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                if set.ascii
                    fprintf(fwid, '%e\n', full(S));
                else
                    fwrite(fwid,S,'double');
                end
            end
        end
        
    end
    
    
    
    fclose(fwid);
end
if ~strcmp(param.analysis,'mechanics')
    fclose(fide);
end
fprintf(1,'vtkexport done in %5.3f s\n',toc);
