function postproVTK3D(filreso,submean,do_error)
if nargin<2,submean=1;end
if nargin<3,do_error=1;end
    nmod=1;

tic

[pp,filres,ext]=fileparts(filreso);
if isempty(ext)
    filreso=[filreso,'.mat'];
ext='.mat';
end
tic
load(filreso,'-mat')


if strcmp(model.basis,'nurbs')
    postproNURBSVTK3D(filreso,submean);
    return
end
if isfield(param,'image_number')
    images=param.image_number;
else
images=1:size(U,2);
end
roi=param.roi;
if isfield(param,'calibration_data')
    roi=0*roi+1;
end
if length(roi)<5
    roi(5)=1;
end
%%
set.ascii=0;
set.remark=' computed by MIC';
nn=prod(Nnodes);
ne=length(elt);

[dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis3D(filreso,[1,1],1,[],'Gauss_points');
phi0=sparse(size(dphidx,1),size(dphidx,2));
epsxx=[dphidx,phi0,phi0];
epsyy=[phi0,dphidy,phi0];
epsxy=[dphidy,dphidx,phi0];
Uxy=[dphidy,phi0,phi0];
Uyx=[phi0,dphidx,phi0];
epsxz=[dphidz,phi0,dphidx];
Uxz=[dphidz,phi0,phi0];
Uzx=[phi0,phi0,dphidx];
epsyz=[phi0,dphidz,dphidy];
Uyz=[phi0,dphidz,phi0];
Uzy=[phi0,phi0,dphidy];
epszz=[phi0,phi0,dphidz];


% load(fullfile('TMP',[num2str(nmod),'_3d_epsxx_0']),'epsxx');
% load(fullfile('TMP',[num2str(nmod),'_3d_epsyy_0']),'epsyy');
% load(fullfile('TMP',[num2str(nmod),'_3d_epszz_0']),'epszz');
% load(fullfile('TMP',[num2str(nmod),'_3d_epsxy_0']),'Uyx','Uxy','epsxy');
% load(fullfile('TMP',[num2str(nmod),'_3d_epsxz_0']),'Uxz','Uzx','epsxz');
% load(fullfile('TMP',[num2str(nmod),'_3d_epsyz_0']),'Uyz','Uzy','epsyz');


rint=0;
if isfield(model,'reduced_integration')
    rint=model.reduced_integration;
end
ngq=1+7*(~rint);
ngp=1+(~rint);
ngt=1;
 foundt3=find(elt==6);
 foundt4=find(elt==4);
 foundq4=find(elt==8);
 foundt=[foundt3;foundt4;foundq4];
indi=repmat((1:length(elt))',1,max(ngp,ngq));
val=[repmat((elt==8),1,ngq)/ngq,repmat(0,length(elt),max(0,ngp+ngt-ngq))]+...
    [repmat((elt==6),1,ngp),repmat(0,length(elt),max(0,ngq-ngp))]+...
    [repmat((elt==4),1,ngt),repmat(0,length(elt),max(0,ngq-ngt))];
indj=ones(size(indi'));
found=find(val'>0);
indj(found)=1:length(found);
indj=indj';
 val=val(foundt,:);
 indi=indi(foundt,:);
 indi=indi(foundt,:);
gp2cell=sparse(indi,indj,val);
selected_elt=1;
neo=ne;
if exist('face_elts')
    selected_elt=ones(length(elt),1);
    selected_elt(face_elts)=0;
    indj=find(selected_elt);
    indi=1:length(indj);
    selected_elt=sparse(indi,indj,1,length(indj),length(elt));
    elt=full(selected_elt*elt);
    conn=full(selected_elt*conn);
    foundt3=find(elt==6);
   foundt4=find(elt==4);
    foundq4=find(elt==8);
    ne=length(elt);
    postproVTKX3D(filres,submean);
end

    gp2cell=selected_elt*gp2cell;
    epsxx=gp2cell*epsxx;
epsxy=gp2cell*epsxy;
epsyy=gp2cell*epsyy;
epszz=gp2cell*epszz;
epsxz=gp2cell*epsxz;
epsyz=gp2cell*epsyz;

Uxx=epsxx;
Uxy=gp2cell*Uxy;
Uxz=gp2cell*Uxz;
Uyx=gp2cell*Uyx;
Uyy=epsyy;
Uyz=gp2cell*Uyz;
Uzx=gp2cell*Uzx;
Uzy=gp2cell*Uzy;
Uzz=epszz;




nt4=sum(elt==4);
np6=sum(elt==6);
nh8=sum(elt==8);
if do_error&&strcmp(param.analysis,'correlation')&&(~isfield(model,'extrusion_parameters'))
    if exist([filres,'-error',ext],'file')
    fide=fopen([filres,'-error',ext]);
    else
    fide=fopen(fullfile('TMP',[num2str(nmod),'_error_0']));
    end
    dynamic=fread(fide,1);
end
delete([fullfile('VTK',filres),'-0*.vtk']);
switch submean
    case {0,1}
L=[1+0*xo,0*xo,0*xo,0*xo,-zo,yo;...
    0*yo,1+0*yo,0*yo,zo,0*yo,-xo;...
    0*zo,0*zo,1+0*zo,-yo,xo,0*zo];
    case 2
 L=[1+0*xo,0*xo,0*xo,0*xo,-zo,yo,xo;...
    0*yo,1+0*yo,0*yo,zo,0*yo,-xo,yo;...
    0*zo,0*zo,1+0*zo,-yo,xo,0*zo,zo];
    case 3
 L=[1+0*xo,0*xo,0*xo,0*xo,-zo,yo,xo,0*xo,0*xo;...
    0*yo,1+0*yo,0*yo,zo,0*yo,-xo,0*yo,yo,0*yo;...
    0*zo,0*zo,1+0*zo,-yo,xo,0*zo,0*zo,0*zo,zo];
       
end
display(sprintf('Result file : %s',filres));
images=[0,images];
U=[zeros(size(U,1),1),U];
for iim=1:size(U,2)    
set.vtkname=[ filres , sprintf('-%06d.vtk',images(iim))];
%%
%fwid = fopen(set.vtkname,'w','b'); % IMPORTANT: big endian
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
data=[xo+roi(1)-1,yo+roi(3)-1,zo+roi(5)-1]';


if set.ascii
   fprintf(fwid, '%e %e %e \n', data);
else  
    fwrite(fwid, data,'double');
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
% 
Ui=full(U(1:3*nn,iim));
 count = fprintf(fwid,'POINT_DATA %u\n',nn);
 count = fprintf(fwid,['VECTORS Displacement double\n']);
 UU=[Ui(1:nn)';Ui(nn+(1:nn))';Ui(2*nn+(1:nn))'];
if set.ascii
   fprintf(fwid, '%e %e %e \n', UU);
else
 fwrite(fwid,UU,'double');     
end
    A=L\Ui;
    Ui=Ui-L*A;
 count = fprintf(fwid,['VECTORS U double\n']);
 UU=[Ui(1:nn)';Ui(nn+(1:nn))';Ui(2*nn+(1:nn))'];
if set.ascii
   fprintf(fwid, '%e %e %e \n', UU);
else
 fwrite(fwid,UU,'double');     
end


  count = fprintf(fwid,'CELL_DATA %u\n',ne);
 count = fprintf(fwid,['TENSORS Strain double\n']);
 Exx=epsxx*Ui;  
 Eyy=epsyy*Ui;  
 Ezz=epszz*Ui;  
 Exy=0.5*epsxy*Ui;
Exz=0.5*epsxz*Ui;
Eyz=0.5*epsyz*Ui;


 E=[ Exx,Exy,Exz,...
     Exy,Eyy,Eyz,...
     Exz,Eyz,Ezz]';

if set.ascii
   fprintf(fwid, '%e %e %e \n', E);
else
      fwrite(fwid, E,'double');
end

    count = fprintf(fwid,['TENSORS Green-Lagrange float\n']);
    Fxx=1+Uxx*Ui;Fxy=Uxy*Ui;Fxz=Uxz*Ui;
    Fyy=1+Uyy*Ui;Fyx=Uyx*Ui;Fyz=Uyz*Ui;
    Fzz=1+Uzz*Ui;Fzx=Uzx*Ui;Fzy=Uzy*Ui;
    FTFxx=Fxx.*Fxx+Fyx.*Fyx+Fzx.*Fzx;
    FTFxy=Fxx.*Fxy+Fyx.*Fyy+Fzx.*Fzy;
    FTFxz=Fxx.*Fxz+Fyx.*Fyz+Fzx.*Fzz;
    FTFyy=Fxy.*Fxy+Fyy.*Fyy+Fzy.*Fzy;
    FTFyz=Fxy.*Fxz+Fyy.*Fyz+Fzy.*Fzz;
    FTFzz=Fxz.*Fxz+Fyz.*Fyz+Fzz.*Fzz;
    EGL=0.5*[(FTFxx-1),FTFxy,FTFxz,FTFxy,(FTFyy-1),FTFyz,FTFxz,FTFyz,(FTFzz-1)]';

    if set.ascii
        fprintf(fwid, '%f %f %f \n', EGL);
    else
        fwrite(fwid, EGL,'float');
    end
%     count = fprintf(fwid,['TENSORS Hencky float\n']);
%     if iim==1
%         H=0*EGL;
%     else
%         H=2*EGL;
%         H(1,:)=H(1,:)+1;
%         H(5,:)=H(5,:)+1;
%         H(9,:)=H(9,:)+1;
%         H=0.5*log(H);
%     end
%     if set.ascii
%         fprintf(fwid, '%f %f %f \n', H);
%     else
%         fwrite(fwid, H,'float');
%     end






if isfield(model,'vtk_export')

    for iex=1:length(model.vtk_export)
        go=0;
        switch model.vtk_export{iex}
            case 'S'
                go=1;
                if iim==1
                    S=sparse(size(gp2cell,1),6);
                    vm=sparse(size(gp2cell,1),1);
                else
                    load(sprintf('%s_%04d',filres,iim-1),'S');
                    vm=gp2cell*sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2 ...
                        -S(:,1).*S(:,3)-S(:,2).*S(:,3)-S(:,1).*S(:,2) ...
                        +3*(S(:,4).^2+S(:,5).^2+S(:,6).^2));
                    S=gp2cell*S;
                end
                    count = fprintf(fwid,['SCALARS VM double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%f\n', full(vm'));
                    else
                        fwrite(fwid, full(vm),'double');
                    end
                Exx=S(:,1);Eyy=S(:,2);Ezz=S(:,3);
                Exy=S(:,4);Eyz=S(:,5);Exz=S(:,6);
                count = fprintf(fwid,['TENSORS Stress double\n']);
            case 'LE'
                go=1;
                if iim==1
                    S=sparse(size(gp2cell,1),6);
                else
                    load(sprintf('%s_%04d',param.result_file,iim-1),'E');
                    S=gp2cell*E;
                end
                Exx=S(:,1);Eyy=S(:,2);Ezz=S(:,3);
                Exy=S(:,4);Eyz=S(:,5);Exz=S(:,6);
                count = fprintf(fwid,['TENSORS LE double\n']);
            case 'EP'
                go=1;
                if iim==1
                    S=sparse(size(gp2cell,1),6);
                else
                    load(sprintf('%s_%04d',param.result_file,iim-1),'Ep');
                    S=gp2cell*Ep;
                end
                Exx=S(:,1);Eyy=S(:,2);Ezz=S(:,3);
                Exy=S(:,4);Eyz=S(:,5);Exz=S(:,6);
                count = fprintf(fwid,['TENSORS EP double\n']);

        end
        if go
        E=[ Exx,Exy,Exz,...
            Exy,Eyy,Eyz,...
            Exz,Eyz,Ezz]';
        if set.ascii
            fprintf(fwid, '%e %e %e \n', full(E));
        else
            fwrite(fwid, full(E),'double');
        end
        end
    end

end



if do_error&&strcmp(param.analysis,'correlation')&&(~isfield(model,'extrusion_parameters'))
    if images(iim)==0
        disc=zeros(neo,1);
    else
    disc=fread(fide,neo);
    end
    disc=100*(selected_elt*disc)/dynamic;
 count = fprintf(fwid,['SCALARS Error double, 1\n']);
count = fprintf(fwid,'LOOKUP_TABLE default\n');
if set.ascii
   fprintf(fwid, '%e\n', disc);
else
 fwrite(fwid,disc,'double');     
end
end  
if isfield(model,'vtk_export')

    for iex=1:length(model.vtk_export)
        go=0;
        switch model.vtk_export{iex}
            case 'PEEQ'
                go=1;
                load(sprintf('%s_%04d',param.result_file,iim),'Eeqp');
                S=gp2cell*abs(Eeqp);
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
if do_error&&strcmp(param.analysis,'correlation')&&(~isfield(model,'extrusion_parameters'))
fclose(fide);
end
fprintf(1,'vtkexport done in %5.3f s\n',toc);
