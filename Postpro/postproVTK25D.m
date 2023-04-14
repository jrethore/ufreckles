function postproVTK25D(filreso,submean)
if nargin<2,submean=0;end
nmod=1;
[pp,filres,ext]=fileparts(filreso);
if isempty(ext)
    filreso=[filreso,'.mat'];
ext='.mat';
end

tic
load(filreso,'-mat')
if ~exist('model1')
model1=model;
end
if strcmp(model1.basis,'nurbs')
    postproNURBSVTK25D(filres,submean);
    return
end

dotopo=(~isfield(model1,'topography'))&&isfield(param,'calibration_data');

if isfield(param,'image_number')
    images=param.image_number;
    if dotopo
    images(1)=[];
    end
else
images=1:size(U,2);
end
if isfield(param,'pixel_size')
    pix2m=param.pixel_size;
else 
    pix2m=1;
end
if strcmp(param.analysis,'mechanics')
    dfac=floor(log10(pix2m));
    fac=pix2m/dfac;
    fac=pix2m;
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
set.remark=' computed by UFRECKLES';
nn=prod(Nnodes);
ne=length(elt);
 foundt3=find(elt==3);
 foundq4=find(elt==4);
 foundt=[foundt3;foundq4];
nt3=sum(elt==3);
 nq4=sum(elt==4);
% load(fullfile('TMP',[num2str(nmod),'_25d_epsxx_0']),'epsxx');
% load(fullfile('TMP',[num2str(nmod),'_25d_epsyy_0']),'epsyy');
% load(fullfile('TMP',[num2str(nmod),'_25d_epszz_0']),'epszz');
% load(fullfile('TMP',[num2str(nmod),'_25d_epsxy_0']),'Uyx','Uxy','epsxy');
% load(fullfile('TMP',[num2str(nmod),'_25d_epsxz_0']),'Uxz','Uzx','epsxz');
% load(fullfile('TMP',[num2str(nmod),'_25d_epsyz_0']),'Uyz','Uzy','epsyz');
rint=0;
if isfield(model,'reduced_integration')
    rint=model.reduced_integration;
end
ngq=1+3*(~rint);

save(filreso,'rint','-append');
[dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis25D(filreso,[1,1],1,[],'Gauss_points');
        phi0=sparse(size(dphidx,1),size(dphidx,2));
        epsxx=[dphidx,phi0,phi0];
        epsyy=[phi0,dphidy,phi0];
        epsxy=[dphidy,dphidx,phi0];
            Uxy=[dphidy,phi0,phi0];
            Uyx=[phi0,dphidx,phi0];
        epsxz=[phi0,phi0,dphidx];
        Uxz=[phi0,phi0,phi0];
        Uzx=[phi0,phi0,dphidx];
        epsyz=[phi0,phi0,dphidy];
        Uyz=[phi0,phi0,phi0];
        Uzy=[phi0,phi0,dphidy];
        epszz=[phi0,phi0,phi0];


ngp=1;
indi=repmat((1:length(elt))',1,max(ngp,ngq));
val=[repmat((elt==4),1,ngq)/ngq,zeros(length(elt),max(0,ngp-ngq))]+[repmat((elt==3),1,ngp),repmat(0,length(elt),max(0,ngq-ngp))];
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
    foundt3=find(elt==3);
    foundq4=find(elt==4);
    ne=length(elt);
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

np6=sum(elt==3);
nh8=sum(elt==4);
nt3=sum(elt==3);
nq4=sum(elt==4);

unix(['rm ',fullfile('VTK',filres),'-0*.vtk']);
L=[1+0*Xo,0*Xo,0*Xo,0*Xo,-Zo,Yo;...
    0*Yo,1+0*Yo,0*Yo,Zo,0*Yo,-Xo;...
    0*Zo,0*Zo,1+0*Zo,-Yo,Xo,0*Zo];
images=[0,images];
U=[zeros(size(U,1),1),U];
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
count = fprintf(fwid,'POINTS %u float\n',nn);
%data=[xo+roi(1)-1;yo+roi(3)-1;0*yo]';
data=[Xo,Yo,Zo]';


if set.ascii
   fprintf(fwid, '%f %f %f \n', fac*data);
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
Ui=U(1:(3*length(Xo)),iim);
 count = fprintf(fwid,'POINT_DATA %u\n',nn);
 count = fprintf(fwid,['VECTORS Displacement float\n']);
 UU=fac*[Ui(1:nn)';Ui(nn+(1:nn))';Ui(2*nn+(1:nn))'];
if set.ascii
   fprintf(fwid, '%f %f %f \n', UU);
else
 fwrite(fwid,UU,'float');     
end
    A=L\Ui;
    Ui=Ui-L*A;
 count = fprintf(fwid,['VECTORS U float\n']);
 UU=fac*[Ui(1:nn)';Ui(nn+(1:nn))';Ui(2*nn+(1:nn))'];
if set.ascii
   fprintf(fwid, '%f %f %f \n', UU);
else
 fwrite(fwid,UU,'float');     
end
    

  count = fprintf(fwid,'CELL_DATA %u\n',ne);
 count = fprintf(fwid,['TENSORS Strain float\n']);
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
   fprintf(fwid, '%f %f %f \n', E);
else
      fwrite(fwid, E,'float');
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
%         fprintf(fwid, '%e %e %e \n', H);
%     else
%         fwrite(fwid, H,'float');
%     end

fclose(fwid);
end
% if ~strcmp(param.analysis,'mechanics')
% fclose(fide);
% end
fprintf(1,'vtkexport done in %5.3f s\n',toc);
