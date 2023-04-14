function postproCSV3D(filreso,submean,do_error)
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
[phi]=CreateFiniteElementBasis3D(filreso,[1,1],1,[],'Gauss_points');
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

xc=gp2cell*(phi*xo);
yc=gp2cell*(phi*yo);
zc=gp2cell*(phi*zo);



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
delete([filres,'-0*.vtk']);
switch submean
    case {0,1}
L=[1+0*xo,0*xo,0*xo,0*xo,-zo,yo;...
    0*yo,1+0*yo,0*yo,zo,0*yo,-xo;...
    0*zo,0*zo,1+0*zo,-yo,xo,0*zo];
    case 2
 L=[1+0*xo,0*xo,0*xo,0*xo,-zo,yo,xo;...
    0*yo,1+0*yo,0*yo,zo,0*yo,-xo,yo;...
    0*zo,0*zo,1+0*zo,-yo,xo,0*zo,zo];
       
end
display(sprintf('Result file : %s',filres));
images=[images];
for iim=1:size(U,2)    
E=[xo+roi(1)-1,yo+roi(3)-1,zo+roi(5)-1,...
    U(0*prod(Nnodes)+(1:prod(Nnodes)),iim),U(1*prod(Nnodes)+(1:prod(Nnodes)),iim),U(2*prod(Nnodes)+(1:prod(Nnodes)),iim)];
filname=[ filres , sprintf('-disp-%05d.csv',images(iim))];
   dlmwrite(filname, ['X','Y','Z','U','V','W'],';');
dlmwrite(filname,E,'-append','delimiter',';');
Ui=full(U(1:3*nn,iim));

 Exx=epsxx*Ui;  
 Eyy=epsyy*Ui;  
 Ezz=epszz*Ui;  
 Exy=0.5*epsxy*Ui;
Exz=0.5*epsxz*Ui;
Eyz=0.5*epsyz*Ui;


 E=[Exx,Eyy,Ezz,Exy,Exz,Eyz];
filname=[ filres , sprintf('-strain-%05d.csv',images(iim))];
   dlmwrite(filname, ['X','Y','Z','Exx','Eyy','Ezz','Exy','Exz','Eyz'],';');
dlmwrite(filname,E,'-append','delimiter',';');
error


if isfield(model,'vtk_export')

    for iex=1:length(model.vtk_export)
        go=0;
        switch model.vtk_export{iex}
            case 'S'
                go=1;
                    load(sprintf('%s_%04d',param.result_file,iim),'S');
                    S=gp2cell*S;
                Exx=S(:,1);Eyy=S(:,2);Ezz=S(:,3);
                Exy=S(:,4);Eyz=S(:,5);Exz=S(:,6);
filname=[ filres , sprintf('-stress-%05d.csv',images(iim))];
   dlmwrite(filname, ['X','Y','Z','Sxx','Syy','Szz','Sxy','Sxz','Syz'],';');

        end
        if go
        E=[ xc,yc,zc,Exx,Eyy,Ezz,Exy,Exz,Eyz];
dlmwrite(filname,E,'-append','delimiter',';');
        end
    end

end




end
fprintf(1,'csvexport done in %5.3f s\n',toc);
