function postproNURBSVTK25D(filres,submean)
if nargin<2,submean=0;end
nmod=1;

tic
load([filres,'.mat'])
if ~exist('model1')
    model1=model;
end
assert(strcmp(model1.basis,'nurbs'));
dotopo=(~isfield(model1,'topography'))&&isfield(param,'calibration_data');

if isfield(param,'image_number')
    images=param.image_number;
    if dotopo
        images(1)=[];
    end
else
    images=1:size(U,2);
end

%%
set.ascii=0;
set.remark=' computed by MIC';
nn=prod(Nnodes);
ne=length(elt);
foundt3=find(elt==3);
foundq4=find(elt==4);
foundt=[foundt3;foundq4];
nt3=sum(elt==3);
nq4=sum(elt==4);

[dphidx,dphidy,dphidz]=CreateGradNURBSBasis25D(filres,model1.degree,'nodes',1,'physical');

phi0=sparse(size(dphidx,1),size(dphidx,2));
Uxx=[dphidx,phi0,phi0];
Uxy=[dphidy,phi0,phi0];
Uxz=[dphidz,phi0,phi0];
Uyx=[phi0,dphidx,phi0];
Uyy=[phi0,dphidy,phi0];
Uyz=[phi0,dphidz,phi0];
Uzx=[phi0,phi0,dphidx];
Uzy=[phi0,phi0,dphidy];
Uzz=[phi0,phi0,dphidz];


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
UP=[zeros(size(UP,1),1),UP];
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
    data=[Xo,Yo,Zo]';
    
    
    if set.ascii
        fprintf(fwid, '%e %e %e \n', data);
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
    Ui=U(1:(3*length(Xo)),iim);
    UPi=UP(:,iim);
    if submean
        A=L\Ui;
        Ui=Ui-L*A;
    end
    count = fprintf(fwid,'POINT_DATA %u\n',nn);
    count = fprintf(fwid,['VECTORS Displacement double\n']);
    UU=[Ui(1:nn)';Ui(nn+(1:nn))';Ui(2*nn+(1:nn))'];
    if set.ascii
        fprintf(fwid, '%e %e %e \n', UU);
    else
        fwrite(fwid,UU,'double');
    end
    
    count = fprintf(fwid,['TENSORS Strain double\n']);
    Exx=Uxx*UPi;
    Eyy=Uyy*UPi;
    Ezz=Uzz*UPi;
    Exy=0.5*(Uxy+Uyx)*UPi;
    Exz=0.5*(Uxz+Uzx)*UPi;
    Eyz=0.5*(Uyz+Uzy)*UPi;
    
    
    E=[ Exx,Exy,Exz,...
        Exy,Eyy,Eyz,...
        Exz,Eyz,Ezz]';
    
    if set.ascii
        fprintf(fwid, '%e %e %e \n', E);
    else
        fwrite(fwid, E,'double');
    end
    count = fprintf(fwid,['TENSORS Green-Lagrange double\n']);
    Fxx=1+Uxx*UPi;Fxy=Uxy*UPi;Fxz=Uxz*UPi;
    Fyy=1+Uyy*UPi;Fyx=Uyx*UPi;Fyz=Uyz*UPi;
    Fzz=1+Uzz*UPi;Fzx=Uzx*UPi;Fzy=Uzy*UPi;
    FTFxx=Fxx.*Fxx+Fyx.*Fyx+Fzx.*Fzx;
    FTFxy=Fxx.*Fxy+Fyx.*Fyy+Fzx.*Fzy;
    FTFxz=Fxx.*Fxz+Fyx.*Fyz+Fzx.*Fzz;
    FTFyy=Fxy.*Fxy+Fyy.*Fyy+Fzy.*Fzy;
    FTFyz=Fxy.*Fxz+Fyy.*Fyz+Fzy.*Fzz;
    FTFzz=Fxz.*Fxz+Fyz.*Fyz+Fzz.*Fzz;
    EGL=0.5*[(FTFxx-1),FTFxy,FTFxz,FTFxy,(FTFyy-1),FTFyz,FTFxz,FTFyz,(FTFzz-1)]';
    
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
%         H(5,:)=H(5,:)+1;
%         H(9,:)=H(9,:)+1;
%         H=0.5*log(H);
%     end
%     if set.ascii
%         fprintf(fwid, '%e %e %e \n', H);
%     else
%         fwrite(fwid, H,'double');
%     end
    
    fclose(fwid);
end

fprintf(1,'vtkexport done in %5.3f s\n',toc);
