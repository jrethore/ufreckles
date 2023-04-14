function [H,A,bulk,HH]=ComputeEffectiveProperties(matmod,periodic)
if nargin<2,periodic=0;end
clear functions
load('cell-mesh','n','zp')
[~,indc]=sort(angle(zp));
lbox=max(abs(zp))+10;
zp=zp(indc);
nc=2;
nmod=5;
filres='cell-hom'; 
param.analysis='mechanics';
param.result_file=filres;
param.pixel_size=1e-6;
pix2m=param.pixel_size;
param.roi=[1,2*lbox+10,1,2*lbox+10];

model.basis='fem';
model.mesh_file='cell-mesh.vtk'; 

xc=lbox;
yc=lbox;
model.gluing_parameters{1}={'scale',[0;0;0],1};
model.gluing_parameters{2}={'translate',[xc-1;yc-1;0]+1};


model.material_model='elastic_homogeneous_isotropic';
model.material_parameters=matmod;
model.vtk_export={'S'};
 
iscale=1;
LoadParameters(param);
LoadParameters(model,nmod);
sizeim=[diff(param.roi(1:2)),diff(param.roi(3:4))]+1;
save(fullfile('TMP','sample0_0'),'sizeim')
save(fullfile('TMP','sample0'),'sizeim')

LoadMask(nmod);   % ???
LoadMeshes(nmod); % generate the meshes
LoadMat(nmod);    % define some materials constants

%%
meshfile=fullfile('TMP',sprintf('%d_mesh_%d',nmod,10*(iscale-1)));
load(meshfile);
CreateGradBasisFunction(1,nmod);
AssembleMechanicalOperator(1,nmod);
%%

load(fullfile('TMP',sprintf('%d_k_operator_%d.mat',nmod,10*(iscale-1))),'K');

load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'wdetJ');

phig=CreateFiniteElementBasis(meshfile,sizeim,1,[],'Gauss_points');


xg=phig*((xo-xc)*pix2m);  
yg=phig*((yo-yc)*pix2m);

      zo=(xo-xc)+1i*(yo-yc);
      zp=[zp(:);zp(1)];
    selected=GetSignedDistanceToZone([real(zp),imag(zp)],[min(xo-xc)-10,max(xo-xc)+10,min(yo-xc)-10,max(yo-xc)+10],xo-xc,yo-xc);


    selected=selected<-1;
if periodic
    edges=cell(n/2,2);onedges=0*selected;
    Ce=sparse(prod(Nnodes),n/2);
    for ie=1:n/2
%              figure
%          plot(zo,'x')
%          hold on
%          plot(zp(ie+[0,1]),'ro')
        rot=exp(-1i*(angle(mean(zp(ie+[0,1])))));
%          figure
%          plot(zo*rot,'x')
%          hold on
%          plot(zp(ie+[0,1])*rot,'ro')
         along=find(abs(real(zo*rot-zp(ie)*rot))<2);
        [~,idm]=min(abs(zo(along)-zp(ie)));
        a1=angle(zo(along(idm))*rot);
        [~,idm]=min(abs(zo(along)-zp(ie+1)));
        a2=angle(zo(along(idm))*rot);
% a1=angle(zp(ie)*rot);
% a2=angle(zp(ie+1)*rot);
       
        on1=((angle(zo*rot)>=min(a1,a2))&(angle(rot*zo)<=max(a1,a2)))&~selected&~onedges;
%       plot(zo(on1)*rot,'ms')
on1=find(on1);
[~,id]=sort(angle(zo(on1)*rot));
on1=on1(id);
edges{ie,1}=on1;
Ce(on1,ie)=1;
       rot=exp(-1i*(angle(mean(zp(ie+n/2+[0,1])))));
%          figure
%          plot(zo*rot,'x')
%          hold on
%          plot(zo(on1)*rot,'ms')
%          plot(zp(ie+n/2+[0,1])*rot,'ro')
         
         along=find(abs(real(zo*rot-zp(ie+n/2)*rot))<2);
                 [~,idm]=min(abs(zo(along)-zp(ie+n/2)));
        a1=angle(zo(along(idm))*rot);
        [~,idm]=min(abs(zo(along)-zp(ie+1+n/2)));
        a2=angle(zo(along(idm))*rot);

%   a1=angle(zp(ie+n/2)*rot);
% a2=angle(zp(ie+1+n/2)*rot);
       
         
        on2=((angle(zo*rot)>=min(a1,a2))&(angle(rot*zo)<=max(a1,a2)))&~selected&~onedges;
%        plot(zo(on2)*rot,'gs')
%        on2=((angle(zo)>=angle(zp(ie+n/2)))&(angle(zo)<=angle(zp(ie+1+n/2))))&~selected&~onedges;

on2=find(on2);
[~,id]=sort(-angle(zo(on2)*rot));
on2=on2(id);
edges{ie,2}=on2;
onedges(on1)=1;
onedges(on2)=1;
 figure
 triplot(conn(:,1:3),xo,yo)
 hold on
 axis equal
 plot(xo(on1),yo(on1),'ro')
 plot(xo(on2),yo(on2),'mo')
 on12=[on1,on2];
 plot(xo(on12'),yo(on12'),'k-')
 %pause


    end
    
    medges=cell2mat(edges);
    medges=unique(medges,'rows');
    C1=sparse(medges(:,1),1:size(medges,1),1,prod(Nnodes),size(medges,1));
    C1=blkdiag(C1,C1);
    C2=sparse(medges(:,2),1:size(medges,1),1,prod(Nnodes),size(medges,1));
    C2=blkdiag(C2,C2);
    C=C1-C2;
    Ub=[(xo-xc),0*xo,(yo-yc)/sqrt(2),1/2*(xo-xc).*(xo-xc),1/2*(xo-xc).*(yo-yc),1/2*(yo-yc).*(xo-xc),1/2*(yo-yc).*(yo-yc),0*xo                ,0*xo                ,0*xo                ,0*xo;... % (x,y) x: vector of all x-coordinates
        0*yo,(yo-yc),(xo-xc)/sqrt(2),0*yo                ,0*yo                ,0*yo                ,0*yo                ,1/2*(xo-xc).*(xo-xc),1/2*(xo-xc).*(yo-yc),1/2*(yo-yc).*(xo-xc),1/2*(yo-yc).*(yo-yc)];   %  and y: vector of all y-coordinates
    Up=[Ub(medges(:,1),:)-Ub(medges(:,2),:);Ub(prod(Nnodes)+medges(:,1),:)-Ub(prod(Nnodes)+medges(:,2),:)];
    if 1
        Ce=blkdiag(Ce,Ce);
    C=[C,Ce];
    Ue=zeros(n,size(Up,2));
    for ie=1:n/2
        Ue(ie,:)=sum(Ub(edges{ie,1},:),1);
        Ue(ie+n/2,:)=sum(Ub(edges{ie,1}+prod(Nnodes),:),1);
    end
    Up=[Up;Ue];
    
    end
    
        Fo=zeros(2*prod(Nnodes),size(Up,2));
    Kc=[K,C;C',sparse(size(C,2),size(C,2))];
    Fc=[Fo;Up];

X=Kc\Fc;
U=X(1:2*prod(Nnodes),:);



else
    
    Ub=[(xo-xc),0*xo,(yo-yc)/sqrt(2),1/2*(xo-xc).*(xo-xc),1/2*(xo-xc).*(yo-yc),1/2*(yo-yc).*(xo-xc),1/2*(yo-yc).*(yo-yc),0*xo                ,0*xo                ,0*xo                ,0*xo;... % (x,y) x: vector of all x-coordinates
        0*yo,(yo-yc),(xo-xc)/sqrt(2),0*yo                ,0*yo                ,0*yo                ,0*yo                ,1/2*(xo-xc).*(xo-xc),1/2*(xo-xc).*(yo-yc),1/2*(yo-yc).*(xo-xc),1/2*(yo-yc).*(yo-yc)];   %  and y: vector of all y-coordinates
    bnd=~selected;
    indi=[find(bnd)];
    indi=[indi;indi+prod(Nnodes)];
    Up=Ub(indi,:);

    C=sparse(indi,1:length(indi),1,2*prod(Nnodes),length(indi));
    Fo=zeros(2*prod(Nnodes),size(Up,2));

    Kc=[K,C;C',sparse(size(C,2),size(C,2))];
    Fc=[Fo;Up];
 X=Kc\Fc;
U=X(1:2*prod(Nnodes),:);
 
end


for ii=1:size(U,2)
UpdateInternalState(nmod,U(:,ii),ii);
end
save(filres,'U','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
    postproVTK(filres,0,0);
    
%%
xo=xo(~selected);
yo=yo(~selected);
conn=delaunay(xo,yo);
Nnodes=[length(xo),1,1];
Nelems=[size(conn,1),1,1];
elt=3*ones(Nelems);
save('plein','Nnodes','Nelems','xo','yo','conn','elt','rint');
wp=GetWeigthDetJ('plein',[1,1],1,'Gauss_points');
So=sum(diag(wp))*pix2m*pix2m;  % w_g*det(J)
wdetJ=diag(wdetJ)*pix2m*pix2m;  % w_g*det(J)

Se=sum(wdetJ);  % w_g*det(J)
rho=Se/So;
A=zeros(size(U,2)-3);
H=zeros(3);
HH=zeros(6);
E=zeros(3);
iii=[1,2,4];
for ii=1:3
  load(sprintf('%s_%04d',filres,ii),'S');
  
  St=wdetJ'*S ; % S: stress
  H(:,ii)=St([1,2,4]).*[1,1,sqrt(2)];
  HH(:,iii(ii))=St.*[1,1,1,sqrt(2),sqrt(2),sqrt(2)];
end
H=H/So;
HH(1:2,3)=HH(3,1:2)';
HH(3,3)=sum(wdetJ)*model.material_parameters.young;
HH(5,5)=HH(4,4);
HH(6,6)=HH(4,4);
HH=HH/So;
bulk=eigs(HH,1)/3;

    %uxxx,uxyx,uxxy,uxyy,uyxx,uyyx,uyxy,uyyy
for ii=4:size(U,2)
  load(sprintf('%s_%04d',filres,ii),'S');
  
  St=S;  % S: stress
  
  A(1,ii-3)=wdetJ'*(St(:,1).*xg); % A 111 ... ~ tau 111
  A(2,ii-3)=wdetJ'*(1/2*(St(:,1).*yg+St(:,4).*xg)); % A 112 ... ~ tau 112
  A(3,ii-3)=wdetJ'*(1/2*(St(:,1).*yg+St(:,4).*xg)); % A 121 ... ~ tau 121
  A(4,ii-3)=wdetJ'*(St(:,4).*yg);  % A 122 ... ~ tau 122
  A(5,ii-3)=wdetJ'*(St(:,4).*xg);  % A 211 ... ~ tau 211
  A(6,ii-3)=wdetJ'*(1/2*(St(:,4).*yg+St(:,2).*xg));  % A 212 ... ~ tau 212
  A(7,ii-3)=wdetJ'*(1/2*(St(:,2).*xg+St(:,4).*yg)); % A 221 ... ~ tau 221
  A(8,ii-3)=wdetJ'*(St(:,2).*yg);  % A 211 ... ~ tau 222
  
end
%% xxx,xxy,xyx,xyy,yyx,yyy
A=A*pix2m/So;
inds=[1,2,4,5,6,8];
A(inds,inds);
%
A=transItoII(A);
% 
A(5:6,:)=[];
A(:,5:6)=[];
H=0.5*(H'+H);
HH=0.5*(HH'+HH);
A=0.5*(A'+A);
Ao=10^(round(log10(max(A(:)))));
A=round(A*10^nc/Ao)/10^nc*Ao
Ho=10^(round(log10(max(H(:)))));
H=round(H*10^nc/Ho)/10^nc*Ho
ct=sqrt(H(1,1)/(1e3*rho))
save(filres,'A','H','HH','bulk','-append')
end
