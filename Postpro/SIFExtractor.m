%function []=SIFExtractor(filres,parsif)
close all
clear all
filexp='-opti-none-rbt';
filres=['xfem',filexp];
load(filres);
parsif.box_size=[100,100,25];
parsif.nb_points_box=[10,10,5];
parsif.nb_points_front=10;
parsif.levelset_file=model.levelset_file;

nmod=1;
filfis=parsif.levelset_file;
load(filfis,'crack','front','s','zone');
sizezone=size(crack);
Nbox=parsif.nb_points_box-1;
hbox=mean(parsif.box_size./Nbox);
xi=(0:hbox:parsif.box_size(1));xi=xi-mean(xi);
yi=(0:hbox:parsif.box_size(2));yi=yi-mean(yi);
zi=(0:hbox:parsif.box_size(3));zi=zi-mean(zi);
[Yi,Xi,Zi]=meshgrid(yi,xi,zi);
figure
plot3(Xi(:),Yi(:),Zi(:),'x')

mesh_size=model.mesh_size;
roi=param.roi;
xo=xo+roi(1)-1;
yo=yo+roi(3)-1;
zo=zo+roi(5)-1;
M=ones(Nnodes);

xl=1:round(hbox):sizezone(1);
yl=1:round(hbox):sizezone(2);
zl=1:round(hbox):sizezone(3);

[Yl,Xl,Zl]=meshgrid(yl,xl,zl);
crackl=crack(xl,yl,zl);
frontl=front(xl,yl,zl);
sl=s(xl,yl,zl);


for iim=1:1%size(U,2)
    Ui=U(:,iim);
    Uglo=reshape(Ui((1:prod(Nnodes))),Nnodes);
    Vglo=reshape(Ui(prod(Nnodes)+(1:prod(Nnodes))),Nnodes);
    Wglo=reshape(Ui(2*prod(Nnodes)+(1:prod(Nnodes))),Nnodes);
  
    
    
on=find((crack<=0)&(crack>-1)&(front<=0)&(front>-1));
[xon,yon,zon]=ind2sub(size(crack),on);
xon=xon+zone(1)-1;
yon=yon+zone(3)-1;
zon=zon+zone(5)-1;


found=find((xon>min(xo))&(xon<max(xo))&(yon>min(yo))&(yon<max(yo))&(zon>min(zo))&(zon<max(zo)));
xon=xon(found);
yon=yon(found);
zon=zon(found);

son=s(on(found));
npts=parsif.nb_points_front;

hfront=(max(son)-min(son))/npts;
so=min(son):hfront:max(son);
ido=0*so;
no=zeros(length(so),3);
to=zeros(length(so),3);
oo=zeros(length(so),3);

for ip=1:1%npts
    found=find(abs(son-so(ip))==min(abs(son-so(ip))));
    so(ip)=son(found);
    ido(ip)=found;
    xi=xon(found)-zone(1)+1;
    xp=xon(found)-zone(1)+1+1;
    xm=xon(found)-zone(1)+1-1;
    yi=yon(found)-zone(3)+1;
    yp=yon(found)-zone(3)+1+1;
    ym=yon(found)-zone(3)+1-1;
    zi=zon(found)-zone(5)+1;
    zp=zon(found)-zone(5)+1+1;
    zm=zon(found)-zone(5)+1-1;
% n=grad(crack);t=grad(front);tn=t*n;
%     nx=0.5*(crack(xp,yi,zi)-crack(xm,yi,zi));
%     ny=0.5*(crack(xi,yp,zi)-crack(xi,ym,zi));
%     nz=0.5*(crack(xi,yi,zp)-crack(xi,yi,zm));
%     tx=0.5*(front(xp,yi,zi)-front(xm,yi,zi));
%     ty=0.5*(front(xi,yp,zi)-front(xi,ym,zi));
%     tz=0.5*(front(xi,yi,zp)-front(xi,yi,zm));
%     tnx=ty*nz-tz*ny;
%     tny=tz*nx-tx*nz;
%     tnz=tx*ny-ty*nx;

%t=grad(front);tn=grad(s);n=tn*t;
    tx=0.5*(front(xp,yi,zi)-front(xm,yi,zi));
    ty=0.5*(front(xi,yp,zi)-front(xi,ym,zi));
    tz=0.5*(front(xi,yi,zp)-front(xi,yi,zm));
    tnx=0.5*(s(xp,yi,zi)-s(xm,yi,zi));
    tny=0.5*(s(xi,yp,zi)-s(xi,ym,zi));
    tnz=0.5*(s(xi,yi,zp)-s(xi,yi,zm));
    nx=tny*tz-tnz*ty;
    ny=tnz*tx-tnx*tz;
    nz=tnx*ty-tny*tx;
    
    Xiglo=tx*Xi+nx*Yi+tnx*Zi;
    Yiglo=ty*Xi+ny*Yi+tny*Zi;
    Ziglo=tz*Xi+nz*Yi+tnz*Zi;
    
    
    Xiglo=Xiglo+xon(found);
    Yiglo=Yiglo+yon(found);
    Ziglo=Ziglo+zon(found);
    
    figure
    plot3(xo,yo,zo,'x')
    hold on
    plot3(xon,yon,zon,'ks')
    plot3(xon(found),yon(found),zon(found),'rs','MarkerSize',10)
    plot3(Xiglo(:),Yiglo(:),Ziglo(:),'r+');

    
    Xiglos=((Xiglo(:))-min(xo))/mesh_size(1)+1;
    Yiglos=((Yiglo(:))-min(yo))/mesh_size(2)+1;
    Ziglos=((Ziglo(:))-min(zo))/mesh_size(3)+1;
    
    figure
    plot3(((xo(:))-min(xo))/mesh_size(1)+1,((yo(:))-min(yo))/mesh_size(2)+1,((zo(:))-min(zo))/mesh_size(3)+1,'x')
    hold on
    plot3(Xiglos,Yiglos,Ziglos,'r+')
    
    
%     Uiglo=mexInterpLinear3D(Xiglos,Yiglos,Ziglos,Uglo);
%     Viglo=mexInterpLinear3D(Xiglos,Yiglos,Ziglos,Vglo);
%     Wiglo=mexInterpLinear3D(Xiglos,Yiglos,Ziglos,Wglo);
    Uiglo=interp3(Uglo,Yiglos,Xiglos,Ziglos,'*linear',0);
    Viglo=interp3(Vglo,Yiglos,Xiglos,Ziglos,'*linear',0);
    Wiglo=interp3(Wglo,Yiglos,Xiglos,Ziglos,'*linear',0);
    mask=interp3(M,Yiglos,Xiglos,Ziglos,'*linear',0);
    
    figure
    plot3(xo,yo,zo,'x')
    hold on
    plot3(xon,yon,zon,'ks')
    plot3(xon(found),yon(found),zon(found),'rs','MarkerSize',10)
    scatter3(Xiglo(:),Yiglo(:),Ziglo(:),0*mask+10,mask)
    
    Uiloc=Uiglo*tx+Viglo*ty+Wiglo*tz;
    Viloc=Uiglo*nx+Viglo*ny+Wiglo*nz;
    Wiloc=Uiglo*tnx+Viglo*tny+Wiglo*tnz;

    xn=max(2,min(sizezone(1)-1,Xiglo(:)-zone(1)+1));
    yn=max(2,min(sizezone(2)-1,Yiglo(:)-zone(3)+1));
    zn=max(2,min(sizezone(3)-1,Ziglo(:)-zone(5)+1));
    crackn=mexInterpLinear3D(xn,yn,zn,crack);
    frontn=mexInterpLinear3D(xn,yn,zn,front);
 crackn=reshape(crackn,Nbox);
 frontn=reshape(frontn,Nbox);
crackn=LSReinit(crackn,round(0.5*max(parsif.box_size)),hbox);
frontn=LSOrtho(frontn,crackn,round(0.5*max(parsif.box_size)),hbox);
frontn=LSReinit(frontn,round(0.5*max(parsif.box_size)),hbox);

[decn,A]=FissureFit(Uiloc,Viloc,Wiloc,crackn,frontn,mask,parsif);




    

end


    



end








%end