%%
clear all
close all
siz=128;

[Yi,Xi]=meshgrid(1:(4*siz),1:(2*siz));

%toto=2*(Xi<(siz+8*sin(2*pi*Yi/siz)))-1;
toto=2*(Xi<(siz+(2*siz-Yi)*0.125))-1;
titi=2*(Yi>2*siz)-1;

%%
tic;
tata=toto;
dt=0.5;
maxiter=floor(2*siz/dt);
for iter=1:maxiter
    tata=mexLSReinit0(tata,1,dt);
end

toc
tic;
tutu=mexLSReinit(toto,1,2*siz);
toc
% tic;
% tete=nvmexLSReinit(toto,1,2*siz);
% toc

%%
figure
imagesc(toto);colorbar;title('toto init');axis image

dt=0.5;
maxiter=floor(2*siz/dt);
hold on

   [c,h1] = contour(toto,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
for iter=1:maxiter
    tata=mexLSReinit0(toto,1,dt);
    toto=tata;
end
imagesc(toto);colorbar;title('toto reinit');axis image
   [c,h1] = contour(toto,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   pause(.250)



figure
imagesc(titi);colorbar;title('titi init');axis image

maxiter=floor(4*siz/dt);
for iter=1:maxiter
    titi=mexLSReinit0(titi,1,dt);
end
figure
imagesc(titi);colorbar;title('titi reinit');axis image
hold on
   [c,h1] = contour(titi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
%%
maxiter=2*floor(siz/dt);
%maxiter=100;
tata=titi;
figure
imagesc(tata);colorbar;title('titi reinit reortho');axis image
hold on
   [c,h2] = contour(toto,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
   [c,h1] = contour(tata,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
for iter=1:maxiter
    tata=mexLSReOrtho(tata,toto,1,dt);
end
    delete(h1)
imagesc(tata);colorbar;
   [c,h2] = contour(toto,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
[c,h1] = contour(tata,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
    pause(.1250)
%%
dt=0.8
vini=Xi/50;
maxiter=.1*floor(siz/dt);
figure
imagesc(vini);colorbar;title('vini');axis image
hold on
   [c,h1] = contour(tata,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
    pause(.1250)

for iter=1:maxiter
    vini=mexLSExtend(vini,tata,1,dt);
    
    imagesc(vini);colorbar;title('vini');axis image
   [c,h1] = contour(tata,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
    pause(.1250)

end
%%

tutu=tata;
maxiter=1*floor(siz/dt);
figure
h0=imagesc(tutu);colorbar;title('tutu');axis image
hold on
   [c,h1] = contour(tutu,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
    pause(.1250)
ndt=ceil(max(vini(:))/0.25);
dt=1/ndt;
tt=0;
    while tt<=ndt
tutu=mexLSPropagate(tutu,vini,1,dt);
delete(h0);delete(h1);
    h0=imagesc(tutu);colorbar;title('tutu');axis image
   [c,h1] = contour(tutu,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(tata,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
    pause(.1250)
tt=tt+1;
    end
%%



clear all
close all
siz=24;

[Yi,Xi]=meshgrid(1:(4*siz),1:(2*siz));

%toto=2*(Xi<(siz+8*sin(2*pi*Yi/siz)))-1;
phi=2*(Xi<siz)-1;%crack
psi=2*(Yi>2*siz)-1;%front

figure
imagesc(phi);colorbar;title('phi init');axis image

dt=0.5;
maxiter=floor(2*siz/dt);
hold on

   [c,h1] = contour(phi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
for iter=1:maxiter
    phi=mexLSReinit0(phi,1,dt);
imagesc(phi);colorbar;title('phi reinit');axis image
   [c,h1] = contour(phi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   pause(.1250)
end



figure
imagesc(psi);colorbar;title('psi init');axis image

maxiter=floor(4*siz/dt);
for iter=1:maxiter
    psi=mexLSReinit0(psi,1,dt);
end
figure
imagesc(psi);colorbar;title('psi reinit');axis image
hold on
   [c,h1] = contour(psi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);

%%

vphi=0.25*siz+0*Xi;
vpsi=0.5*siz+0*Xi;

maxiter=floor(siz/dt);
for iter=1:maxiter
    vphi=mexLSExtend(vphi,psi,1,dt);
end
for iter=1:maxiter
    vpsi=mexLSExtend(vpsi,psi,1,dt);
end

for iter=1:maxiter
    vphi=mexLSExtend(vphi,phi,1,dt);
end
for iter=1:maxiter
    vpsi=mexLSExtend(vpsi,phi,1,dt);
end


vphi=vphi.*psi./(vpsi*1);
vphi=vphi.*double(psi>=0);
%%
ndt=ceil(max(abs(vphi(:)))/0.25);
dtp=1/ndt;
figure
h0=imagesc(phi);colorbar;title('phi');axis image
hold on
   [c,h1] = contour(phi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
    pause(.1250)
tt=0;
    while tt<=ndt
phi=mexLSPropagate(phi,vphi,1,dtp);
delete(h0);delete(h1);
    h0=imagesc(phi);colorbar;title('phi');axis image
   [c,h1] = contour(phi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(psi,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
    pause(.1250)
tt=tt+1;
    end
    
%%    
    
ndt=ceil(max(abs(vpsi(:)))/0.25);
dtp=1/ndt;
figure
h0=imagesc(psi);colorbar;title('psi');axis image
hold on
   [c,h1] = contour(psi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
    pause(.1250)
tt=0;
    while tt<=ndt
psi=mexLSPropagate(psi,vpsi,1,dtp);
delete(h0);delete(h1);
    h0=imagesc(psi);colorbar;title('psi');axis image
   [c,h1] = contour(psi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(phi,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
    pause(.1250)
tt=tt+1;
    end

%%    
    figure
imagesc(phi);colorbar;title('phi fini');axis image

dt=0.5;
maxiter=floor(0.25*siz/dt);
hold on

   [c,h1] = contour(phi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(psi,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
for iter=1:maxiter
    phi=mexLSReinit0(phi,1,dt);
imagesc(phi);colorbar;title('phi fini');axis image
   [c,h1] = contour(phi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(psi,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
   pause(.250)
end
%%
    figure
imagesc(psi);colorbar;title('psi fini');axis image

dt=0.5;
maxiter=floor(siz/dt);
hold on

   [c,h1] = contour(psi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(phi,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
for iter=1:maxiter
    psi=mexLSReOrtho(psi,phi,1,dt);
imagesc(phi);colorbar;title('psi fini');axis image
   [c,h1] = contour(psi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(phi,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
   pause(.1250)
end
for iter=1:maxiter
    psi=mexLSReinit0(psi,1,dt);
imagesc(phi);colorbar;title('psi fini');axis image
   [c,h1] = contour(psi,[0. 0.]);
   set(h1,'EdgeColor','white','LineWidth',2);
   [c,h2] = contour(phi,[0. 0.]);
   set(h2,'EdgeColor','black','LineWidth',2);
   pause(.1250)
end

    
    