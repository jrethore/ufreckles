function postpro(filreso,comp)
if nargin<2, comp=false;end
disp(['Result file:',filreso]);
lim=[];
epsflag=false;
nmod=1;
load([filreso,'.mat'])
filexp='';
if ~exist('model1')
model1=model;
end
roi=param.roi;
if comp
    filreso=['comp-',filreso];
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
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
if isfield(model1,'crack_id')
    ncrack=model1.crack_id;
else
   ncrack=1; 
end
if strcmp(param.analysis,'mechanics')
    fac=pix2m;
else
    fac=1;
end
submean=0;

unstruc=0;
if isfield(model1,'mesh_file')
    unstruc=1;
end
if strcmp(param.analysis,'correlation')
    if exist([filreso,'-error.mat'],'file')
        fide=fopen([filreso,'-error.mat']);
    else
        fide=fopen(fullfile('TMP',[num2str(nmod),'_error_0.mat']));
    end
    dynamic=fread(fide,1);
end
%% mask

load(fullfile('TMP',[num2str(nmod),'_mask_0']),'mask');
mask=double(diag(mask)>0);
found=find(mask==0);
mask=full(mask);
mask(found)=NaN;
if unstruc
    xn=xo;
    yn=yo;
xo=[ceil(min(xo));floor(max(xo))];
yo=[ceil(min(yo));floor(max(yo))];
else
xo=reshape(xo,Nnodes)+0.5;
yo=reshape(yo,Nnodes)+0.5;
xo=xo(1:Nnodes(1),1)';
yo=yo(1,1:Nnodes(2));
xn=(xo);
xn(length(xn))=xo(length(xo))-1;
yn=(yo);
yn(length(yn))=yo(length(yo))-1;
end



Ut=U;
load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix','sizeim');
load(fullfile('TMP',[num2str(nmod),'_phiy_0']),'phiy','sizeim');
if unstruc
   hasval=(sum(phix,2)>0)|(sum(phiy,2)>0);
   found=find(~hasval);
   if numel(mask)==1
       mask=repmat(1,[prod(sizeim),1]);
   end
   mask(found)=NaN;
end
for iim=1:length(images)
    U=Ut(:,iim);
    if length(images)>1
close all
pause(0.1)
        filexp=sprintf('-%03d',images(iim));
    end
    filres=[filreso,filexp];
%% displacement map
umap=fac*reshape((phix*U).*mask,sizeim);
umap=umap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
if submean
    umap=umap-mean(umap(:));
end
%[limx]=FindMinMaxDisp(min(umap(:)),max(umap(:)));
if comp
    limx=[-10,10];
else
    limx=[];
end
PlotMap(umap,'Ux',limx,[filres,'-u'],epsflag);
vmap=fac*reshape((phiy*U).*mask,sizeim);
vmap=vmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
if submean
    vmap=vmap-mean(vmap(:));
end
%[limy]=FindMinMaxDisp(min(umap(:)),max(umap(:)));
if comp
    limy=[-1,1];
else
    limy=[];
end
PlotMap(vmap,'Uy',limy,[filres,'-v'],epsflag);
%VTKExportVectorMap(filres,'disp',xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1,umap,vmap);
[dudy,exx]=gradient(umap);
[eyy,dvdx]=gradient(vmap);
exy=0.5*(dudy+dvdx);
if comp
limyy=[-3.e-3,3e-3];
else
    limyy=[];
end
figure
cb=colormap(hot);
cb=cb(:,[3,2,1]);
colormap(cb);
if isempty(limyy)
axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
    'XTick',zeros(1,0),...
    'FontSize',20);
else
axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
    'XTick',zeros(1,0),...
    'FontSize',20,'CLim',limyy);
end
box('on');
hold 'all';
image(eyy','Parent',axes1,'CDataMapping','scaled');
        [c,h1] = contour(eyy',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
        set(h1,'EdgeColor','black','LineWidth',2);
axis off;
axis xy;
axis image;
colorbar('peer',axes1,'FontSize',20);
title('Eyy','FontSize',24);

        print ('-djpeg', fullfile('FIG',[filres,'-eyy']));
if strcmp(model1.basis,'fem')||strcmp(model1.basis,'nurbs')||strcmp(model1.basis,'beam')||strcmp(model1.basis,'nurbs-beam')||strcmp(model1.basis,'affine')
    if comp
        facx=15;%85
        %fac=0.2*max(sizeim)/max(abs(unmap(:)+i*vnmap(:)));
    else
        facx=0.2*max(sizeim)/max(abs(umap(:)+i*vmap(:)));
    end
    if comp
        limd=[0,12];%85
        %limd=[5,15];%60
    else
        limd=[min(abs(umap(:)+i*vmap(:))),max(abs(umap(:)+i*vmap(:)))];
    end
    DeformedMesh(xn,yn,umap,vmap,limd,facx,[filres,'-deformed-mesh'],epsflag,1);
end
clear umap vmap    
%% error map
if strcmp(param.analysis,'correlation')
    disc=fread(fide,prod(sizeim));
    disc=reshape(disc,sizeim);
    disc=100*abs(disc(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1))/dynamic;
    if comp
        limd=[0,10];
    else
    limd=[0,0.5*100*max(disc(:))/dynamic];
    end
   PlotMap(disc,'Discrepancy wrt dyn (%)',limd,[filres,'-error'],epsflag,0);
    disp(sprintf('Mean error wrt image dyn.: %6.2f %%',mean(disc(:))));
end

end
%% fleche...
U=Ut;
        beamtype=0;
if ~strcmp(model1.basis,'nurbs-beam')
    load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix','sizeim');
    for iim=1:size(U,2)
        utmp=reshape((phix*U(:,iim)).*mask,sizeim);
    utmp=utmp(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    umap(:,iim)=mean(utmp,1);
    clear phix
    dumap(:,iim)=gradient(utmp);
    ddumap(:,iim)=gradient(dumap(:,iim));
    end
else
    if isfield(model1,'beam_type')
        beamtype=strcmp(model1.beam_type,'timoshenko');
    end
    type_nurbs=~strcmp(model1.continuity,'c0');

    mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',1,1-1));
    pscale=1;
    p=model1.degree;
    load(mesh_file,'xo','yo','uo','vo','Nnodes','Nelems','Smesh');

    xo=(xo-0.5)*pscale+1;
yo=(yo-0.5)*pscale+1;
xo=reshape(xo,Nnodes);
yo=reshape(yo,Nnodes);
xo=xo(1:Nnodes(1),1)';
yo=yo(1,1:Nnodes(2));

    
    uo=(uo-1)*pscale+1;
    vo=(vo-1)*pscale+1;

    indp=zeros((p+1)*Smesh(2)*pscale,1);
    indn=zeros((p+1)*Smesh(2)*pscale,1);
    val=zeros((p+1)*Smesh(2)*pscale,1);
    dval=zeros((p+1)*Smesh(2)*pscale,1);

    ddval=zeros((p+1)*Smesh(2)*pscale,1);
    if type_nurbs
        Nny=length(yo)+p-1;
    else
        Nny=length(yo)*p-(p-1);
    end
    nel=0;
    toto=1;
    for iy=1:Nelems(2)
        yp=(yo(iy):(yo(iy+1)-1));

        if type_nurbs
            [N]=NURBSBasisFunc(iy+p,p,yp',vo-0.5,2);
        else
            toto=toto+p;
            [N]=NURBSBasisFunc(toto,p,yp',vo-0.5,2);
        end

        Sel=length(yp);
        for ip=1:(p+1)
            if type_nurbs
                indn(nel+(1:Sel))=iy+ip-1;
            else
                indn(nel+(1:Sel))=iy*p-p+ip;
            end
            indp(nel+(1:Sel))=yp;
            val(nel+(1:Sel))=N(:,ip,1);
            dval(nel+(1:Sel))=N(:,ip,2);
            ddval(nel+(1:Sel))=N(:,ip,3);
            nel=nel+Sel;
        end
    end
    phi1y=sparse(indp,indn,val,sizeim(2)*pscale,Nny);
    dphi1y=sparse(indp,indn,dval,sizeim(2)*pscale,Nny);
    ddphi1y=sparse(indp,indn,ddval,sizeim(2)*pscale,Nny);
    if beamtype
        umap=(phi1y*U(size(U,1)-2*size(dphi1y,2)+(1:size(dphi1y,2)),:))';
    else
        umap=(phi1y*U(size(U,1)-size(dphi1y,2)+(1:size(dphi1y,2)),:))';
    end
    dumap=(dphi1y*U(size(U,1)-size(dphi1y,2)+(1:size(dphi1y,2)),:))';
    ddumap=(ddphi1y*U(size(U,1)-size(dphi1y,2)+(1:size(dphi1y,2)),:))';
     [toto,dumap0]=gradient(umap);
     [toto,ddumap0]=gradient(dumap0);

 end
% figure1 = figure('XVisual',...
%     '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%     'PaperSize',[20.98 29.68]);
% axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
%     'FontName','Times');
% box('on');
% hold('all');
% plot(umap(2:size(umap,1)-1,:),'b-','DisplayName','Fleche  ','Parent',axes1,'LineWidth',2);
% title('Fleche','FontSize',20,'FontName','Times');
% xlabel('Position [pixel]','FontSize',20,'FontName','Times');
% ylabel('y [pixel]','FontSize',20,'FontName','Times');
% if comp
%     ylim([-1,10])
% end
% legend(axes1,'show','Location','NorthEast');
% print ('-djpeg', fullfile('FIG',[filres,'-fleche']));
% print ('-depsc', fullfile('FIG',[filres,'-fleche']));
% 
% figure1 = figure('XVisual',...
%     '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%     'PaperSize',[20.98 29.68]);
% axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
%     'FontName','Times');
% box('on');
% hold('all');
% plot(dumap(2:size(dumap,1)-1,:),'r-','DisplayName','Rotation','Parent',axes1,'LineWidth',2);
% title('Rotation','FontSize',20,'FontName','Times');
% xlabel('Position [pixel]','FontSize',20,'FontName','Times');
% ylabel('dy/dx []','FontSize',20,'FontName','Times');
% if comp
%     ylim([-0.03,0.03])
% end
% legend(axes1,'show','Location','NorthEast');
% print ('-djpeg', fullfile('FIG',[filres,'-rot']));
% print ('-depsc', fullfile('FIG',[filres,'-rot']));
% 
% figure1 = figure('XVisual',...
%     '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%     'PaperSize',[20.98 29.68]);
% axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
%     'FontName','Times');
% box('on');
% hold('all');
% plot(ddumap(3:size(dumap,1)-2,:),'k-','DisplayName','Courbure','Parent',axes1,'LineWidth',2);
% title('','FontSize',20,'FontName','Times');
% xlabel('Position [pixel]','FontSize',20,'FontName','Times');
% ylabel('d^2y/dx^2 [pixel^{-1}]','FontSize',20,'FontName','Times');
% if comp
%     ylim([-1.5e-4,0.5e-4])
% end
% legend(axes1,'show','Location','NorthEast');
% print ('-djpeg', fullfile('FIG',[filres,'-curv']));
% print ('-depsc', fullfile('FIG',[filres,'-curv']));
% 
% %resm=[(1:length(umap))',umap(:),dumap(:),ddumap(:)];
% 
% if beamtype
% %    resm=[(1:length(umap))',umap(:),dumap(:),ddumap(:),dumap0(:),ddumap0(:)];
%     figure1 = figure('XVisual',...
%         '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%         'PaperSize',[20.98 29.68]);
%     axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
%         'FontName','Times');
%     box('on');
%     hold('all');
%     plot(dumap0(2:size(dumap,1)-1,:),'r-','DisplayName','Euler-Bernoulli','Parent',axes1,'LineWidth',2);
%     plot(dumap(2:size(dumap,1)-1,:),'k-','DisplayName','Timoshenko','Parent',axes1,'LineWidth',2);
%     title('Rotation','FontSize',20,'FontName','Times');
%     xlabel('Position [pixel]','FontSize',20,'FontName','Times');
%     ylabel('dy/dx []','FontSize',20,'FontName','Times');
%     if comp
%         ylim([-0.03,0.03])
%     end
%     legend(axes1,'show','Location','NorthEast');
%     print ('-djpeg', fullfile('FIG',[filres,'-comp-rot']));
%     print ('-depsc', fullfile('FIG',[filres,'-comp-rot']));
% 
%     figure1 = figure('XVisual',...
%         '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%         'PaperSize',[20.98 29.68]);
%     axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
%         'FontName','Times');
%     box('on');
%     hold('all');
%     plot(ddumap(3:size(dumap,1)-2,:),'r-','DisplayName','Euler-Bernoulli','Parent',axes1,'LineWidth',2);
%     plot(ddumap0(3:size(dumap,1)-2,:),'k-','DisplayName','Timoshenko','Parent',axes1,'LineWidth',2);
%     title('','FontSize',20,'FontName','Times');
%     xlabel('Position [pixel]','FontSize',20,'FontName','Times');
%     ylabel('d^2y/dx^2 [pixel^{-1}]','FontSize',20,'FontName','Times');
%     if comp
%         ylim([-1.5e-4,0.5e-4])
%     end
%     legend(axes1,'show','Location','NorthEast');
%     print ('-djpeg', fullfile('FIG',[filres,'-comp-curv']));
%     print ('-depsc', fullfile('FIG',[filres,'-comp-curv']));
% 
% 
% end
% save([filres,'-curv'],'umap','dumap','ddumap');
yi=1:size(umap,2);
figure1 = figure('XVisual',...
    '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98 29.68]);
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
    'FontName','Times');
box('on');
hold('all');
for iim=1:length(images)
    plot(pix2m*(yi)',pix2m*umap(iim,:)','DisplayName',num2str(images(iim)),'Parent',axes1,'LineWidth',2);
end
title('Contour','FontSize',20,'FontName','Times');
if pix2m==1
    xlabel('Position [pixel]','FontSize',20,'FontName','Times');
ylabel('y [pixel]','FontSize',20,'FontName','Times');
else
    xlabel('Position [m]','FontSize',20,'FontName','Times');
ylabel('y [m]','FontSize',20,'FontName','Times');
end
if comp
    ylim([-1,10])
end
legend(axes1,'show','Location','NorthEast');
print ('-djpeg', fullfile('FIG',[filreso,'-contour']));
print ('-depsc', fullfile('FIG',[filreso,'-contour']));

figure1 = figure('XVisual',...
    '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98 29.68]);
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
    'FontName','Times');
box('on');
hold('all');
for iim=1:length(images)
    plot(pix2m*(yi)',dumap(iim,:)','DisplayName',num2str(images(iim)),'Parent',axes1,'LineWidth',2);
end
title('Rotation','FontSize',20,'FontName','Times');
if pix2m==1
xlabel('Position [pixel]','FontSize',20,'FontName','Times');
ylabel('dy/dx []','FontSize',20,'FontName','Times');
else
    xlabel('Position [m]','FontSize',20,'FontName','Times');
ylabel('dy/dx []','FontSize',20,'FontName','Times');
end
if comp
    ylim([-0.03,0.03])
end
legend(axes1,'show','Location','NorthEast');
print ('-djpeg', fullfile('FIG',[filreso,'-rot']));
print ('-depsc', fullfile('FIG',[filreso,'-rot']));

figure1 = figure('XVisual',...
    '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'PaperSize',[20.98 29.68]);
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
    'FontName','Times');
box('on');
hold('all');
for iim=1:length(images)
    plot(pix2m*(yi)',ddumap(iim,:)'/pix2m,'DisplayName',num2str(images(iim)),'Parent',axes1,'LineWidth',2);
end
title('Curvature','FontSize',20,'FontName','Times');
if pix2m==1
xlabel('Position [pixel]','FontSize',20,'FontName','Times');
ylabel('d^2y/dx^2 [pixel^{-1}]','FontSize',20,'FontName','Times');
else
    xlabel('Position [m]','FontSize',20,'FontName','Times');
ylabel('d^2y/dx^2 [m^{-1}]','FontSize',20,'FontName','Times');
end
if comp
    ylim([-1.5e-4,0.5e-4])
end
legend(axes1,'show','Location','NorthEast');
print ('-djpeg', fullfile('FIG',[filreso,'-curv']));
print ('-depsc', fullfile('FIG',[filreso,'-curv']));
save([filreso,'-data.mat'],'umap','dumap','ddumap','yo');
%% ROI
if strcmp(param.analysis,'correlation')
    load(fullfile('TMP','sample0'),'im0')
    PlotZOI(roi,(xo-1)*psample+roi(1),(yo-1)*psample+roi(3),im0,'ZOI',[filres,'-zoi'],epsflag);
    clear im0

end
if strcmp(param.analysis,'correlation')
    fclose(fide);
end
end


