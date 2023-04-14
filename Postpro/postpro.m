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
if isfield(model1,'mesh_file')||isfield(model1,'adapt_mesh')
    unstruc=1;
end
if isfield(model1,'element_type')
    if model1.element_type==3
        unstruc=1;
    end
end
%unstruc=1;
if strcmp(param.analysis,'correlation')
    if exist([filreso,'-error.mat'],'file')
        fide=fopen([filreso,'-error.mat']);
    else
        fide=fopen(fullfile('TMP',[num2str(nmod),'_error_0.mat']));
    end
    erroronelt=fread(fide,1);
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
   hasval=(sum(abs(phix),2)>0)|(sum(abs(phiy),2)>0);
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
if strcmp(model1.basis,'fem')||strcmp(model1.basis,'nurbs')||strcmp(model1.basis,'beam')||strcmp(model1.basis,'nurbs-beam')||strcmp(model1.basis,'affine')
    if comp
        facx=15;%85
        %fac=0.2*max(sizeim)/max(abs(unmap(:)+i*vnmap(:)));
    else
        facx=0.2*max(sizeim)/max(abs(umap(:)-mean(umap(:))+i*(vmap(:)-mean(vmap(:)))));
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
if strcmp(param.analysis,'correlation')&&~erroronelt
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
%% displacement jump
if    (strcmp(model1.basis,'fem')&&isfield(model1,'enrichment'))||strcmp(model1.basis,'KM+')||strcmp(model1.basis,'KM')||strcmp(model1.basis,'CZ')
    load(fullfile('TMP',[num2str(ncrack),'_levelsets']),'onfaces','front');
    load(fullfile('TMP',[num2str(ncrack),'_levelsets_cylco']),'dist');
    if isfield(model1,'mask_radius')
        rmin=model1.mask_radius;
    else
        rmin=0;
    end
    mask1=double((dist(:)>rmin));
    found1=find(~mask1(:));
    mask1(found1)=NaN;
    mask1(found)=NaN;

    onfaces=onfaces(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    onfaces=find(onfaces(:));
    front=front(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    xfaces=fac*front(onfaces);
    [xfaces,ind]=sort(xfaces,'descend');
    onfaces=onfaces(ind);
    clear front

    load(fullfile('TMP',[num2str(nmod),'_dphit_0']),'dphit');
    uxmap=fac*reshape((dphit*U).*mask1,sizeim);
    uxmap=uxmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    uxmap=uxmap(onfaces);
    load(fullfile('TMP',[num2str(nmod),'_dphin_0']),'dphin');
    vxmap=fac*reshape((dphin*U).*mask1,sizeim);
    vxmap=vxmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    vxmap=vxmap(onfaces);

    nxmap=abs(uxmap+i*vxmap);
    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(-xfaces(:),nxmap(:),'r-','DisplayName','Norm    ','Parent',axes1,'LineWidth',3);
    plot(-xfaces(:),vxmap(:),'b-','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    plot(-xfaces(:),uxmap(:),'g-','DisplayName','Tangent ','Parent',axes1,'LineWidth',2);
    title('Displacement jump','FontSize',20,'FontName','Times');
    xlabel('Position','FontSize',20,'FontName','Times');
    ylabel('[U]','FontSize',20,'FontName','Times');
    if comp
%        ylim([-1,6])
    end
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-jump']));
    print ('-depsc', fullfile('FIG',[filres,'-jump']));
    resm=full([xfaces(:),uxmap(:),vxmap(:),nxmap(:)]);
    save([filres,'-jump.dat'],'resm','-ASCII');

    if exist('W')
    if isfield(param,'zeta')
        if param.zeta==0
            zeta=0;
        else
        zeta=1/param.zeta;
        end
    else
        zeta=0;
    end
    load(fullfile('TMP',sprintf('%d_iphi_0',nmod)),'ztip','zn','iphin','iphino');
    load(fullfile('TMP',[num2str(nmod),'_dpsit_0']),'dpsitcz','dpsitczo');
    ut=fac*dpsitcz*W-dpsitczo*S*zeta;
    ut=[ut(:);fac*iphin*W(length(W)/2+(1:length(W)/2))-zeta*iphino*S(length(S)/2+(1:length(S)/2))];
    load(fullfile('TMP',[num2str(nmod),'_dpsin_0']),'dpsincz','dpsinczo');
    un=fac*dpsincz*W-dpsinczo*S*zeta;
    un=[un(:);fac*iphin*W((1:length(W)/2))-zeta*iphino*S((1:length(S)/2))];
    xfaces=ztip*fac;
    xfaces=[xfaces(:);fac*zn(:)];
    [xfaces,ind]=sort(xfaces,'ascend');
    ut=ut(ind);un=un(ind);
    ns=abs(ut+i*un);
    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(-xfaces,ns(:),'r-','DisplayName','Norm    ','Parent',axes1,'LineWidth',3);
    plot(-xfaces,un(:),'b-','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    plot(-xfaces,ut(:),'g-','DisplayName','Tangent ','Parent',axes1,'LineWidth',2);
    title('Displacement jump (interface)','FontSize',20,'FontName','Times');
    xlabel('Position','FontSize',20,'FontName','Times');
    ylabel('[U]','FontSize',20,'FontName','Times');
    %    ylim([-0.25,0.5])
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-jump-cz']));
    print ('-depsc', fullfile('FIG',[filres,'-jump-cz']));
    resm=full([xfaces(:),ut(:),un(:),ns(:)]);
    save([filres,'-jump-cz.dat'],'resm','-ASCII');


        
        
        
    end
    
    
    
end
%% cohesive tractions
if    (strcmp(model1.basis,'fem')&&isfield(model1,'enrichment'))&&exist('S')

%     load(sprintf('TMP/%_interface',nmod),'fn','ft');
%     uxmap=reshape(ft,sizeim);
%     uxmap=uxmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
%     uxmap=uxmap(onfaces);
%     vxmap=reshape(fn,sizeim);
%     vxmap=vxmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
%     vxmap=vxmap(onfaces);

    load(fullfile('TMP',[num2str(ncrack),'_levelsets']),'onfaces','front');
    load(fullfile('TMP',[num2str(ncrack),'_levelsets_cylco']),'dist');
    if isfield(model1,'mask_radius')
        rmin=model1.mask_radius;
    else
        rmin=0;
    end
    if isfield(param,'zeta')
        zeta=param.zeta;
    else
        zeta=0;
    end
    if exist('W') zeta=0;end
    mask1=double((front(:)>0));
    found1=find(~mask1(:));
    mask1(found1)=NaN;
onfaces=onfaces.*double((front>0));
    onfaces=onfaces(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    onfaces=find(onfaces(:));
    front=front(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    xfaces=fac*front(onfaces);
    [xfaces,ind]=sort(xfaces,'descend');
    onfaces=onfaces(ind);
    clear front

    load(fullfile('TMP',[num2str(nmod),'_dphit_0']),'dphit');
    uxmap=fac*reshape((dphit*U).*mask1,sizeim);
    uxmap=uxmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    uxmap=uxmap(onfaces);
    load(fullfile('TMP',[num2str(nmod),'_dphin_0']),'dphin');
    vxmap=fac*reshape((dphin*U).*mask1,sizeim);
    vxmap=vxmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    vxmap=vxmap(onfaces);
    
    if size(S,1)==size(U,1)
        
     qt=reshape((dphit*S).*mask1,sizeim);
    qn=reshape((dphin*S).*mask1,sizeim);
       
    else
    load(fullfile('TMP',[num2str(nmod),'_dpsit_0']),'dpsito');
    qt=reshape((dpsito*S).*mask1,sizeim);
    load(fullfile('TMP',[num2str(nmod),'_dpsin_0']),'dpsino');
    qn=reshape((dpsino*S).*mask1,sizeim);
    end
    qt=qt(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    qt=qt(onfaces);
    qn=qn(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
    qn=qn(onfaces);
    
b=2*(1+0.3);
    sn=qn-zeta*vxmap;
    st=qt-zeta*uxmap/b;
ns=abs(st+i*sn);

    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(-xfaces(:),ns(:),'r-','DisplayName','Norm    ','Parent',axes1,'LineWidth',3);
    plot(-xfaces(:),sn(:),'b-','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    plot(-xfaces(:),st(:),'g-','DisplayName','Tangent ','Parent',axes1,'LineWidth',2);
    title('Interface tractions','FontSize',20,'FontName','Times');
    xlabel('Position','FontSize',20,'FontName','Times');
    ylabel('T','FontSize',20,'FontName','Times');
    %    ylim([-0.25,0.5])
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-traction']));
    print ('-depsc', fullfile('FIG',[filres,'-traction']));
    resm=full([xfaces(:),st(:),sn(:),ns(:)]);
    save([filres,'-traction.dat'],'resm','-ASCII');
if isfield(model1,'tip_step')||isfield(model1,'interface_type')
    
                load(fullfile('TMP',sprintf('%d_iphi_0',nmod)),'ztip','zn','dz','iphin','iphino');
        load(fullfile('TMP',[num2str(nmod),'_dphit_0']),'dphitcz');
    load(fullfile('TMP',[num2str(nmod),'_dphin_0']),'dphincz');
    
    
    load(fullfile('TMP',[num2str(nmod),'_dpsit_0']),'dpsitczo');
    qt=dpsitczo*S;
    load(fullfile('TMP',[num2str(nmod),'_dpsin_0']),'dpsinczo');
    qn=dpsinczo*S;
    qt=[qt(:);iphino*S(length(S)/2+(1:length(S)/2))];
    qn=[qn(:);iphino*S((1:length(S)/2))];

%     dz=diag(sparse(dz));
%     Mint=dpsinczo'*dz*dpsinczo+dpsitczo'*dz*dpsitczo;
% %    full(Mint)
%    Fint=fac*dpsinczo'*dz*dphincz*U+fac*dpsitczo'*dz*dphitcz*U;
%    Uint=Mint\Fint;
%    
%    uxmap=dpsitczo*Uint;
%    uxmap=[uxmap(:);iphino*Uint(length(Uint)/2+(1:length(Uint)/2))];
%    vxmap=dpsinczo*Uint;
%    vxmap=[vxmap(:);iphino*Uint((1:length(Uint)/2))];
   uxmap=fac*dphitcz*U;
    uxmap=[uxmap(:);interp1(ztip,uxmap,zn(:))];
   vxmap=fac*dphincz*U;
    vxmap=[vxmap(:);interp1(ztip,vxmap,zn(:))];
    xfaces=ztip*fac;    
        xfaces=[xfaces(:);fac*zn(:)];
    [xfaces,ind]=sort(xfaces,'ascend');
    qt=qt(ind);
    qn=qn(ind);
    uxmap=uxmap(ind);vxmap=vxmap(ind);
    
%     
b=2*(1+0.3);
    sn=qn-zeta*vxmap;
    st=qt-zeta*uxmap/b;



ns=abs(st+i*sn);

    
    
    
    
    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(-xfaces,ns(:),'r-','DisplayName','Norm    ','Parent',axes1,'LineWidth',3);
    plot(-xfaces,sn(:),'b-','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    plot(-xfaces,st(:),'g-','DisplayName','Tangent ','Parent',axes1,'LineWidth',2);
    title('Interface tractions (interface)','FontSize',20,'FontName','Times');
    xlabel('Position','FontSize',20,'FontName','Times');
    ylabel('T','FontSize',20,'FontName','Times');
    %    ylim([-0.25,0.5])
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-traction-cz']));
    print ('-depsc', fullfile('FIG',[filres,'-traction-cz']));
    resm=full([xfaces(:),st(:),sn(:),ns(:)]);
    save([filres,'-traction-cz.dat'],'resm','-ASCII');

    
end
elseif strcmp(model1.basis,'KM+')&&(param.cz_length>0)

    load(fullfile('TMP',[num2str(ncrack),'_levelsets']),'onfaces');
    lcz=param.cz_length;
    load(fullfile('TMP',sprintf('%d_iphi_%d',nmod,0)),'ztip');


    clear front
    load(fullfile('TMP',[num2str(nmod),'_dpsit_0']),'dpsitcz');
    txmap2=dpsitcz*U;
    load(fullfile('TMP',[num2str(nmod),'_dpsin_0']),'dpsincz');
    nxmap2=dpsincz*U;
    load(fullfile('TMP',[num2str(nmod),'_dphit_0']),'dphitcz');
    uxmap2=pix2m*dphitcz*U;
    load(fullfile('TMP',[num2str(nmod),'_dphin_0']),'dphincz');
    vxmap2=pix2m*dphincz*U;

    nnxmap2=abs(txmap2+i*nxmap2);
    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(ztip/max(ztip),nnxmap2(:),'r--','DisplayName','Norm    ','Parent',axes1,'LineWidth',2);
    plot(ztip/max(ztip),nxmap2(:),'b--','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    plot(ztip/max(ztip),txmap2(:),'g--','DisplayName','Tangent ','Parent',axes1,'LineWidth',2);
    title('Interface tractions','FontSize',20,'FontName','Times');
    xlabel('Position','FontSize',20,'FontName','Times');
    ylabel('T','FontSize',20,'FontName','Times');
    %    ylim([-0.25,0.5])
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-traction']));
    print ('-depsc', fullfile('FIG',[filres,'-traction']));


    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(vxmap2(:),nxmap2(:),'b-','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    title('Cohesive law','FontSize',20,'FontName','Times');
    xlabel('[U]','FontSize',20,'FontName','Times');
    ylabel('T','FontSize',20,'FontName','Times');
    %    ylim([-0.25,0.5])
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-czm']));
    print ('-depsc', fullfile('FIG',[filres,'-czm']));
    
    resm=[ztip',uxmap2(:),vxmap2(:),txmap2(:),nxmap2(:)];
         save([filres,'-jump-trac-cz.dat'],'resm','-ASCII');

elseif strcmp(model1.basis,'CZ')
    load(fullfile('TMP',[num2str(nmod),'_dphit_0']),'dphitcz');
    uxmap2=pix2m*dphitcz*U;
    load(fullfile('TMP',[num2str(nmod),'_dphin_0']),'dphincz');
    vxmap2=pix2m*dphincz*U;


    load(fullfile('TMP',[num2str(nmod),'_iphi_0']),'iphi','ztip');

    P=iphi*U(4+(1:size(iphi,2)));
    Q=iphi*U(4+size(iphi,2)+(1:size(iphi,2)));
        figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(ztip/max(ztip),abs(P(:)+i*Q(:)),'r--','DisplayName','Norm    ','Parent',axes1,'LineWidth',2);
    plot(ztip/max(ztip),P(:),'b--','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    plot(ztip/max(ztip),Q(:),'g--','DisplayName','Tangent ','Parent',axes1,'LineWidth',2);
    title('Interface tractions','FontSize',20,'FontName','Times');
    xlabel('Position','FontSize',20,'FontName','Times');
    ylabel('T','FontSize',20,'FontName','Times');
    %    ylim([-0.25,0.5])
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-traction']));
    print ('-depsc', fullfile('FIG',[filres,'-traction']));


    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(vxmap2(:),P(:),'b-','DisplayName','Normal  ','Parent',axes1,'LineWidth',2);
    title('Cohesive law','FontSize',20,'FontName','Times');
    xlabel('[U]','FontSize',20,'FontName','Times');
    ylabel('T','FontSize',20,'FontName','Times');
    %    ylim([-0.25,0.5])
    legend(axes1,'show','Location','NorthWest');
    print ('-djpeg', fullfile('FIG',[filres,'-czm']));
    print ('-depsc', fullfile('FIG',[filres,'-czm']));
end
%% SIF
if strcmp(model1.basis,'KM')

    ind=model1.km_indices;
    modes=model1.modes;
    kappa=model1.kolosov;
    mu=model1.mu;
    pix2m=param.pixel_size;
    scal=2*mu*sqrt(2*pi);
    scalamp=scal*(pix2m.^(1-(ind)*.5));
    fid=fopen([filres,'-ks.dat'],'w');

    found=find(ind==1);
    if found
        for m=1:length(modes)
            if modes(m)==1
                K1=scalamp(found)*U(found+(modes(m)-1)*length(ind));
                disp(sprintf('Mode I SIF : %10.3e Pa.sqrt(m)',K1));
                fprintf(fid,'%d %d %12.5e\n',1,1,K1);
            elseif modes(m)==2
                K2=scalamp(found)*U(found+(modes(m)-1)*length(ind));
                disp(sprintf('Mode II SIF : %10.3e Pa.sqrt(m)',K2));
                fprintf(fid,'%d %d %12.5e\n',2,1,K2);

            end
        end
    end
    found=find(ind==-1);
    if found
        for m=1:length(modes)
            if modes(m)==1
                SK1=scalamp(found)*U(found+(modes(m)-1)*length(ind));
                disp(sprintf('Crack tip shift : %10.3f pixel %10.3e m',-2*SK1/K1/pix2m,-2*SK1/K1));
                fprintf(fid,'%d %d %12.5e\n',1,-1,SK1);

            end
        end
    end

    found=find(ind==-3);
    if found
        for m=1:length(modes)
            if modes(m)==1
                SK3=scalamp(found)*U(found+(modes(m)-1)*length(ind));
                disp(sprintf('Process zone width indicator : %10.3f pixel %10.3e m',abs(sqrt(-8*SK3/K1))/pix2m,abs(sqrt(-8*SK3/K1))));
                fprintf(fid,'%d %d %12.5e\n',1,-3,SK3);
            end
        end
    end
    fclose(fid);
elseif strcmp(model1.basis,'KM+')
    lmin=param.cz_length;
    modes=model1.modes;
    kappa=model1.kolosov;
    mu=model1.mu;
    pix2m=param.pixel_size;
    scal=2*mu*sqrt(2*pi);
    scalamp=scal*(pix2m.^(1-(1)*.5));
    fid=fopen([filres,'-ks.dat'],'w');
    load(fullfile('TMP',[num2str(nmod),'_iphi_0']),'iphi','ztip');

    K1=scalamp*iphi*U(4+(1:size(iphi,2)));
    K2=scalamp*iphi*U(4+size(iphi,2)+(1:size(iphi,2)));
    fprintf(fid,'%d %d %12.5e\n',ztip',K1,K2);
    fclose(fid);
    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(ztip/max(ztip),K1(:)*1.e-6,'b-','DisplayName','K_{I}','Parent',axes1,'LineWidth',2);
    plot(ztip/max(ztip),K2(:)*1.e-6,'r-','DisplayName','K_{II}','Parent',axes1,'LineWidth',2);
    title('Stress intensity factors distribution','FontSize',20,'FontName','Times');
    xlabel('Position [pixel]','FontSize',20,'FontName','Times');
    ylabel('K_{I,II} MPa.m^{1/2}','FontSize',20,'FontName','Times');
    legend(axes1,'show','Location','NorthEast');
    print ('-djpeg', fullfile('FIG',[filres,'-fics']));
    print ('-depsc', fullfile('FIG',[filres,'-fics']));
elseif  strcmp(model1.basis,'fem')&&isfield(model1,'enrichment')

%         load(sprintf('TMP/%d_enrichment_%d',nmod,ncrack));
    if strcmp(model1.enrichment,'near_tip')
        load(fullfile('TMP',[num2str(nmod),'_mask_0']),'unmasked_nodes');

        if isempty(unmasked_nodes)
            Nddl=2*prod(Nnodes);
        else
            Nddl=2*length(unmasked_nodes);
        end

        ind=model1.km_indices;
        modes=model1.modes;
        kappa=model1.kolosov;
        mu=model1.mu;
        pix2m=param.pixel_size;
        scal=2*mu*sqrt(2*pi);
        scalamp=scal*(pix2m.^(1-(ind)*.5));
        found=find(ind==1);
        UK=repmat(0,length(U),1);
        found=find(ind==1);
        if ~isempty(found)
            if strcmp(model1.enrichment_type,'zone')
                K1n=scalamp(found)*U(Nddl+found);
                K2n=scalamp(found)*U(Nddl+found+length(ind));
            elseif strcmp(model1.enrichment_type,'node')
                K1n=scalamp(found)*U(Nddl+(found:(length(modes)*length(ind)):(length(modes)*length(tip_nodes)*length(ind))));
                K2n=scalamp(found)*U(Nddl+((found+length(ind)):(length(modes)*length(ind)):(length(modes)*length(tip_nodes)*length(ind))));
            end
            load(fullfile('TMP',[num2str(ncrack),'_levelsets_cylco']),'dist');
            dist=dist(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
            tip=find(dist(:)==min(dist(:)));
            clear dist
            load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix','sizeim');
            K1=0*U;
            K1(tip_nodes)=K1n(:);
            kmap=reshape((phix*K1),sizeim);
            kmap=kmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
            if strcmp(model1.enrichment_type,'zone')
                K1o=K1n;
            elseif strcmp(model1.enrichment_type,'node')
                K1o=mean(kmap(tip));
            end
            PlotMap(kmap,'K1 map Pa.sqrt(m)',[],[filres,'-k1-map'],epsflag);
            K2=0*U;
            K2(tip_nodes)=K2n(:);
            kmap=reshape((phix*K2),sizeim);
            kmap=kmap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
            if strcmp(model1.enrichment_type,'zone')
                K2o=K2n;
            elseif strcmp(model1.enrichment_type,'node')
                K2o=mean(kmap(tip));
            end
            PlotMap(kmap,'K2 map Pa.sqrt(m)',[],[filres,'-k2-map'],epsflag);
            clear phix


            [Yo Xo]=meshgrid(yo,xo);
            fid=fopen([filres,'-k1.dat'],'w');
            fprintf(fid,'%12.5e %12.5e %12.5e \n',0,0,K1o);
            fprintf(fid,'%12.5e %12.5e %12.5e \n',Xo(tip_nodes),Yo(tip_nodes),K1n(:));
            fclose(fid);
            fid=fopen([filres,'-k2.dat'],'w');
            fprintf(fid,'%12.5e %12.5e %12.5e \n',0,0,K2o);
            fprintf(fid,'%12.5e %12.5e %12.5e \n',Xo(tip_nodes),Yo(tip_nodes),K2n(:));
            fclose(fid);

            disp(sprintf('Mode I SIF : %10.3e Pa.sqrt(m)',K1o));
            disp(sprintf('Mode II SIF : %10.3e Pa.sqrt(m)',K2o));
        end
        %         UK((Nddl+1):length(U))=U((Nddl+1):length(U));
        %         load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix','sizeim');
        %         umap=-reshape((phix*UK).*mask,sizeim);
        %         umap=umap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
        %         PlotMap(umap,'U_{enr}',[],[filres,'-u-enr']);
        %         umap=-reshape((phix*(U-UK)).*mask,sizeim);
        %         umap=umap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
        %         PlotMap(umap,'U_{cla}',[],[filres,'-u-cla']);
        %         clear phix
        %
        %         load(fullfile('TMP',[num2str(nmod),'_phiy_0']),'phiy','sizeim');
        %         umap=-reshape((phiy*UK).*mask,sizeim);
        %         umap=umap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
        %         PlotMap(umap,'V_{enr}',[],[filres,'-v-enr']);
        %         umap=-reshape((phiy*(U-UK)).*mask,sizeim);
        %         umap=umap(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1);
        %         PlotMap(umap,'V_{cla}',[],[filres,'-v-cla']);
        %         clear phiy
        found=find(ind==-1);
        if any(found)
            assert(strcmp(model1.enrichment_type,'zone'))
            SK1=scalamp(found)*U(Nddl+found);
            %            K1o=3.e6;
            disp(sprintf('Crack tip shift : %6.2f pixel %10.3e m',-2*SK1/K1o/pix2m,-2*SK1/K1o));
        end

    end
end
end
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


