function postproVIC(filreso,comp)
if nargin<2, comp=false;end
disp(['Result file:',filreso]);
lim=[];
epsflag=true;
nmod=1;
load([filreso])
filexp='';
if ~exist('model1')
model1=model;
end
roi=param.roi;
reverse=0;
if isfield(param,'reverse_image')
    reverse=param.reverse_image;
end
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
fac=1;
submean=0;
% if strcmp(param.analysis,'correlation')
%     if exist([filreso,'-error.mat'],'file')
%     load([filreso,'-error.mat'],'disc','dynamic');
%     else
%     load(fullfile('TMP',[num2str(nmod),'_error_0']),'disc','dynamic');
%     end
% end
%% mask
tau=param.transition_length;

load(fullfile('TMP',[num2str(nmod),'_mask_0']),'mask');
load(fullfile('TMP',[num2str(nmod),'_phi_0']),'phi','sizeim');
mask=double(diag(mask)>0);
found=find(mask==0);
mask=full(mask);
mask(found)=NaN;
% xo=reshape(xo,Nnodes)+0.5;
% yo=reshape(yo,Nnodes)+0.5;
% xo=xo(1:Nnodes(1),1)';
% yo=yo(1,1:Nnodes(2));
% xi= xo(1):xo(length(xo))-1;
% yi=yo(1):yo(length(yo))-1;
xi= 1:sizeim(1);
yi=1:sizeim(2);
Ut=U;
load(fullfile('TMP',sprintf('%d_phi_%d',nmod,(1-1))),'phii','dphii','ddphii','ls1i');
[xon,yon]=ind2sub(sizeim,nband(on));
nxon=-nx(nband(on));
nyon=-ny(nband(on));
xono=xon;
yono=yon;
for iim=1:length(images)
    U=Ut(:,iim);
    if length(images)>1
        filexp=sprintf('-%03d',images(iim));
    end
    filres=[filreso,filexp];
%% displacement map
            if length(images)==1
                fildef=param.deformed_image;
                im1=double(readim(fildef));
           else
                fildef=param.deformed_image{iim};
                im1=double(readim(fildef));
            end
                        if length(size(im1))==3
                im1=mean(im1,3);
            end
    if reverse
                    im1=im1';
     end
 ux=nxon.*(phii*U+0.5*(1*tau-tau));
 uy=nyon.*(phii*U+0.5*(1*tau-tau));
 DeformedImageAndContour(ux+xono+roi(1)-1,uy+yono+roi(3)-1,im1,[filres,'-cont-im'],epsflag,1);
clear umap vmap
xon=xono+ux+roi(1)-1;
yon=yono+uy+roi(3)-1;
save([filres,'-xyon'],'xon','yon');
%% errormap
% if strcmp(param.analysis,'correlation')
% 
%     load(fullfile('TMP',[num2str(nmod),'_error_0']),'disc','dynamic');
%     disci=reshape(disc(:,iim),sizeim);
%     disci=100*abs(disci(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1))/dynamic;
%     if comp
%         limd=[0,10];
%     else
%     limd=[0,0.5*100*max(disc(:))/dynamic];
%     end
%    PlotMap(disci,'Discrepancy wrt dyn (%)',limd,[filres,'-error'],epsflag,0);
%     disp(sprintf('Mean error wrt image dyn.: %6.2f %%',mean(disci(:))));
% end
if length(images)>1&&iim<length(images)
close all
pause(0.1)
end




 end
%     type_nurbs=~strcmp(model1.continuity,'c0');
%     pscale=1;
%     p=model1.degree;
%     dy=(yo(length(yo))-yo(1))/100;
% yi=yo(1):dy:yo(length(yo));
%     [phi1y,dphi1y,ddphi1y]=CreateNURBSBasis0D(filreso,p,yi,type_nurbs,model1.closed);

%     figure
%     plot(yi,ddphi1y)
[yi,ind]=sort(ls1(nband(on)));
        umap=(phii(ind,:)*Ut)'+0.5*(1*tau-tau);
    dumap=(dphii(ind,:)*Ut)';
    ddumap=(ddphii(ind,:)*Ut)';
    
if submean
  M=[repmat(1,size(umap,2),1),(1:size(umap,2))'];
  Um=M\umap';
umap=umap-(M*Um)';
end
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
end


