function UpdateMesh(Ui,iscale,nmod,ijm)
ncam=1;
assert(iscale==1);
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
roi=param0.roi;
load(fullfile('TMP','sample0'),'sizeim');
pscale=2^(iscale-1);
mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
load(mesh_file,'xo','yo');
xo=xo+Ui((1:length(xo)));
yo=yo+Ui(length(xo)+(1:length(xo)));


    filref=param0.deformed_image{ncam,ijm};
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end

    im0=double(readim(filref));
    sizeim=size(im0);
    if length(sizeim)==3
        im0=mean(im0,3);
    end        
    if reverse
        im0=im0';
    end
    sizeim=size(im0);
% figure
% imagesc(im0')
%         colormap('gray')
% 
%         hold on
% plot(xo+roi(1)-1,yo+roi(3)-1,'bo')    
% 
% pause

roi(2)=min(sizeim(1),roi(1)+ceil(max(xo))+5);
roi(1)=max(1,roi(1)+floor(min(xo))-5);
roi(4)=min(sizeim(2),roi(3)+ceil(max(yo))+5);
roi(3)=max(1,roi(3)+floor(min(yo))-5);
    save(fullfile('TMP',sprintf('sample%d',ncam-1)),'im0','sizeim','roi');
    im0=(im0(roi(1):roi(2),roi(3):roi(4)));
sizeim=size(im0);
    xo=max(xo,1);
xo=min(xo,sizeim(1));
yo=max(yo,1);
yo=min(yo,sizeim(2));
save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','roi','-append');
save(fullfile('TMP',sprintf('sample%d_%d',(ncam-1),0)),'im0','sizeim','roi');
param=param0;
param.roi=roi;
save(fullfile('TMP','params'),'param');
CreateBasisFunction(iscale,nmod);
ComputeGradFPhi(iscale,nmod);
AssembleCorrelationOperator(iscale,nmod);




end