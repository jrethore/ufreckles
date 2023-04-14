function run_uni_job(paramo,modelo)
global phix phiy
nmod=0;
roio=paramo.roi;
param=paramo;
filres=paramo.result_file;
filreso=strrep(filres,'.res','');
gages=modelo.zone(2,:);
param.onflight=1;
model.nscale=modelo.nscale;
model.basis='uni';
xo=zeros(4*length(gages),1);
yo=zeros(4*length(gages),1);
elt=4*ones(1*length(gages),1);
conn=zeros(length(elt),4);
Nnodes=[length(xo),1,1];
Nelems=[length(elt),1,1];
nim=1;
if iscell(param.deformed_image)
nim=size(param.deformed_image,2);
end
U=zeros(2*prod(Nnodes),nim);
for ig=1:length(gages)
    filresi=[filreso,sprintf('-%02d',ig)];
    param.result_file=filresi;
    gage=gages{ig};
    xp=gage(1:4,1)+roio(1)-1;
    yp=gage(1:4,2)+roio(3)-1;
    roi(1:2) = [(floor(min(xp))-1),(ceil(max(xp))+1)];
    roi(3:4) = [(floor(min(yp))-1),(ceil(max(yp))+1)];
    sizeim=[roi(2)-roi(1),roi(4)-roi(3)]+1;
    param.roi=roi;
    LoadParameters(param);
    LoadParameters(model,nmod);
    ReferenceImage(nmod);
    
    LoadMask(nmod)
    [Yi, Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
    mask0=inpolygon(Xi,Yi,xp([1:length(xp),1])-roi(1)+1,yp([1:length(xp),1])-roi(3)+1);
    for iscale=1:model.nscale
        if iscale>1
            mask0=CoarseImage(mask0);
            mask0=mask0==1;
        end
        mask=diag(sparse(mask0(:)));
        save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask','-append');
    end
    
    U1=[];
    
    for iscale=model.nscale:-1:1
        disp(sprintf('Pre-processing scale %d...',iscale));
        CreateBasisFunction(iscale,nmod);
        ComputeGradFPhi(iscale,nmod);
        AssembleCorrelationOperator(iscale,nmod);
        
        Uini=InitializeSolution(U1,iscale,nmod);
        
        
        
        [U1]=Solve(Uini,iscale,nmod);
    end
    copyfile(fullfile('TMP','0_error_0.mat'), [filresi,'-error.res']);
    delete([fullfile('VTK',[filresi,'-error']),'-0*.vtk']);
    images=1:size(U1,2);
    for iim=1:size(U1,2)
        movefile(fullfile('VTK',sprintf('camr-1-camd-1-scale-1-%d-error.vtk',iim)),fullfile('VTK',sprintf('%s-error-%04d.vtk',filresi,images(iim))));
    end
    VTKExportScalarMap([filresi,sprintf('-%d',0)],'error',roi(1)-1+(1:sizeim(1)),roi(3)-1+(1:sizeim(2)),zeros(sizeim),1);
    movefile(fullfile('VTK',[filresi,sprintf('-%d',0),'-error.vtk']),fullfile('VTK',sprintf('%s-error-%04d.vtk',filresi,0)));
    delete(fullfile('VTK',['camr*','-error.vtk']));
%     xi=round([xp(:);mean(xp)]);
%     yi=round([yp(:);mean(yp)]);
%     xo((1:5)+(ig-1)*5)=xi;
%     yo((1:5)+(ig-1)*5)=yi;
%     for ie=1:4
%         conn(ie+(ig-1)*4,1:3)=[mod(ie-1,4)+1,mod(ie+1-1,4)+1,5]+(ig-1)*5;
%     end
    xi=round([xp(:)]);
    yi=round([yp(:)]);
    xo((1:4)+(ig-1)*4)=xi;
    yo((1:4)+(ig-1)*4)=yi;
   conn(ig,1:4)=(1:4)+(ig-1)*4;
    indp=xi-roi(1)+1+(yi-roi(3)+1-1)*sizeim(1);
    phixo=phix(indp,:);
    Ux=phixo*U1;
    phiyo=phiy(indp,:);
    Uy=phiyo*U1;
    U((1:4)+(ig-1)*4,:)=Ux;
    U((1:4)+(ig-1)*4+prod(Nnodes),:)=Uy;
                exx=U1(4,:)/mean(sizeim);
                eyy=U1(5,:)/mean(sizeim);
                exy=U1(6,:)/mean(sizeim);
    
    Up{ig}=[U1;exx;eyy;exy];
    
    l1=abs(diff(gage(1:2,:)*[1;1i]));
    l2=abs(diff(gage(2:3,:)*[1;1i]));
    angl=angle(diff(gage(1:2,:)*[1;1i]));
    if l2<l1,angl=angl+pi/2;end
    angl=exp(1i*angl);
    
    R=[real(angl*exp(1i*pi/2)),real(angl);...
        imag(angl*exp(1i*pi/2)),imag(angl)];
El=R(1,1)*(R(1,1)*exx+R(2,1)*exy)+R(2,1)*(R(1,1)*exy+R(2,1)*eyy);
Et=R(1,2)*(R(1,2)*exx+R(2,2)*exy)+R(2,2)*(R(1,2)*exy+R(2,2)*eyy);
Es=R(1,2)*(R(1,2)*exx+R(2,1)*exy)+R(2,1)*(R(1,1)*exy+R(2,1)*eyy);

                    fid=fopen([filresi,'-data.csv'],'w');
fprintf(fid,'Result file;%s\n',param.result_file);
fprintf(fid,'Reference image;%s\n',param.reference_image);
fprintf(fid,'Strain gage;%d\n',ig);
fprintf(fid,'Position X [pixel];Position Y [pixel]; Orientation; Length [pixel]; Width [pixel]\n');
fprintf(fid,'%f;%f;%f;%f;%f\n',mean(xp),mean(yp),(angle(angl)+pi/2)*180/pi,max(l1,l2),min(l1,l2));
                fprintf(fid,'Filename;Image;Ux [pixel];Uy [pixel]; Rot [];Exx [];Eyy [];Exy[];El [];Et [];Es[]\n');
if size(U1,2)==1
    iim=1;
                    fprintf(fid,'"%s";%d;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e\n',param.deformed_image,iim,U1(1:2,iim),U1(3:6,iim)/mean(sizeim),El(iim),Et(iim),Es(iim));
else
    for iim=1:size(U1,2)
                    fprintf(fid,'"%s";%d;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e\n',param.deformed_image{iim},iim,U1(1:2,iim),U1(3:6,iim)/mean(sizeim),El(iim),Et(iim),Es(iim));
    end
end
                fclose(fid);

end
    rflag=1;
    rint=false;
    ng=0;
    param=paramo;
    model=modelo;
            save(filres,'U','Up','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
%ExportImageToVTK(filres)
postproVTK(filres,0,0);


end