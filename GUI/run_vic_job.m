function run_vic_job(paramo,modelo)
nmod=0;
param=paramo;
filres=paramo.result_file;
filreso=strrep(filres,'.res','');
nscale=modelo.nscale;
filo=param.reference_image;
[~,~,ext]=fileparts(filo);
if iscell(param.deformed_image)
    fildefs=param.deformed_image;
    fildefs=[filo,fildefs];
    images=1:size(fildefs,2);
else
    if strcmp(filo,param.deformed_image)
        fildefs=filo;
        images=1;
    else
        fildefs={filo,param.deformed_image};
        images=1:size(fildefs,2);
    end
end
param.deformed_image=fildefs;
filim=strrep(param.result_file,'.res','-vi.mat');
param.reference_image=filim;
[Yi,Xi]=meshgrid(1:param.sizeim(2),1:param.sizeim(1));
xo=[];
yo=[];
conn=[];
Ux=[];
Uy=[];
dU=[];
ddU=[];
Up=[];
nmax=0;
curv=cell(size(modelo.zone,2),1);
fleche=cell(size(modelo.zone,2),1);
s=cell(size(modelo.zone,2),1);
for iz=1:size(modelo.zone,2)
    
    zone=modelo.zone(:,iz);
    switch zone{4}
        case 1
            xyp=zone{2};
            xyp=xyp*[1;1i];
            roi=round([min(real(xyp)),max(real(xyp)),min(imag(xyp)),max(imag(xyp))]);
            try
                model.starting_point=zone{10};
            catch
                model.starting_point=[real(xyp(1)),imag(xyp(1))];
            end
            mx=1;
            my=1;
            
            
            model.continuity='periodic';
        case 3
            xyp=zone{5};
            Zi=Xi+1i*Yi-(xyp(1)+1i*xyp(2));
            lso=abs(Zi);
            ls1=angle(Zi)*xyp(3);
            nx=-real(Zi)./(max(lso,1));
            ny=-imag(Zi)./(max(lso,1));
            lso=-(lso-xyp(3));
            roi=[1,param.sizeim(1),1,param.sizeim(2)];
            mx=1;my=1;
            model.continuity='periodic';
        case 4
            
            model.continuity='open';
            xyp=zone{2};
            xyp=xyp*[1;1i];
            Zi=Xi+1i*Yi-xyp(1);
            seg=diff(xyp([1;numel(xyp)]));
            n=(-1i*seg)/abs(seg);
            t=(seg)/abs(seg);
            lso=(real((Zi)'*n))';
            ls1=(real((Zi)'*t))';
            nx=ones(param.sizeim(1:2))*real(n);
            ny=ones(param.sizeim(1:2))*imag(n);
            roi=round([min(real(xyp)),max(real(xyp)),min(imag(xyp)),max(imag(xyp))]);
            mx=(abs(angle(n))<pi/4)||(abs(angle(n))>3*pi/4);
            my=(abs(angle(n))>pi/4)&&(abs(angle(n))<3*pi/4);
            
    end
    if zone{4}>1
        if zone{1}<0
            im=readim(filo);
            mim=mean(im((lso>0)&(lso<0.5*max(lso(:)))));
            mam=mean(im((lso<0)&(lso>-0.5*max(lso(:)))));
            if mim<mam
                lso=-lso;
                nx=-nx;ny=-ny;
            end
            model.contour_type='edge';
        else
            model.contour_type='line';
            model.line_thickness=zone{9};
        end
        save(filim,'lso','ls1','nx','ny','-v7.3');
    else
        model.contour_type='line';
        model.line_thickness=zone{9};
        param.reference_image=filo;
        
    end
    thickness=zone{9};
    tau=zone{8};
    roi([1])=max(roi([1])-mx*tau*2^(nscale-1),1);
    roi([3])=max(roi([3])-my*tau*2^(nscale-1),1);
    roi([2])=min(roi([2])+mx*tau*2^(nscale-1),param.sizeim(1));
    roi([4])=min(roi([4])+my*tau*2^(nscale-1),param.sizeim(2));
    param.roi=roi;
    param.transition_length=tau;
    param.onflight=0;
    LoadParameters(param);
    model.nscale=nscale;
    model.degree=zone{7};
    model.basis='vic-nurbs';
    model.Nelems=zone{6};
    model.closed=zone{4}==3;
    LoadParameters(model,nmod);
    
    if zone{4}==1
        InitializeVICContour(nmod,1,50,0);
    else
        if zone{4}==4&&size(xyp,1)>2
            
            %pour Rana
            si=cumsum(abs(diff(xyp)));
            si=[0;si(:)];
            
            ls1i=linspace(0,max(si(:)),1000);
            if 0
            xcoefs=polyfit(si,real(xyp),3);
            ycoefs=polyfit(si,imag(xyp),3);
            xi=polyval(xcoefs,ls1i);
            yi=polyval(ycoefs,ls1i);
            dxcoefs=(numel(xcoefs):-1:1).*xcoefs;
            dxcoefs=dxcoefs(1:end-1);
            dycoefs=(numel(ycoefs):-1:1).*ycoefs;
            dycoefs=dycoefs(1:end-1);
            dxi=polyval(dxcoefs,ls1i);
            dyi=polyval(dycoefs,ls1i);
            else
                xi=interp1(si,real(xyp),ls1i,'spline');
                 yi=interp1(si,imag(xyp),ls1i,'spline');
                dxi=gradient(xi);
                dyi=gradient(yi);
            end
            
            
            ni=(-1i*(dxi+1i*dyi));
            ni=ni(:)./abs(ni(:));
            
            xyon=xi(:)+1i*yi(:);
            figure
            plot(xyon)
            hold on
            plot(xyp)
            ini=inpolygon(real(xyon),imag(xyon),roi([1,2,2,1,1]),roi([3,3,4,4,3]));
            xyon=xyon(ini);
            ni=-ni(ini);
            [lson,ls1n,nmesh]=ComputeLevelSetFromPoints(nmod,1,real(xyon(:))-roi(1)+1,imag(xyon(:))-roi(3)+1,real(ni(:)),imag(ni(:)),0);
            lson(~nmesh)=100000;
            figure 
            imagesc(lson)
             lso(roi(1):roi(2),roi(3):roi(4))=lson;
            ls1(roi(1):roi(2),roi(3):roi(4))=ls1n;
            nxi=FDgradient(lson,1);
            nyi=FDgradient(lson,2);
            nx(roi(1):roi(2),roi(3):roi(4))=nxi;
            ny(roi(1):roi(2),roi(3):roi(4))=nyi;
            
            save(filim,'lso','ls1','nx','ny','-v7.3');
            
        end
        VirtualImage(nmod);
    end
    U1=[];
    VICMeshes(nmod);
    LoadMask(nmod);
    for iscale=(nscale):-1:1
        CreateVICBasisFunction(iscale,nmod);
        ComputeGradFPhi(iscale,nmod);
        AssembleCorrelationOperator(iscale,nmod);
        Uini=InitializeVICSolution(U1,iscale,nmod);
        [U1]=SolveVIC(Uini,iscale,nmod);
    end
    Up=[Up;U1];
    load(fullfile('TMP',sprintf('sample0_%d',0)),'lso','ls1','nx','ny','nband','on');
    load(fullfile('TMP',[num2str(nmod),'_phi_0']),'phii','dphii','ddphii','sizeim');
    [xon,yon]=ind2sub(sizeim,nband(on));
    ls1on=ls1(nband(on));
    nxon=-nx(nband(on));
    nyon=-ny(nband(on));
    [si,ind]=sort(ls1(nband(on)));
    ls1on=ls1on(ind);
    s{iz}=si;
    phii=phii(ind,:);
    dphii=dphii(ind,:);
    ddphii=ddphii(ind,:);
    fleche{iz}=phii*U1;
    curv{iz}=ddphii*U1;
    xon=xon(ind)+param.roi(1)-1;
    yon=yon(ind)+param.roi(3)-1;
    nxon=nxon(ind);
    nyon=nyon(ind);
    uxi=diag(sparse(nxon))*(phii*U1);
    uyi=diag(sparse(nyon))*(phii*U1);
    nn=length(xon);
    %xon=xon+uxi(:,1);
    %yon=yon+uyi(:,1);
    %uxi(:,1)=0;
    %uyi(:,1)=0;
    ddui=ddphii*U1;
    
    
    
    fid=fopen([filreso,sprintf('-zone-%02d.csv',iz)],'w');
    fprintf(fid,'Result file;%s\n',param.result_file);
    fprintf(fid,'Zone gage;%d\n',iz);
    fprintf(fid,'Position s [pixel];X [pixel];Y [pixel];Curvature [1/pixel]\n');
    fprintf(fid,'Image file;');
    if size(U1,2)>1
        for iim=1:size(U1,2)
            fprintf(fid,'"%s";;;',param.deformed_image{iim});
        end
    else
        fprintf(fid,'"%s";;;',param.deformed_image);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Image number;');
    for iim=1:size(U1,2)
        fprintf(fid,'%d;;;',iim);
    end
    fprintf(fid,'\n');
    for ix=1:length(xon)
        fprintf(fid,'%f;',ls1on(ix));
        for  iim=1:size(U1,2)
            fprintf(fid,'%f;%f;%f;',xon(ix)+uxi(ix,iim),yon(ix)+uyi(ix,iim),ddui(ix,iim));
        end
        fprintf(fid,'\n');
    end
    
    
    switch model.contour_type
        case 'line'
            xon=[xon-0.5*(tau+thickness)*nxon;xon+0.5*(tau+thickness)*nxon];
            yon=[yon-0.5*(tau+thickness)*nyon;yon+0.5*(tau+thickness)*nyon];
        case 'edge'
            xon=[xon-(tau+thickness)*nxon;xon+0*(tau+thickness)*nxon];
            yon=[yon-(tau+thickness)*nyon;yon+0*(tau+thickness)*nyon];
    end
    uxi=[uxi;uxi];
    uyi=[uyi;uyi];
    Ux=[Ux;uxi];
    Uy=[Uy;uyi];
    ddU=[ddU;ddui;ddui];
    xo=[xo;xon];
    yo=[yo;yon];
    conni=repmat([1,2,2+nn,1+nn],nn-1,1)+repmat((0:nn-2)',1,4);
    conn=[conn;conni+nmax];
    nmax=max(conn(:));
    
    
end
paramo.deformed_image=param.deformed_image;
param=paramo;
model=modelo;
Nnodes=[length(xo),1,1];
Nelems=[size(conn,1),1,1];
elt=repmat(4,prod(Nelems),1);
U=[Ux;Uy];
delete([fullfile('VTK',filreso),'-0*.vtk']);
ExportVTKField(filreso,images,xo,yo,elt,conn,{U,ddU},{'Displacement','Curvature'});
rint=false;
ng=0;
rflag=0;
save(filres,'U','Up','s','curv','fleche','ddU','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');

end
