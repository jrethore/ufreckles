function InitializeVICContour(nmod,iscale,h,refine)
check=1;
if nargin<3,h=50;end
if nargin<4,refine=0;end
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
if iscell(param0.reference_image)
    nim=length(param0.reference_image);
else
    nim=1;
end
if nim==1
    fildef=param0.reference_image;
else
    fildef=param0.reference_image{1};
end
im1=double(readim(fildef));
if numel(size(im1))==3
    im1=mean(im1,3);
end
if reverse
    im1=im1';
end
mim=0;
if isfield(param,'min_grey_level')
    mim=param.min_grey_level;
end
dim=255;
if isfield(param,'max_grey_level')
    dim=(param.max_grey_level)-mim;
end
roi=param0.roi;
sizeim=[roi(2)-roi(1)+1,roi(4)-roi(3)+1];
tau=param0.transition_length;
tau=tau*2^(iscale-1);
thickness=0;
if isfield(param,'line_thickness')
    thickness=round(0.5*(param.line_thickness));
end
[ls1,lso]=meshgrid(0:h,(-tau-thickness):(tau+thickness));
xc=reshape(param.starting_point,2,1);
yo=xc(2)-roi(3)+1+h*[0;1];
xo=xc(1)*ones(2,1)-roi(1)+1;
    switch param.contour_type
        case 'edge'
            im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
        case 'line'
            im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
    end
sizeim0=size(im0);
[gt,gn]=gradient(im0);
gn=diag(sparse(gn(:)));

npt=1;
if check
figure
hold on
imagesc(im1');
colormap(gray)
axis xy;
axis image;
end
continu=1;
while continu
    npt=length(xo);
    tx=xo(npt)-xo(npt-1);
    ty=yo(npt)-yo(npt-1);
    nnorm=abs(tx+i*ty);
    tx=tx/nnorm;ty=ty/nnorm;
    nx=-ty;ny=tx;
    erm=[];
    zo=ls1(:)+i*lso(:);
    thet=(-135:5:135)*pi/180;
    for it=1:length(thet)
        zoi=zo*exp(i*thet(it));
        Xi=xo(npt-1)+real(zoi)*tx+imag(zoi)*nx;
        Yi=yo(npt-1)+real(zoi)*ty+imag(zoi)*ny;
        disc=mexInterpLinear(Xi-1+roi(1),Yi-1+roi(3),im1);
        disc=mean(abs(im0(:)-disc(:)));
        erm=[erm,disc];
    end
    thet=thet(erm==min(erm));
    xo(npt)=xo(npt-1)+h*(cos(thet)*tx+sin(thet)*nx);
    yo(npt)=yo(npt-1)+h*(cos(thet)*ty+sin(thet)*ny);
    zoi=zo*exp(i*thet);
    Xi=xo(npt-1)+real(zoi)*tx+imag(zoi)*nx;
    Yi=yo(npt-1)+real(zoi)*ty+imag(zoi)*ny;
    tx=xo(npt)-xo(npt-1);
    ty=yo(npt)-yo(npt-1);
    nnorm=abs(tx+i*ty);
    tx=tx/nnorm;ty=ty/nnorm;
    nx=-ty;ny=tx;

    if npt==2
        phi=[ones(numel(im0),1),ls1(:)];
    else
        phi=[ls1(:)];
    end
    U=zeros(size(phi,2),1);
    phix=phi*nx;
    phiy=phi*ny;
    phidf=gn*phi;
    M=phidf'*phidf;
    res=1;ii=1;
    while (res>1.e-3)&&(ii<20)
        Ux=phix*U;
        Uy=phiy*U;
        disc=mexInterpLinear(Xi-1+roi(1)+Ux,Yi-1+roi(3)+Uy,im1);
        disc=(im0(:)-disc(:));
        F=phidf'*disc;
        dU=M\F;
        U=U+dU;
        res=norm(dU)/norm(U);
        ii=ii+1;


    end

    if npt==2
        xo=xo+U(1)*nx;
        yo=yo+U(1)*ny;
        xo(2)=xo(2)+U(2)*nx;
        yo(2)=yo(2)+U(2)*ny;
    else
        xo(npt)=xo(npt)+U*nx;
        yo(npt)=yo(npt)+U*ny;

    end
    if param.closed
        continu=(abs(xo(npt)-xo(1)+i*(yo(npt)-yo(1))))>h;
        if ~continu
            xo(npt)=xo(1);yo(npt)=yo(1);
        end
    else
        xnew=xo(npt)+h*tx;
        ynew=yo(npt)+h*ty;
        continu=(ynew>1)&&(xnew>1)&&(ynew<sizeim(2))&&(xnew<sizeim(1));
    end
    if check
    plot(xo+roi(1)-1,yo+roi(3)-1,'bx-','LineWidth',3);
    pause(0.01)
    end
    if continu
        xo=[xo;xo(npt)+h*tx];
        yo=[yo;yo(npt)+h*ty];
    end
end
if check
  axis off
  pause(0.01)
  print ('-djpeg', fullfile('FIG',sprintf('%s-init.jpg',param0.result_file)));
   
end
npt=length(xo);

si=(0:(h*(npt-1)-1))'+0.5;
so=0:h:(h*(npt-1));
yon=interp1(so,yo,si,'linear');
xon=interp1(so,xo,si,'linear');


yo=so;
type_nurbs='periodic';
if isfield(param,'continuity')
    type_nurbs=param.continuity;
end
degree=2;
if isfield(param,'degree')
degree=param.degree;
end
    vo=yo;
    switch type_nurbs
        % uniform open knot vector
        case 'open'
            for id=1:degree
                vo=[vo(1),vo,vo(length(vo))];
            end
            % periodic knot vector
        case 'periodic'
            for id=1:degree
                vo=[vo(1)-h,vo,vo(length(vo))+h];
            end

        otherwise
            error('Invalid type of nurbs basis')
    end
mesh_size=[h,h];
    Nnodes=[length(yo),1,1];
    Nelems=max(Nnodes-1,1);
    meshfile=fullfile('TMP',sprintf('%d_vicmesh_%d',nmod,iscale-1));
    save(meshfile,'yo','vo','Nnodes','Nelems','mesh_size');
[phi,dphi]=CreateNURBSBasis0D(meshfile,degree,si,type_nurbs,param.closed);
M=phi'*phi;
    L=phi'*xon;
    xs=M\L;
    L=phi'*yon;
    ys=M\L;
    xon=phi*xs;
    yon=phi*ys;
    dxon=dphi*xs;
    dyon=dphi*ys;
    nxon=dyon;
    nyon=-dxon;
    nnorm=abs(nxon+i*nyon);
    nxon=nxon./nnorm;
    nyon=nyon./nnorm;
[lso,ls1,nmesh]=ComputeLevelSetFromPoints(nmod,iscale,xon,yon,nxon,nyon,refine);
    nx=FDgradient(lso,1);
    ny=FDgradient(lso,2);
    tau=param0.transition_length;

    switch param.contour_type
        case 'edge'
            im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
            nband=find((lso(:)>=0)&(lso(:)<=tau)&nmesh);
        case 'line'
            im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
            nband=find((abs(lso(:))<=(tau+thickness))&nmesh);
    end

    on=find((lso(nband)<1)&(lso(nband)>=0));
            im1=(im1(roi(1):roi(2),roi(3):roi(4)));
mean0=mean(im1(:));
std0=1.5*std(im1(:));
        save(fullfile('TMP','sample0_0'),'im0','nx','ny','lso','ls1','sizeim','nband','on','mean0','std0');
    for iscale=2:param.nscale
        tau=2*tau;

    switch param.contour_type
        case 'edge'
            im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
            nband=find((lso(:)>=0)&(lso(:)<=tau)&nmesh);
        case 'line'
            im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
            nband=find((abs(lso(:))<=(tau+thickness))&nmesh);
    end

        on=find((lso(nband)<1)&(lso(nband)>=0));
        save(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim','nband','on','mean0','std0');
    end

    
    
    
    
end