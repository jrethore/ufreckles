function [disc,dynamic,edisc]=GetDiscrepancy2D(U,nmod,zone,ijm,scal)
if nargin<4,ijm=0;end
if nargin<5,scal=0;end
iscale=1;
load(fullfile('TMP','params'),'param');
param0=param;
roi=param0.roi;

load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
inorm=false;
if isfield(param0,'normalize_grey_level')
    inorm=param0.normalize_grey_level;
end

load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
dynamic=max(im0(:))-min(im0(:));

[mean0,std0]=mexImScalling2D(im0);
mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
load(mesh_file,'rflag','xo','yo','Nnodes','Nelems','Smesh','conn','elt','ng','ns');
invmap=~rflag;
Nn=prod(Nnodes);
edisc=zeros(prod(Nelems),1);
disc=zeros([zone(2)-zone(1)+1,zone(4)-zone(3)+1]);
videok=~isfield(param0,'deformed_image');
if videok
    reader=VideoReader(param0.reference_image);
    if ijm==0
        error('I NEDD THE IMAGE NUMBER');
    else
        if size(U,2)>1
            U=U(:,ijm);
        end
    end
    nbf=reader.NumberOfFrames-1;
    if isfield(param0,'number_of_frames')
        nbf=param0.number_of_frames;
    end
    dim=1;
    if isfield(param0,'video_sampling')
        dim=param0.video_sampling;
    end
    frames=2:dim:nbf;
    im1=readim(reader,frames(ijm));
else
    fildef=param0.deformed_image;
    if iscell(fildef)
        if ijm==0
            error('I NEDD THE IMAGE NUMBER');
        else
            fildef=fildef{1,ijm};
            if size(U,2)>1
                U=U(:,ijm,1);
            end
        end
    end
    im1=double(readim(fildef));
end
if length(size(im1))==3
    im1=mean(im1,3);
end
sizeim1=size(im1);

[mean1,std1]=mexImScalling2D(im1(roi(1):roi(2),roi(3):roi(4)));


for i1=1:prod(Nelems)
    
    inods=conn(i1,1:elt(i1));
    xn=xo(inods);
    yn=yo(inods);
    inbox=any(xn+roi(1)-1>zone(1))&&any(xn+roi(1)-1<zone(2))&&any(yn+roi(3)-1>zone(3))&&any(yn+roi(3)-1<zone(4));
    
    if inbox
        Un=U(inods+0*Nn);
        Vn=U(inods+1*Nn);
        ipix=max(1,ceil(min(xn)-2)):min(sizeim(1),floor(max(xn)+2));
        jpix=max(1,ceil(min(yn)-2)):min(sizeim(2),floor(max(yn)+2));
        im0e=im0(ipix,jpix);
        [ypix,xpix]=meshgrid(jpix,ipix);
        if invmap
            [xg,yg,wg]=GetGaussPointsPixels(elt(i1),xn,yn,xpix(:),ypix(:));
        else
            xg=-1+2*(xpix-min(xn))/(max(xn)-min(xn));
            yg=-1+2*(ypix-min(yn))/(max(yn)-min(yn));
            switch elt(i1)
                case 3
                    wg=~(xg(:)<0|yg(:)<0|1-xg(:)-yg(:)<0);
                    
                case 4
                    
                    wg=~(abs(xg(:))>1|abs(yg(:))>1);
            end
        end
        
        N=GetFiniteElementShapeFunctions(elt(i1),xg(:),yg(:));
        
        if elt(i1)<4
            N=N(:,1:elt(i1));
        end
        xpix=round(N*xn);
        ypix=round(N*yn);
        xi=xpix+roi(1)-1+N*(Un);
        yi=ypix+roi(3)-1+N*(Vn);
        
        im1e=mexInterpSpline(xi,yi,im1);
        maske=(im1e<0)|abs(im1e)>100*mean0|isnan(im1e);
        im1e(maske)=im0e(maske);
        
        if scal
            if inorm
                if any(maske)
                    im1e=im1e-mean(im0e(:));
                    sc=1;
                else
                    im1e=im1e-mean(im1e(:));
                    sc=max(1,std(im0e(:)))/max(1,std(im1e(:)));
                end
                im0e=im0e-mean(im0e(:));
                im1e=sc*im1e;
            else
                
                im0e=im0e-mean0;
                im1e=im1e-mean1;
                sc=max(1,std0)/max(1,std1);
                im1e=sc*im1e;
            end
        end
        
        
        
        xi=xpix(:)+roi(1)-1;
        yi=ypix(:)+roi(3)-1;
        inboxi=(wg)&(xi>=zone(1))&(xi<=zone(2))&(yi>=zone(3))&(yi<=zone(4));
        ind=sub2ind(size(disc),xi(inboxi)-zone(1)+1,yi(inboxi)-zone(3)+1);
        disc(ind)=abs(im1e(inboxi)-im0e(inboxi));
        edisc(i1)=sqrt(mean((im1e(inboxi)-im0e(inboxi)).^2));
        %        edisc(i1)=(mean(abs(im1e(inboxi)-im0e(inboxi))));
    end
    
end

disc=100*disc/dynamic;
edisc=100*edisc/dynamic;

end