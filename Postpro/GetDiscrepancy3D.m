function [disc,dynamic]=GetDiscrepancy3D(U,nmod,zone,ijm,err)
if nargin<4,ijm=0;end
if nargin<5,err=1;end
iscale=1;
load(fullfile('TMP','params'),'param');
param0=param;
roi=param0.roi;
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
dynamic=255;
if err
    switch param0.stack_format
        case 'bin'
            fid=fopen(fullfile('TMP',sprintf('dsample0_%d',iscale-1)));
            im0=fread(fid,prod(sizeim));
            fclose(fid);
            im0=reshape(im0,sizeim);
        case 'mat'
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
    end
    
    dynamic=double(max(im0(:))-min(im0(:)));
end
mesh_file=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1));
load(mesh_file,'rflag','xo','yo','zo','Nnodes','Nelems','Smesh','conn','elt','ng','ns');
invmap=~rflag;
Nn=prod(Nnodes);
disc=zeros([zone(2)-zone(1)+1,zone(4)-zone(3)+1,zone(6)-zone(5)+1]);
fildef=param0.deformed_image;
if iscell(fildef)
    if ijm==0
        error('I NEED THE IMAGE NUMBER');
    else
        if isfield(param0,'stack_size')
            imsiz0=param0.stack_size;
            fid=fopen(fildef{ijm},'r');
            jm3=fread(fid,prod(imsiz0));
            jm3=reshape(jm3,imsiz0);
            fclose(fid);
        else
            [~, ~, ext] = fileparts(fildef{ijm});
            switch ext
                case '.mat'
                    load(fildef{ijm},'jm3');
                case {'.tif','.tiff'}
                    jm3=readTIFFasRAW(fildef{ijm});
            end
        end
        if size(U,2)>1
            U=U(:,ijm);
        end
    end
else
    if isfield(param0,'stack_size')
        imsiz0=param0.stack_size;
        fid=fopen(fildef,'r');
        jm3=fread(fid,prod(imsiz0));
        jm3=reshape(jm3,imsiz0);
        fclose(fid);
    else
        [~, ~, ext] = fileparts(fildef);
        switch ext
            case '.mat'
                load(fildef,'jm3');
            case {'.tif','.tiff'}
                jm3=readTIFFasRAW(fildef);
        end
    end
end
im1=(jm3);
sizeim1=size(im1);
clear jm3

for i1=1:prod(Nelems)
    
    inods=conn(i1,1:elt(i1));
    xn=xo(inods);
    yn=yo(inods);
    zn=zo(inods);
    inbox=any(xn+roi(1)-1>zone(1))&&any(xn+roi(1)-1<zone(2))&&any(yn+roi(3)-1>zone(3))&&any(yn+roi(3)-1<zone(4))&&any(zn+roi(5)-1>zone(5))&&any(zn+roi(5)-1<zone(6));
    
    if inbox
        Un=U(inods+0*Nn);
        Vn=U(inods+1*Nn);
        Wn=U(inods+2*Nn);
        
        ipix=max(1,ceil(min(xn)-1)):min(sizeim(1),floor(max(xn)+1));
        jpix=max(1,ceil(min(yn)-1)):min(sizeim(2),floor(max(yn)+1));
        kpix=max(1,ceil(min(zn)-1)):min(sizeim(3),floor(max(zn)+1));
        [ypix,xpix,zpix]=meshgrid(jpix,ipix,kpix);
        if invmap
            [xg,yg,zg,wg]=GetGaussPointsVoxels(elt(i1),xn,yn,zn,xpix(:),ypix(:),zpix(:));
        else
            xg=-1+2*(xpix-min(xn))/(max(xn)-min(xn));
            yg=-1+2*(ypix-min(yn))/(max(yn)-min(yn));
            zg=-1+2*(zpix-min(zn))/(max(zn)-min(zn));
            switch elt(i1)
                case 6
                    wg=~(xg(:)<0|yg(:)<0|1-xg(:)-yg(:)<0|abs(zg(:))>1);
                case 4
                    wg=~(xg(:)<0|yg(:)<0|zg(:)<0|1-xg(:)-yg(:)-zg(:)<0);
                case 8
                    wg=~(abs(xg(:))>1|abs(yg(:))>1|abs(zg(:))>1);
            end
        end
        N=GetFiniteElementShapeFunctions3D(elt(i1),xg(:),yg(:),zg(:));
        if elt(i1)<8
            N=N(:,1:elt(i1));
        end
        xi=xpix(:)+(N*Un)+roi(1)-1;
        yi=ypix(:)+(N*Vn)+roi(3)-1;
        zi=zpix(:)+(N*Wn)+roi(5)-1;
        xmin=max(floor(min(xi))-2,1);xmax=min(ceil(max(xi))+2,size(im1,1));
        ymin=max(floor(min(yi))-2,1);ymax=min(ceil(max(yi))+2,size(im1,2));
        zmin=max(floor(min(zi))-2,1);zmax=min(ceil(max(zi))+2,size(im1,3));
        im1el=double(im1(xmin:xmax,ymin:ymax,zmin:zmax));
        im1e=mexInterpLinear3D(xi-xmin+1,yi-ymin+1,zi-zmin+1,im1el);
        
        
        xi=xpix(:)+roi(1)-1;
        yi=ypix(:)+roi(3)-1;
        zi=zpix(:)+roi(5)-1;
        inboxi=(wg)&(xi>=zone(1))&(xi<=zone(2))&(yi>=zone(3))&(yi<=zone(4))&(zi>=zone(5))&(zi<=zone(6));
        found=find(inboxi);
        ind=sub2ind(size(disc),xi(found)-zone(1)+1,yi(found)-zone(3)+1,zi(found)-zone(5)+1);
        if err
            im0e=double(im0(ipix,jpix,kpix));
            disc(ind)=abs(im1e(found)-im0e(found));
        else
            disc(ind)=im1e(found);
        end
    end
    
end
if err
    disc=100*disc/dynamic;
end
end