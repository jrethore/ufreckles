function [S]=GetEntropy3D(meshfile,im0,nflag)
if nargin<3,nflag=0;end
[pp,filname,ext]=fileparts(meshfile);
if isempty(ext)
    meshfile=[meshfile,'.mat'];
end
if nargin<2
    load(fullfile('TMP',sprintf('sample0_%d',1-1)),'sizeim');
    fid=fopen(fullfile('TMP',sprintf('dsample0_%d',1-1)));
    if fid>-1
        im0=fread(fid,prod(sizeim));
        fclose(fid);
        im0=reshape(im0,sizeim);
    else
        load(fullfile('TMP',sprintf('sample0_%d',1-1)),'im0');
    end
end
rflag=0;
if max(im0(:))<256
    nvg=0:255;
else
    nvg=0:(2^16-1);
end
load(meshfile,'-mat','rflag','xo','yo','zo','Nnodes','Nelems','conn','elt');


S=zeros(prod(Nelems),1);

for i1=1:prod(Nelems)
    inods=conn(i1,1:elt(i1));
    
    xn=xo(inods);
    yn=yo(inods);
    zn=zo(inods);
    jpix=ceil(min(yn)):floor(max(yn));
    ipix=ceil(min(xn)):floor(max(xn));
    kpix=ceil(min(zn)):floor(max(zn));
    [ypix,xpix,zpix]=meshgrid(jpix,ipix,kpix);
    im0e=im0(ipix,jpix,kpix);
    if ~rflag||~(elt(i1)==8)
        [ypix,xpix,zpix]=meshgrid(jpix,ipix,kpix);
        [xg,yg,zg,wg]=GetGaussPointsVoxels(elt(i1),xn,yn,zn,xpix(:),ypix(:),zpix(:));
    else
        xg=-1+2*(xpix(:)-min(xn))/(max(xn)-min(xn));
        yg=-1+2*(ypix(:)-min(yn))/(max(yn)-min(yn));
        zg=-1+2*(zpix(:)-min(zn))/(max(zn)-min(zn));
        switch elt(i1)
            case 6
                wg=~(xg<0|yg<0|1-xg-yg<0|abs(zg)>1);
            case 4
                wg=~(xg<0|yg<0|zg<0|1-xg-yg-zg<0);
            case 8
                wg=~(abs(xg)>1|abs(yg)>1|abs(zg)>1);
        end
    end
    
    npix=sum(wg>0);
    im0e=im0e(wg>0);
    
    p=(histc(im0e,nvg))/npix;
    p=p(p>0);
    if nflag
        smax=log2(npix)/max(nvg);
        S(i1)=-sum(p.*log2(p))/max(nvg)/smax;
    else
        S(i1)=-sum(p.*log2(p))/max(nvg);
    end
end
S=reshape(S,Nelems);
end