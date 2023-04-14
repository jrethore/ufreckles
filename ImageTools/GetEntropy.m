function [S]=GetEntropy(meshfile,im0,nflag)
if nargin<3,nflag=0;end
[pp,filname,ext]=fileparts(meshfile);
if isempty(ext)
    meshfile=[meshfile,'.mat'];
end
if nargin<2
    load(meshfile,'-mat','param');
    
    if ~isfield(param,'deformed_image')
        reader=VideoReader(param.reference_image);
        im0=read(reader,1);
    else
        im0=readim(param.reference_image);
    end
    if length(size(im0))==3
        im0=mean(im0,3);
    end
    roi=param.roi;
    im0=(im0(roi(1):roi(2),roi(3):roi(4)));
    
end
rflag=0;
if max(im0(:))<256
    nvg=0:255;
else
    nvg=0:(2^16-1);
end
load(meshfile,'-mat','rflag','xo','yo','Nnodes','Nelems','conn','elt');


S=zeros(prod(Nelems),1);

for i1=1:prod(Nelems)
    inods=conn(i1,1:elt(i1));
    
    xn=xo(inods);
    yn=yo(inods);
    jpix=ceil(min(yn)):floor(max(yn));
    ipix=ceil(min(xn)):floor(max(xn));
    [ypix,xpix]=meshgrid(jpix,ipix);
    im0e=im0(ipix,jpix);
    if ~rflag
        [xg,yg,wg]=GetGaussPointsPixels(elt(i1),xn,yn,xpix(:),ypix(:));
    else
        xg=-1+2*(xpix(:)-min(xn))/(max(xn)-min(xn));
        yg=-1+2*(ypix(:)-min(yn))/(max(yn)-min(yn));
        switch elt(i1)
            case 3
                wg=~(xg<0|yg<0|1-xg-yg<0);
            case 4
                wg=~(abs(xg)>1|abs(yg)>1);
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

end