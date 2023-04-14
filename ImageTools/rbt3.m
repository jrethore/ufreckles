function [dec]=rbt3(im0,im1,inorm)
if nargin<3,inorm=false;end
subpix=1;
nsiz=size(im0);
nsiz2=2.^floor(log2(nsiz+.5));
off=floor((nsiz-nsiz2)/2);
im0=im0(1+off(1):nsiz2(1)+off(1),1+off(2):nsiz2(2)+off(2),1+off(3):nsiz2(3)+off(3));
im1=im1(1+off(1):nsiz2(1)+off(1),1+off(2):nsiz2(2)+off(2),1+off(3):nsiz2(3)+off(3));
if inorm
    im0=(im0-mean(im0(:)))/std(im0(:));
    im1=(im1-mean(im1(:)))/std(im1(:));
    
end

ffim0=fftn(im0);ffim0(1,1,1)=0;
ffim1=fftn(im1);ffim1(1,1,1)=0;
ffinter=ffim0.*conj(ffim1);
inter=real(ifftn(ffinter));
shinter=fftshift(inter);
[mx,ind]=max(shinter(:));
[i1,i2,i3]=ind2sub(nsiz2,ind);
ii=[i1,i2,i3];
dec=-(ii-nsiz2/2-1)';
if (subpix==1)
    xi=min(max(i1+(-2:2),1),size(shinter,1));
    yi=min(max(i2+(-2:2),1),size(shinter,2));
    zi=min(max(i3+(-2:2),1),size(shinter,3));
    T=shinter(xi,yi,zi);
    [Y,X,Z]=meshgrid(yi-i2,xi-i1,zi-i3);
    One=ones(length(xi),length(yi),length(zi));
    X2=X.^2;Y2=Y.^2;Z2=Z.^2;
    XY=X.*Y;YZ=Y.*Z;ZX=Z.*X;
    M=[X2(:) Y2(:) Z2(:) XY(:) YZ(:) ZX(:)  X(:) Y(:) Z(:) One(:)];
    Coeff=M\T(:);
    P=zeros(3,3);
    A=zeros(3,1);

    P(1,1)=2*Coeff(1);
    P(2,2)=2*Coeff(2);
    P(3,3)=2*Coeff(3);
    P(1,2)=Coeff(4);P(2,1)=Coeff(4);
    P(2,3)=Coeff(5);P(3,2)=Coeff(5);
    P(3,1)=Coeff(6);P(1,3)=Coeff(6);
    A(:)=-Coeff(7:9);
    B=P\A;

    dec=dec-B;
end


end