function [U,V]=rbt(im0,im1)
subpix=1;
nsiz=size(im0);
nsiz2=2.^floor(log2(nsiz+.5));
off=floor((nsiz-nsiz2)/2);
im0=im0(1+off(1):nsiz2(1)+off(1),1+off(2):nsiz2(2)+off(2));
im1=im1(1+off(1):nsiz2(1)+off(1),1+off(2):nsiz2(2)+off(2));
ffim0=fft2(im0);ffim0(1,1,1)=0;
ffim1=fft2(im1);ffim1(1,1,1)=0;
ffinter=ffim0.*conj(ffim1);
inter=real(ifft2(ffinter));
shinter=fftshift(inter);
[mx,ind]=max(shinter(:));
[i1,i2]=ind2sub(nsiz2,ind);
ii=[i1,i2];
dec=-(ii-nsiz2/2-1)';
if (subpix==1)
    xi=min(max(i1+(-2:2),1),size(shinter,1));
    yi=min(max(i2+(-2:2),1),size(shinter,2));
    T=shinter(xi,yi);
    [Y,X]=meshgrid(yi-i2,xi-i1);
    One=ones(length(xi),length(yi));
    X2=X.^2;Y2=Y.^2;XY=X.*Y;
    M=[X2(:) Y2(:) XY(:) X(:) Y(:) One(:)];
    Coeff=M\T(:);
    P=zeros(2,2);
    A=zeros(2,1);

    P(1,1)=2*Coeff(1);
    P(2,2)=2*Coeff(2);
    P(1,2)=Coeff(3);P(2,1)=Coeff(3);
    A(:)=-Coeff(4:5);
    B=P\A;
if ~any(isnan(B))
    dec=dec-B;
else
    warning('Sub-pixel interpolation failed');
end
end
U=dec(1);
V=dec(2);

end