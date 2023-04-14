function  [im]=MCoarseImage(jm,scale)

        imsiz0=size(jm);
        imsiz1=floor(imsiz0/scale);
        nn=scale*imsiz1;

ndim=length(nn);
        if ndim==2
            im=jm(1:nn(1),1:nn(2));
        elseif ndim==3
            im=jm(1:nn(1),1:nn(2),1:nn(3));
        end

im=reshape(im,scale,prod(nn)/scale);
im=mean(im,1);
nn(1)=nn(1)/scale;
im=reshape(im,nn);

if ndim==2
im=im';
im=reshape(im,scale,prod(nn)/scale);
im=mean(im,1);
nn(2)=nn(2)/scale;
im=reshape(im,nn([2,1]));
im=im';
elseif ndim==3
  im=permute(im,[2,3,1]);
im=reshape(im,scale,prod(nn)/scale);
im=mean(im,1);
nn(2)=nn(2)/scale;
im=reshape(im,nn([2,3,1]));
im=permute(im,[3,1,2]);

im=permute(im,[3,1,2]);
im=reshape(im,scale,prod(nn)/scale);
im=mean(im,1);
nn(3)=nn(3)/scale;
im=reshape(im,nn([3,1,2]));
im=permute(im,[2,3,1]);
end

end

