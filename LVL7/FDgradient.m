function [grad]=FDgradient(im,dim)

ndim=length(size(im));
if ndim==2
    if dim==2
    im=im';
    end
    grad=mexFDGradient(im);
  if dim==2
    grad=grad';
    end

    
elseif  ndim==3
    if dim==2
        im=permute(im,[2,3,1]);
    elseif dim==3
     im=permute(im,[3,1,2]);
   
    end
    
   sizeim=size(im);
   grad = im;
    
         grad(1,:,:) = im(2,:,:)-im(1,:,:);
      grad(sizeim(1),:,:) = im(sizeim(1),:,:)-im(sizeim(1)-1,:,:);

        grad(2:sizeim(1)-1,:,:) = (im(3:sizeim(1),:,:) - im(1:sizeim(1)-2,:,:))/2;

    if dim==2
        grad=permute(grad,[3,1,2]);
    elseif dim==3
        grad=permute(grad,[2,3,1]);

    end
end
 
end
