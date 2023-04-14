function [im1]=CoarseImage(im0)


imsiz0=size(im0);
imsiz1=floor(imsiz0/2);
imsiz0=2*imsiz1;
if length(imsiz0)==2
im1=MCoarseImage(im0(1:imsiz0(1),1:imsiz0(2)),2);
elseif length(imsiz0)==3
im1=MCoarseImage(im0(1:imsiz0(1),1:imsiz0(2),1:imsiz0(3)),2);
end

end
