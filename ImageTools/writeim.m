function im=writeim(fil,im)
for id=1:size(im,3)
    im(:,:,id)=fliplr(im(:,:,id));
end

im=(permute(im,[2,1,3]));
imwrite(im,fil);
end