function tiff2mat(fil2)
display(sprintf('Converting %s...',fil2))
tmp=imread(fil2);
if ~strcmp(class(tmp),'uint8')
    error('Only 8-bit encoding is accepted !!!')
end

sizeim=tiffdims(fil2);
jm3=zeros(sizeim,'uint8');


for jj=1:sizeim(3)
    tmp=imread(fil2,jj);
%    tmp=(permute(tmp,[2,1,3]));
%    tmp=fliplr(tmp);
    jm3(:,:,jj)=uint8(tmp);
end
[~,fil2, ext] = fileparts(fil2);
save(fil2,'jm3','sizeim','-v7.3');

end