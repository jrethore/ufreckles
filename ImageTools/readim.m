function im=readim(fil,iim)
[~,c,ext]=fileparts(fil);
if strcmp(ext,'.raw')
    sizeim=[str2num(c(end-8:end-5)),str2num(c(end-3:end))];
fid=fopen(fil,'rb');
im=fread(fid,'uint16');
    fclose(fid);
    im=reshape(im,sizeim);
else
if nargin>1
im=read(fil,iim);    
else
    im=imread(fil);
end

im=(permute(im,[2,1,3]));

for id=1:size(im,3)
    im(:,:,id)=fliplr(im(:,:,id));
end
end
%im=double(im);
end