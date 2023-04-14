function image = readTIFFasRAW(fname)
% 
% image = readTIFFasRAW(fname)
% 
% Read a TIFF file as a RAW file but directly reshape according the
% following shape : [dim_x, dim_y, dim_z].
% 
% It's handle unsigned int 8 bits and unsigned int 16 bits. The data has
% to be contigus and in big endian (the usual case from imagej).
%

image = [];

data = imfinfo(fname, 'TIFF');

shape = [data(1).Width, data(1).Height, length(data)];
offset = data(1).StripOffsets;
is_contigus = true;

if data(1).BitDepth == 8
    slice_size = data(1).Height*data(1).Width;
elseif data(1).BitDepth == 16
    slice_size = data(1).Height*data(1).Width*2;
end

for i = 2:length(data)
    if data(i).StripOffsets == offset + slice_size
        offset = offset + slice_size;
    else
        is_contigus = false;
        break
    end
end

fd = fopen(fname);

if is_contigus
    fseek(fd, data(1).StripOffsets, -1);
    if data(1).BitDepth == 8
        b = fread(fd, prod(shape), '*uint8');
    elseif data(1).BitDepth == 16
        % Warning : TIFF is in BIG ENDIAN
        b = fread(fd, prod(shape), '*uint16', 0, 'b');
    end
else
    disp TIFF file is not contigus
    disp TODO : be able to read that kind of file
    return
end

image = reshape(b, shape);
