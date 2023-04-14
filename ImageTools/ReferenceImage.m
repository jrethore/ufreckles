function ReferenceImage(nmod,ncam,roi)
%filref,nscale,roi,irestart
if nargin<1, nmod=1;end
if nargin<2, ncam=1;end
if nargin<3, roi=[];end

load(fullfile('TMP','params'),'param');
tic;
if isempty(roi)
    roi=param.roi;
end
if iscell(param.reference_image)
    filref=param.reference_image{ncam};
else
    filref=param.reference_image;
end
reverse=0;
if isfield(param,'reverse_image')
    reverse=param.reverse_image;
end
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=0;
end
param0=param;
if ~isfield(param0,'deformed_image')
    param.restart=0;
    save(fullfile('TMP','params'),'param');
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
nscale=param.nscale;
if length(roi)==4
    if ~isfield(param0,'deformed_image')
        reader=VideoReader(filref);
        im0=readim(reader,1);
    else
        im0=double(readim(filref));
    end
    %        im0=(im0-min(im0(:)))/(max(im0(:))-min(im0(:)));
    
    sizeim=size(im0);
    if length(sizeim)==3
        im0=mean(im0,3);
    end
    if reverse
        im0=im0';
    end
    sizeim=size(im0);
    save(fullfile('TMP',sprintf('sample%d',ncam-1)),'im0','sizeim','roi','-v7.3');
    if ~isfield(param0,'deformed_image')
        save(fullfile('TMP',sprintf('sample%d',ncam-1)),'reader','-append');
    end
    if psample
        [Yi,Xi]=meshgrid(roi(3):psample:roi(4),roi(1):psample:roi(2));
        [im00]=mexInterpSpline(Xi,Yi,im0);
        im0=im00;
    else
        im0=(im0(roi(1):roi(2),roi(3):roi(4)));
    end
    sizeim=size(im0);
    save(fullfile('TMP',sprintf('sample%d_%d',(ncam-1),0)),'im0','sizeim','roi','-v7.3');
    
    disp(sprintf('Loading reference image...%6.2f s',toc));
    
    for iscale=2:nscale
        tic();
        NestedCoarseImage();
        sizeim=size(im0);
        save(fullfile('TMP',sprintf('sample%d_%d',(ncam-1),iscale-1)),'im0','sizeim','roi','-v7.3');
        disp(sprintf('   Coarsening level %d...%6.2f s',iscale,toc()));
    end
elseif length(roi)==6
    delete(fullfile('TMP',sprintf('sample%d*',ncam-1)));
    delete(fullfile('TMP',sprintf('dsample%d*',ncam-1)));
    if isfield(param0,'stack_size')
        sizeim=param0.stack_size;
        fid=fopen(filref,'r');
        im0=fread(fid,prod(sizeim));
        im0=reshape(im0,sizeim);
        fclose(fid);
    else
 [~, ~, extref] = fileparts(filref);
 switch extref
     case '.mat'
         load(filref,'jm3');
        im0=(jm3);
        clear jm3
     case {'.tif','.tiff'}
          im0=readTIFFasRAW(filref);
 end
        sizeim=size(im0);
        
    end
    %save(fullfile('TMP',sprintf('sample%d',ncam-1)),'im0','sizeim','roi','-v7.3');
    save(fullfile('TMP',sprintf('sample%d',ncam-1)),'sizeim','roi','-v7.3');
    im0=(im0(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)));
    sizeim=size(im0);
    save(fullfile('TMP',sprintf('sample%d_%d',(ncam-1),0)),'sizeim','roi','-v7.3');
    switch param0.stack_format
        case 'bin'
            fid=fopen(fullfile('TMP',sprintf('dsample%d_%d',(ncam-1),0)),'w');
            fwrite(fid,im0);
            fclose(fid);
        case 'mat'
            save(fullfile('TMP',sprintf('sample%d_%d',(ncam-1),0)),'im0','-append');
    end
    
    disp(sprintf('Loading reference image...%6.2f s',toc));
    
    for iscale=2:nscale
        tic();
        NestedCoarseImage();
        sizeim=size(im0);
        save(fullfile('TMP',sprintf('sample%d_%d',(ncam-1),iscale-1)),'sizeim','roi','-v7.3');
        switch param0.stack_format
            case 'bin'
                fid=fopen(fullfile('TMP',sprintf('dsample%d_%d',(ncam-1),iscale-1)),'w');
                fwrite(fid,im0);
                fclose(fid);
            case 'mat'
                save(fullfile('TMP',sprintf('sample%d_%d',(ncam-1),iscale-1)),'im0','-append');
        end
        disp(sprintf('   Coarsening level %d...%6.2f s',iscale,toc()));
    end
    
end


    function NestedCoarseImage()
        
        scale=2;
        imsiz0=size(im0);
        imsiz1=floor(imsiz0/2);
        nn=2*imsiz1;
        ndim=length(nn);
        if ndim==2
            im0=im0(1:nn(1),1:nn(2));
        elseif ndim==3
            im0=im0(1:nn(1),1:nn(2),1:nn(3));
        end
        
        im0=reshape(im0,scale,prod(nn)/scale);
        im0=mean(im0,1);
        nn(1)=nn(1)/scale;
        im0=reshape(im0,nn);
        
        if ndim==2
            im0=im0';
            im0=reshape(im0,scale,prod(nn)/scale);
            im0=mean(im0,1);
            nn(2)=nn(2)/scale;
            im0=reshape(im0,nn([2,1]));
            im0=im0';
            %toc();
        elseif ndim==3
            im0=permute(im0,[2,3,1]);
            im0=reshape(im0,scale,prod(nn)/scale);
            im0=mean(im0,1);
            nn(2)=nn(2)/scale;
            im0=reshape(im0,nn([2,3,1]));
            im0=permute(im0,[3,1,2]);
            
            im0=permute(im0,[3,1,2]);
            im0=reshape(im0,scale,prod(nn)/scale);
            im0=mean(im0,1);
            nn(3)=nn(3)/scale;
            im0=reshape(im0,nn([3,1,2]));
            im0=permute(im0,[2,3,1]);
        end
        
    end




end

