function LoadLevelSets()
load(fullfile('TMP','params'),'param');
tic();
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=0;
end
if ~iscell(param.levelset_file)
    nbfis=1;
filref{1}=param.levelset_file;
else
    nbfis=numel(param.levelset_file);
    filref=param.levelset_file;
end
roi=param.roi;
fmax=20;
if isfield(param,'cz_length')
    fmax  =fmax+param.cz_length;
end
for ic=1:nbfis
fid=fullfile('TMP',sprintf('%d_levelsets',ic));
fid1=fullfile('TMP',sprintf('%d_levelsets_grads',ic));
fid2=fullfile('TMP',sprintf('%d_levelsets_cylco',ic));
load(filref{ic},'crack')
crack=(permute(crack,[2,1,3]));
for id=1:size(crack,3)
    crack(:,:,id)=fliplr(crack(:,:,id));
end
if psample
      [Yi,Xi]=meshgrid(roi(3):psample:roi(4),roi(1):psample:roi(2));
  crack=interp2(crack,Yi,Xi,'*linear');
else
crack=double(crack(roi(1):roi(2),roi(3):roi(4)));
end
sizeim=size(crack);
save(fid,'crack','sizeim','-v7.3');

nx=FDgradient(crack,1);
save(fid1,'nx');

ny=FDgradient(crack,2);
save(fid1,'ny','-append');
z=i*crack;
clear crack

nn=nx+i*ny;
clear nx ny

nnorm=abs(nn);
%nn=nn./nnorm;
clear nnorm
nn=nn*exp(-i*pi/2);
thetap=angle(nn);
theta=mean(thetap(:));
clear nn
save(fid2,'theta','thetap','-v7.3');

load(filref{ic},'front')
front=(permute(front,[2,1,3]));
for id=1:size(front,3)
    front(:,:,id)=fliplr(front(:,:,id));
end
if psample
  front=interp2(front,Yi,Xi,'*linear');
else
front=double(front(roi(1):roi(2),roi(3):roi(4)));
end
save(fid,'front','-append');

tx=FDgradient(front,1);
save(fid1,'tx','-append');
clear tx

ty=FDgradient(front,2);
save(fid1,'ty','-append');
clear ty
z=z+front;
clear front

dist=max(abs(z),1);
save(fid2,'dist','-append');
clear dist

angl=angle(z);
save(fid2,'angl','-append');
clear angl
onfaces=(imag(z)<=0)&(imag(z)>-1)&(real(z)<=fmax);
clear z
save(fid,'onfaces','-append');
end
disp(sprintf('Loading levelsets...%6.2f s',toc()));


end