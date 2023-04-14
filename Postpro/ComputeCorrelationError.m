function [mdisc,disc,Urbm]=ComputeCorrelationError(nmod,U,xo,yo,zo)
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
roi=param0.roi;
if iscell(param0.deformed_image)
    nim=length(param0.deformed_image);
else
    nim=1;
end

reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
psample=1;
if isfield(param0,'sampling_factor')
    psample=param0.sampling_factor;
end

iscale=1;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
disp(sprintf('Computing correlation error for scale %d...',iscale));
pscale=2^(iscale-1);
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim');
load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim');
load(fullfile('TMP',sprintf('%d_phiy_%d',nmod,10*(iscale-1))),'phiy');

if iscale==1
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
    if strcmp(param.basis,'fem')
        mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
        load(mesh_file,'ng');
        if ng>0||isfield(param,'gp_file')||isfield(param,'extrusion_parameters')
            load(fullfile('TMP','sample0'),'im0','sizeim');
            im00=im0;
            [im00]=mexInterpSpline((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
            im0=im00;
        end
    end
end
on=sum(phix,2)>0;


dynamic=max(im0(:))-min(im0(:));
mean0=mean(im0(:));
std0=std(im0(:));
im0=im0-mean0;
ndim=length(sizeim);

% load(fullfile('TMP',sprintf('%d_operator_%d',nmod,iscale-1)),'M');
% Mc=M;
% load(fullfile('TMP',sprintf('%d_phidf_%d',nmod,iscale-1)),'phidf');
% P=[ones(length(xo(:)),1),sparse(length(xo(:)),1),-yo(:);...
%    sparse(length(xo(:)),1),ones(length(xo(:)),1),xo(:);...
%    sparse(length(xo(:)),1),sparse(length(xo(:)),1),sparse(length(xo(:)),1)];
% Mc=P'*M*P;
% phidf=phidf*P;
% phixrbm=phix*P;
% phiyrbm=phiy*P;
% Urbm=zeros(3,nim);
Ux=phix*U;
%figure
%imagesc(reshape(Ux(:,1),sizeim))
Uy=phiy*U;
NestedResidual();
disc=mask*disc;
merror=0;
for iim=1:nim
    merrori=mean(abs(disc(on,iim)))/dynamic;
    merror=merror+merrori;
    mdisc(iim)=merrori;
end
merror=merror/nim;


disp(sprintf('Discrepancy wrt dyn. =%6.2f %% ',merror*100));

disp(sprintf('Enlapsed time for computation = %6.2f s',cputime-ttic));

fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
fwrite(fid,min(255,dynamic));
if dynamic>255
    disc=255*disc/dynamic;
end
fwrite(fid,abs(disc));
fclose(fid);

    function NestedResidual()
        disc=repmat(0,size(phix,1),nim);
        for jjm=1:nim
            Uxi=Ux(:,jjm);
            Uyi=Uy(:,jjm);
            
            if nim==1
                fildef=param0.deformed_image;
                
            else
                
                fildef=param0.deformed_image{jjm};
                im1=double(readim(fildef));
                if reverse
                    im1=im1';
                end
            end
         disci=mexInterpSpline((Xi-1)*psample+Uxi+roi(1),(Yi-1)*psample+Uyi+roi(3),im1);
        
        mean1=mean(disci(:));
        std1=std(disci(:));
        disci=disci-mean1;
        st=std0/std1;
        disci=(im0(:)-st*disci(:));
%  mean(abs(disci))
%         res=1;iter=1;conv=1.e-3;maxiter=20;
%         Urbmi=zeros(size(Mc,1),1);
%         while res>conv && iter<maxiter
%             F=phidf'*wdetJ*disci;
%             dU=Mc\F;
%             Urbmi=Urbmi+dU;
%             Uxi=Uxi+phixrbm*dU;
%             Uyi=Uyi+phiyrbm*dU;
%          disci=mexInterpSpline((Xi-1)*psample+Uxi+roi(1),(Yi-1)*psample+Uyi+roi(3),im1);
%                     disci=disci-mean1;
%         disci=(im0(:)-st*disci(:));
% res=norm(dU)/norm(Urbmi);
%  mean(abs(disci))
%         end
%         Urbm(:,jjm)=Urbmi;
        disc(:,jjm)=disci;
       end
    end
end

