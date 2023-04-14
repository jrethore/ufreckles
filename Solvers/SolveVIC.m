function [U]=SolveVIC(Uini,iscale,nmod)
check=0;
ttic=cputime;
load(fullfile('TMP',sprintf('sample0_%d',1-1)),'sizeim');
sizeim0=sizeim;
load(fullfile('TMP','params'),'param');
param0=param;
roi=param0.roi;
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
if iscell(param0.deformed_image)
    nim=length(param0.deformed_image);
else
    nim=1;
end
if isfield(param0,'restart')
    restart=param0.restart;
else
    restart=1;
end
maxiter=param0.iter_max;
conv=param0.convergance_limit;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
disp(sprintf('Starting resolution for scale %d...',iscale));

pscale=2^(1-1);
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim','mean0','std0','nband');
dynamic=max(im0(:))-min(im0(:));
ndim=length(sizeim);

if check
    load(fullfile('TMP',sprintf('sample0_%d',1-1)),'lso','ls1','nx','ny');
    load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'on');
load(fullfile('TMP',sprintf('%d_phi_%d',nmod,(iscale-1))),'phii','ls1i');
[xon,yon]=ind2sub(sizeim,nband(on));
nxon=-nx(nband(on));
nyon=-ny(nband(on));
end

load(fullfile('TMP',sprintf('%d_phix_%d',nmod,(iscale-1))),'phix','wdetJ','Xi','Yi');
load(fullfile('TMP',sprintf('%d_phiy_%d',nmod,(iscale-1))),'phiy');
load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
load(fullfile('TMP',sprintf('%d_phidf_%d',nmod,iscale-1)),'phidf');
load(fullfile('TMP',sprintf('%d_operator_%d',nmod,iscale-1)),'M','R');
%R=0*R;
M=M+R;
im0=im0-mean0;
im1=0;
U=Uini;
xi=roi(1)-1+(1:sizeim(1));
yi=roi(3)-1+(1:sizeim(2));
fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
fwrite(fid,min(255,dynamic));

if restart
    nmax=nim;
else
    if iscale==1
        nmax=nim;
    else
        if nim<2
            nmax=nim;
        else
            nmax=2;
        end
    end
end
%nmax=nim;

for ijm=1:nmax

display(sprintf('Image %d/%d',ijm,nmax));
res=1;
    ii=1;
    if restart
        Ui=Uini(:,ijm);
    else
                if ijm<3
                    Ui=Uini(:,ijm);
                else
                            DU=U(:,ijm-1)-U(:,ijm-2);
        
                Ui=U(:,ijm-1)+DU;
                end
%         if iscale==param.nscale
%             Ui=U(:,max(1,ijm-1));
%         else
%             Ui=Uini(:,ijm);
%         end
    end
    if check
    figure
    title(sprintf('Image %d/%d',ijm,nmax));
    hold on
    end
    while ( res>conv && ii< maxiter)
        Ux=phix*Ui;
        Uy=phiy*Ui;
        NestedResidual(ijm);
        disc=mask*disc;
        merror=std(disc(nband))/dynamic;

        F=phidf'*(wdetJ*disc)-R*Ui;
        if isfield(param0,'solver')
            solv=param0.solver;

            if strcmp(solv,'lu')
                if ii==1
                    [LT,UT]=lu(sparse(M));
                end
                dU=UT\(LT\F);
            elseif strcmp(solv,'chol')
                if ii==1
                    R=chol(sparse(M));
                end
                dU=R\(R'\F);
            end
        else
            dU=(M\F);
        end
        disp(sprintf('At iteration # %d',ii));
        disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
%        if ii==1
            nU0=norm(Ui+dU);
%        else
            res=norm(dU)/nU0;
            disp(sprintf('|dU|=%f',res));

%        end
        Ui=Ui+dU;
        if check
if exist('h')&&ii>1
    delete(h)
else
    h=plot(xon+roi(1)-1,yon+roi(3)-1,'r.','LineWidth',2);
%   [C,h]=contour(xi,yi,lso',[0 0],'Color','red','LineWidth',2);  
       print ('-djpeg', fullfile('FIG',sprintf('%s-scale-%d-iter-%d.jpg',param0.result_file,iscale,ii-1)));
%       print ('-depsc2', fullfile('FIG',sprintf('%s-scale-%d-iter-%d.eps',param0.result_file,iscale,ii-1)));
%pause
 delete(h)
 pause(0.1)
end
% uv=(phii*Ui);
% uvi=interp1(ls1i,uv,ls1(:),'linear','extrap');
%  cont=lso+reshape(uvi,sizeim);
%  [C,h]=contour(xi,yi,cont',[0 0],'Color','blue','LineWidth',2);
    h=plot(xon+roi(1)-1+nxon.*(phii*Ui),yon+roi(3)-1+nyon.*(phii*Ui),'r.','LineWidth',2);
        print ('-djpeg', fullfile('FIG',sprintf('%s-scale-%d-iter-%d.jpg',param0.result_file,iscale,ii)));
 %       print ('-depsc2', fullfile('FIG',sprintf('%s-scale-%d-iter-%d.eps',param0.result_file,iscale,ii)));
 pause(0.1)
        end
 %         figure
%         imagesc(reshape(phi*Ui,sizeim))
%         colorbar
%         pause
        ii=ii+1;
    end
    if check 
        close all
 pause(0.1)
    end
    U(:,ijm)=Ui;
    if dynamic>255
        disc=255*disc/dynamic;
    end
    fwrite(fid,abs(disc));
end
disp(sprintf('Enlapsed time for resolution = %6.2f s',cputime -ttic));
fclose(fid);
    function NestedResidual(kkk)
        if ii==1
            if nim==1
                fildef=param0.deformed_image;
                im1=double(readim(fildef));
            else
                fildef=param0.deformed_image{kkk};
                im1=double(readim(fildef));
             end
            if length(size(im1))==3
                im1=mean(im1,3);
%                im1=0.299 * im1(:,:,1)  + 0.587 *   im1(:,:,2)+ 0.114 *  im1(:,:,3);
            end
                if reverse
                    im1=im1';
                end
        end
%         figure
%         imagesc(im1)
        %         if iscale==2
%                  figure
%                  imagesc(reshape(Ux,sizeim*pscale))
%                  pause
        %         end
%        disc=mexInterpSpline((Xi+Ux+roi(1)-1),(Yi+Uy+roi(3)-1),im1);
                    sizeim1=size(im1);
if check
    if ii==1
                    h2=imagesc(im1');
                    colormap(gray)    
axis off;
axis xy;
axis image;
%delete(h2)
    end
end
        disc=(mexInterpSpline(min(sizeim1(1)-2,max(3,(Xi+Ux+roi(1)-1))),min(sizeim1(2)-2,max(3,(Yi+Uy+roi(3)-1))),im1));

        %                  figure
%                  imagesc(reshape(disc,sizeim*pscale))

%         if iscale>1
%             disc=reshape(disc,sizeim*pscale);
%             NestedCoarseImage(pscale);
%         end

        mean1=mean(disc(nband));
        std1=std(disc(nband));
%std1=max(disc(nband))-min(disc(nband));
%std1=std(disc(:));
%std1=sqrt(std(disc(nband))^2*numel(nband)/prod(sizeim));
        disc=disc-mean1;
%std1=max(im1(:))-min(im1(:));
        st=std0/std1;
        disc=(im0(:)-st*disc(:));


        function NestedCoarseImage(scale)

            imsiz0=size(disc);
            imsiz1=floor(imsiz0/2);
            nn=2*imsiz1;
            ndim=length(nn);
            if ndim==2
                disc=disc(1:nn(1),1:nn(2));
            elseif ndim==3
                disc=disc(1:nn(1),1:nn(2),1:nn(3));
            end

            disc=reshape(disc,scale,prod(nn)/scale);
            disc=mean(disc,1);
            nn(1)=nn(1)/scale;
            disc=reshape(disc,nn);

            if ndim==2
                disc=disc';
                disc=reshape(disc,scale,prod(nn)/scale);
                disc=mean(disc,1);
                nn(2)=nn(2)/scale;
                disc=reshape(disc,nn([2,1]));
                disc=disc';
                %toc();
            elseif ndim==3
                disc=permute(disc,[2,3,1]);
                disc=reshape(disc,scale,prod(nn)/scale);
                disc=mean(disc,1);
                nn(2)=nn(2)/scale;
                disc=reshape(disc,nn([2,3,1]));
                disc=permute(disc,[3,1,2]);

                disc=permute(disc,[3,1,2]);
                disc=reshape(disc,scale,prod(nn)/scale);
                disc=mean(disc,1);
                nn(3)=nn(3)/scale;
                disc=reshape(disc,nn([3,1,2]));
                disc=permute(disc,[2,3,1]);
            end

        end
    end
end
