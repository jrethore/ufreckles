function CreateBasisFunction(iscale,nmod)
load(fullfile('TMP','params'));
onflight=0;
if isfield(param,'onflight')
    onflight=param.onflight;
end
if onflight
    global phiy phix Xi Yi wdetJ inde on
end
doms=0;
param0=param;
roi=param0.roi;
tic;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
pscale=2^(iscale-1);
sbasis=param.basis;
if isfield(param0,'reference_image')
    if iscell(param0.reference_image)
        ncams=length(param0.reference_image);
    else
        ncams=1;
    end
else
    ncams=1;
end
for ncam=1:ncams
    
    nmod2=nmod*10^(ncam-1);
    inde=[];
    switch sbasis
        case 'disto'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            wdetJ=1;
            siz=size(Xi);
            Zi=Xi(:)-mean(Xi(:))+i*(Yi(:)-mean(Yi(:)));
            ri=abs(Zi);
            angl=angle(Zi);
            Urbt=param.rbt;
            E11=repmat(0,numel(Yi),3);
            E12=repmat(0,numel(Yi),3);
            E22=repmat(0,numel(Yi),3);
            for in=1:3
                E11(:,in)=cos(angl).*cos(angl).*(in*ri.^(in-1));
                E12(:,in)=sin(angl).*cos(angl).*(in*ri.^(in-1));
                E22(:,in)=sin(angl).*sin(angl).*(in*ri.^(in-1));
            end
            phix=[1+0*ri,0*ri,E11*Urbt(1)+E12*Urbt(2)];
            phiy=[0*ri,1+0*ri,E12*Urbt(1)+E22*Urbt(2)];
        case 'rbt'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            wdetJ=1;
            if iscale>1&doms
                [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                siz=size(Xi);
                phix=zeros(prod(siz), 2);
                phix(:,1)=1;
                phiy=zeros(prod(siz), 2);
                phiy(:,2)=1;
                
                Xi=Xi(:);
                Yi=Yi(:);
                save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','-v7.3');
                save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
            end
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            siz=size(Xi);
            phix=zeros(prod(siz), 2);
            phix(:,1)=1;
            phiy=zeros(prod(siz), 2);
            phiy(:,2)=1;
        case 'rbm'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            wdetJ=1;
            if iscale>1&doms
                [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                siz=size(Xi);
                phix=zeros(prod(siz), 3);
                phix(:,1)=1;
                phix(:,3)=-Yi(:);
                phiy=zeros(prod(siz), 3);
                phiy(:,2)=1;
                phiy(:,3)=Xi(:);
                
                Xi=Xi(:);
                Yi=Yi(:);
                save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','-v7.3');
                save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
            end
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            siz=size(Xi);
            phix=zeros(prod(siz), 3);
            phix(:,1)=1;
            phix(:,3)=-Yi(:);
            phiy=zeros(prod(siz), 3);
            phiy(:,2)=1;
            phiy(:,3)=Xi(:);
        case 'cordin'
            load(fullfile('TMP',sprintf('sample%d',ncam-1)),'sizeim');
            sizeim0=sizeim;
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            if iscale>1&doms
                [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                Xi=Xi(:);
                Yi=Yi(:);
                wdetJ=1;
                X=(param0.roi(1)-1+Xi)/(sizeim0(1))-0.5;
                Y=(param0.roi(3)-1+Yi)/(sizeim0(2))-0.5;
                
                phix=zeros(prod(sizeim*pscale), 16);
                phiy=zeros(prod(sizeim*pscale), 16);
                
                phix(:,1)=1;
                phix(:,3)=X(:);
                phix(:,5)=Y(:);
                phix(:,7)=Y(:).*X(:);
                phix(:,9)=X(:).^2;
                phix(:,11)=Y(:).^2;
                phix(:,13)=(X(:).^2).*Y(:);
                phix(:,15)=(Y(:).^2).*X(:);
                
                phiy(:,1+1)=1;
                phiy(:,3+1)=X(:);
                phiy(:,5+1)=Y(:);
                phiy(:,7+1)=Y(:).*X(:);
                phiy(:,9+1)=X(:).^2;
                phiy(:,11+1)=Y(:).^2;
                phiy(:,13+1)=(X(:).^2).*Y(:);
                phiy(:,15+1)=(Y(:).^2).*X(:);

                save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','-v7.3');
                save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
            end
            [Yi Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
            Xi=Xi(:);
            Yi=Yi(:);
            wdetJ=1;
            X=(param0.roi(1)-1+0.+(Xi-0.)*pscale)/(sizeim0(1))-0.5;
            Y=(param0.roi(3)-1+0.+(Yi-0.)*pscale)/(sizeim0(2))-0.5;
            [phix,phiy]=GetCordinCorrectionFied(X,Y);
        case 'uni'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,1-1)),'sizeim');
            sizeim0=sizeim;
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            if iscale>1&doms
                [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                Xi=Xi(:);
                Yi=Yi(:);
                wdetJ=1;
                X=Xi/mean(sizeim0)-0.5;
                Y=Yi/mean(sizeim0)-0.5;
                
                phix=zeros(prod(sizeim*pscale), 6);
                phiy=zeros(prod(sizeim*pscale), 6);
                
                phix(:,1)=1;
                phiy(:,2)=1;
                phix(:,3)=Y(:);
                phiy(:,3)=-X(:);
                
                phix(:,4)=X(:);
                phiy(:,5)=Y(:);
                
                phix(:,6)=0.5*Y(:);
                phiy(:,6)=0.5*X(:);
                save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','-v7.3');
                save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
            end
            [Yi Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
            Xi=Xi(:);
            Yi=Yi(:);
            wdetJ=1;
            X=Xi*pscale/mean(sizeim0)-0.5;
            Y=Yi*pscale/mean(sizeim0)-0.5;
            
            phix=zeros(prod(sizeim), 6);
            phiy=zeros(prod(sizeim), 6);
            
            phix(:,1)=1;
            phiy(:,2)=1;
            phix(:,3)=Y(:);
            phiy(:,3)=-X(:);
            
            phix(:,4)=X(:);
            phiy(:,5)=Y(:);
            
            phix(:,6)=0.5*Y(:);
            phiy(:,6)=0.5*X(:);
            
            % %             x=(1:sizeim(1))-0.5-0.5*(sizeim(1)+1);
            % %             y=(1:sizeim(2))-0.5;
            %             x=(1:sizeim(1))-0.5;
            %             y=(1:sizeim(2))-0.5;
            %             lx=sizeim(1);
            %             ly=sizeim(2);
            %             [Y,X]=meshgrid(y/ly-0.5,x/lx-0.5);
            %             phix=zeros(prod(sizim), 6);
            %             phiy=zeros(prod(sizim), 6);
            %
            % %             phix(:,1)=1;
            % %             phiy(:,2)=1;
            %
            % %             phix(:,3)=(Y(:)-0.5*ly);
            % %             phiy(:,3)=-(X(:)-0.5*lx);
            % %
            % %             phix(:,4)=(X(:)-0.5*lx);
            % %             phiy(:,5)=(Y(:)-0.5*ly);
            % %
            % %             phix(:,6)=0.5*(Y(:)-0.5*ly);
            % %             phiy(:,6)=0.5*(X(:)-0.5*lx);
            % %             Xi=X+0.5+0.5*(sizeim(1)+1);
            % %             Yi=Y+0.5;
            %
            %             phix(:,3)=Y(:);
            %             phiy(:,3)=-X(:);
            %
            %             phix(:,4)=X(:);
            %             phiy(:,5)=Y(:);
            %
            %             phix(:,6)=0.5*Y(:);
            %             phiy(:,6)=0.5*X(:);
            %             Xi=(X+0.5)*lx;
            %             Yi=(Y+0.5)*ly;
            %             wdetJ=1;
        case 'affine'
            %             load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            %             [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            %             wdetJ=1;
            %             siz=size(Xi);
            %             phix=zeros(prod(siz), 4);
            %             phix(:,1)=1;
            %             phix(:,3)=Xi(:);
            %             phix(:,4)=-Yi(:);
            %             phiy=zeros(prod(siz), 4);
            %             phiy(:,2)=1;
            %             phiy(:,3)=Yi(:);
            %             phiy(:,4)=Xi(:);
            
            
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,1-1)),'sizeim');
            sizeim0=sizeim;
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            if iscale>1&doms
                [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                Xi=Xi(:);
                Yi=Yi(:);
                wdetJ=1;
                X=Xi/mean(sizeim0)-0.;
                Y=Yi/mean(sizeim0)-0.;
                
                phix=zeros(prod(sizeim*pscale), 4);
                phiy=zeros(prod(sizeim*pscale), 4);
                
                phix(:,1)=1;
                phiy(:,2)=1;
                phix(:,4)=Y(:);
                phiy(:,4)=-X(:);
                
                phix(:,3)=X(:);
                phiy(:,3)=Y(:);
                save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','-v7.3');
                save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
            end
            [Yi Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
            Xi=Xi(:);
            Yi=Yi(:);
            wdetJ=1;
            X=Xi*pscale/mean(sizeim0)-0.;
            Y=Yi*pscale/mean(sizeim0)-0.;
            
            phix=zeros(prod(sizeim), 4);
            phiy=zeros(prod(sizeim), 4);
            
            phix(:,1)=1;
            phiy(:,2)=1;
            phix(:,4)=Y(:);
            phiy(:,4)=-X(:);
            
            phix(:,3)=X(:);
            phiy(:,3)=Y(:);
            
            
        case 'beam'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            
            x=(1:sizeim(1))-0.5-0.5*(sizeim(1)+1);
            y=(1:sizeim(2))-0.5;
            [Y,X]=meshgrid(y,x);
            siz=size(Y);
            l=sizeim(2);
            phix=zeros(prod(siz), 6);
            phiy=zeros(prod(siz), 6);
            
            phiy(:,1)=(1-Y(:)/l);
            phiy(:,2)=Y(:)/l;
            
            tmp=(2*Y.^3-3*l*Y.^2+l^3)/l^3;
            phix(:,3)=tmp(:);
            tmp=(6*Y.^2-6*l*Y)/l^3;
            phiy(:,3)=-X(:).*tmp(:);
            
            tmp=(Y.^3-2*l*Y.^2+Y*l^2)/l^3;
            phix(:,4)=tmp(:);
            tmp=(3*Y.^2-4*l*Y+l^2)/l^3;
            phiy(:,4)=-X(:).*tmp(:);
            
            tmp=(-2*Y.^3+3*l*Y.^2)/l^3;
            phix(:,5)=tmp(:);
            tmp=(-6*Y.^2+6*l*Y)/l^3;
            phiy(:,5)=-X(:).*tmp(:);
            
            tmp=(Y.^3-l*Y.^2)/l^3;
            phix(:,6)=tmp(:);
            tmp=(3*Y.^2-2*l*Y)/l^3;
            phiy(:,6)=-X(:).*tmp(:);
            Xi=X+0.5+0.5*(sizeim(1)+1);
            Yi=Y+0.5;
            wdetJ=1;
        case 'vic-beam'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            p=param.degree;
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod2,iscale-1));
            if isfield(param,'beam_type')
                beamtype=strcmp(param.beam_type,'timoshenko');
            else
                beamtype=0;
            end
            type_nurbs=~strcmp(param.continuity,'c0');
            [phix,phiy]=CreateNURBSBasis1D(mesh_file,sizeim,p,pscale,beamtype,type_nurbs);
            phiy=0*phiy(:,3:size(phiy,2));
            phix=phix(:,3:size(phix,2));
            wdetJ=1;
        case 'nurbs-beam'
            warning('COMPATIBILITY WITH NEW NURBS STUFF NOT CHECKED !!!');
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            p=max(param.degree);
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod2,iscale-1));
            if isfield(param,'beam_type')
                beamtype=strcmp(param.beam_type,'timoshenko');
            else
                beamtype=0;
            end
            type_nurbs=~strcmp(param.continuity,'c0');
            [phix,phiy]=CreateNURBSBasis1D(mesh_file,sizeim,p,pscale,beamtype,type_nurbs);
            wdetJ=1;
            
        case 'KM'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            wdetJ=1;
            ind=param.km_indices;
            if isfield(param,'crack_id')
                ic=param.crack_id;
            else
                ic=1;
            end
            
            if isfield(param,'mask_radius')
                rmin=param.mask_radius;
            else
                rmin=0;
            end
            if isfield(param,'zoi_radius')
                rmax=param.zoi_radius;
            else
                rmax=0;
            end
            if isfield(param,'mask_width')
                dmin=param.mask_width;
            else
                dmin=0;
            end
            kappa=param.kolosov;
            modes=param.modes;
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
            sizemask=numel(mask);
            if rmin>0
                load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack','front');
                load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'dist');
                mask1=double(~((dist<rmin)|((abs(crack)<dmin)&(front<0))));
                clear crack front dist
                mask=mask*diag(sparse(mask1(:)));
                clear mask1
                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask','-v7.3');
            end
            if rmax>0
                load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack','front');
                load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'dist');
                mask1=double(dist<rmax);
                clear crack front dist
                mask=mask*diag(sparse(mask1(:)));
                clear mask1
                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask','-v7.3');
            end
            if sizemask==1
                ipix=[];
            else
                ipix=find(diag(mask)>0);
            end
            phic=CreateKMBasis(modes,ind,kappa,ipix,ic);
            phix=real(phic);
            phiy=imag(phic);
            clear phic
        case 'KM+'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            wdetJ=1;
            ind=[0,2];
            if isfield(param,'crack_id')
                ic=param.crack_id;
            else
                ic=1;
            end
            kappa=param.kolosov;
            modes=param.modes;
            if isfield(param,'mask_radius')
                rmin=param.mask_radius;
            else
                rmin=0;
            end
            if isfield(param,'mask_width')
                dmin=param.mask_width;
            else
                dmin=0;
            end
            lmin=param.cz_length;
            %            dz=param.tip_step;
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
            if rmin>0
                load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack','front');
                load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'dist');
                mask1=double(~((dist<rmin)|((abs(crack)<dmin)&(front<0))));
                clear crack front dist
                mask=mask*diag(sparse(mask1(:)));
                clear mask1
                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask','-v7.3');
            end
            sizemask=numel(mask);
            if sizemask==1
                ipix=[];
            else
                ipix=find(diag(mask)>0);
            end
            load(fullfile('TMP',sprintf('%d_iphi_%d',nmod,iscale-1)),'iphi','ztip','dz');
            
            phix=[];
            phiy=[];
            if any(modes==1)
                for iz=1:length(ztip)
                    phic=CreateKMBasis(1,1,kappa,ipix,ic,ztip(iz));
                    phix=[phix,real(phic)*dz(iz)];
                    phiy=[phiy,imag(phic)*dz(iz)];
                    
                    %                     figure
                    %                     imagesc(reshape(imag(phic(:,1)),sizeim))
                    %                     colorbar
                    
                end
            end
            if any(modes==2)
                for iz=1:length(ztip)
                    phic=CreateKMBasis(2,1,kappa,ipix,ic,ztip(iz));
                    phix=[phix,real(phic)*dz(iz)];
                    phiy=[phiy,imag(phic)*dz(iz)];
                end
            end
            
            if length(modes)==2
                iphi=blkdiag(iphi,iphi);
            end
            phix=phix*iphi;
            phiy=phiy*iphi;
            %             tmpx=0;
            %             tmpy=0;
            %             for iz=1:length(ztip)
            %                 phic=CreateKMBasis(1,1,kappa,ipix,-ztip(iz));
            %                 tmpx=tmpx+real(phic);
            %                 tmpy=tmpy+imag(phic);
            %
            %             end
            %             phix=tmpx;
            %             phiy=tmpy;
            %             tmpx=0;
            %             tmpy=0;
            %             for iz=1:length(ztip)
            %                 phic=CreateKMBasis(2,1,kappa,ipix,-ztip(iz));
            %                 tmpx=tmpx+real(phic);
            %                 tmpy=tmpy+imag(phic);
            %
            %             end
            %             phix=[phix,tmpx];
            %             phiy=[phiy,tmpy];
            
            
            
            phic=CreateKMBasis([1,2],ind,kappa,ipix,ic);
            phix=[real(phic),phix];
            phiy=[imag(phic),phiy];
            
            %             if any(modes==2)
            %                     phic=CreateKMBasis(2,1,kappa,ipix,ic);
            %                     phix=[phix,real(phic)];
            %                     phiy=[phiy,imag(phic)];
            %             end
            
            if isfield(param,'sub_indices')
                sind=param.sub_indices;
                if ~isempty(sind)
                    phic=CreateKMBasis(modes,sind,kappa,ipix,ic,lmin);
                    phix=[phix,real(phic)];
                    phiy=[phiy,imag(phic)];
                end
            end
            
            
            clear phic
        case 'CZ'
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            a=norm(sizeim);
            wdetJ=1;
            ind=[0,2];
            if isfield(param,'crack_id')
                ic=param.crack_id;
            else
                ic=1;
            end
            kappa=param.kolosov;
            modes=param.modes;
            if isfield(param,'mask_radius')
                rmin=param.mask_radius;
            else
                rmin=0;
            end
            if isfield(param,'mask_width')
                dmin=param.mask_width;
            else
                dmin=0;
            end
            load(fullfile('TMP',sprintf('%d_iphi_%d',nmod,iscale-1)),'iphi','ztip','dz','zn');
            lmin=max(zn);
            %            dz=param.tip_step;
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
            if rmin>0
                load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack','front');
                load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'dist');
                mask1=double(~((dist<rmin)|((abs(crack)<dmin)&(front<0))));
                clear crack front dist
                mask=mask*diag(sparse(mask1(:)));
                clear mask1
                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask','-v7.3');
            end
            sizemask=numel(mask);
            if sizemask==1
                ipix=[];
            else
                ipix=find(diag(mask)>0);
            end
            
            phix=[];
            phiy=[];
            if any(modes==1)
                phic=CreateKMBasis(1,1,kappa,ipix,ic,lmin);
                for iz=1:length(ztip)
                    phix=[phix,real(phic)*dz(iz)/sqrt(1-(ztip(iz)+a-lmin)/a)];
                    phiy=[phiy,imag(phic)*dz(iz)/sqrt(1-(ztip(iz)+a-lmin)/a)];
                    
                end
            end
            if any(modes==2)
                phic=CreateKMBasis(2,1,kappa,ipix,ic,lmin);
                for iz=1:length(ztip)
                    phix=[phix,real(phic)*dz(iz)/sqrt(1-(ztip(iz)+a-lmin)/a)];
                    phiy=[phiy,imag(phic)*dz(iz)/sqrt(1-(ztip(iz)+a-lmin)/a)];
                end
            end
            
            if length(modes)==2
                iphi=blkdiag(iphi,iphi);
            end
            phix=phix*iphi;
            phiy=phiy*iphi;
            
            phic=CreateKMBasis([1,2],ind,kappa,ipix,ic);
            phix=[real(phic),phix];
            phiy=[imag(phic),phiy];
            
            
            if isfield(param,'sub_indices')
                sind=param.sub_indices;
                if ~isempty(sind)
                    phic=CreateKMBasis(modes,sind,kappa,ipix,ic,lmin);
                    phix=[phix,real(phic)];
                    phiy=[phiy,imag(phic)];
                end
            end
            
            
            clear phic
            
        case {'fem','sfem'}
            
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod2,iscale-1));
            load(mesh_file,'xo','yo','ng')
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
            if iscale>1&doms
                phi1=CreateFiniteElementBasis(mesh_file,sizeim,pscale,unmasked_nodes);
                phi0=sparse(size(phi1,1),size(phi1,2));
                wdetJ=GetWeigthDetJ(mesh_file,sizeim,pscale);
                if ng==0
                    [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                else
                    Xi=phi1*((xo-1)*pscale+1);
                    Yi=phi1*((yo-1)*pscale+1);
                end
                Xi=Xi(:);
                Yi=Yi(:);
                phix=[phi1,phi0];
                on=sum(abs(phi1),2)>0;
                save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','on','-v7.3');
                phiy=[phi0,phi1];
                save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
            end
            if ~isfield(param,'gp_file')||iscale>1
                [wdetJ,inde]=GetWeigthDetJ(mesh_file,sizeim);
                phi1=CreateFiniteElementBasis(mesh_file,sizeim,1,unmasked_nodes);
                if ng==0
                    [Yi Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
                else
                    Xi=phi1*xo(:);
                    Yi=phi1*yo(:);
                end
            else
                toto=dlmread(param.gp_file);
                Xi=toto(:,1);
                Yi=toto(:,2);
                figure
                plot(Xi,Yi,'x')
                if isfield(param,'gluing_parameters')
                    load(fullfile('TMP','sample0'),'im0');
                    figure
                    
                    colormap(gray)
                    imagesc(im0')
                    hold on;
                    rect1=rectangle;
                    set(rect1,'position',[roi(1),roi(3),roi(2)-roi(1),roi(4)-roi(3)],...
                        'EdgeColor','yellow',...
                        'LineStyle','--',...
                        'LineWidth',2)
                    plot(Xi,Yi,'b+','LineWidth',0.5,'MarkerSize',10);
                    axis xy
                    axis image
                    [Xi,Yi]=GlueMesh(nmod,Xi,Yi);
                    plot(Xi,Yi,'rx','LineWidth',0.5,'MarkerSize',10);
                    Xi=Xi-roi(1)+1;
                    Yi=Yi-roi(3)+1;
                    
                end
                wdetJ=diag(sparse(toto(:,4)));
                toto=dlmread(param.basis_function_file);
                indi=toto(:,1);
                indj=toto(:,2);
                val=toto(:,3);
                phi1=sparse(indi,indj,val,max(indi),max(indj));
            end
            if ~onflight
            phi=phi1;
            save(fullfile('TMP',sprintf('%d_phi_%d',nmod2,(iscale-1))),'phi','-v7.3');
            end
            phi0=sparse(size(phi1,1),size(phi1,2));
            phix=[phi1,phi0];
            phiy=[phi0,phi1];
        case 'nurbs'
            if iscale==1
                p=param.degree;
                type_nurbs=param.continuity;
                load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
                mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod2,iscale-1));
                wdetJ=1;
                [phi1,Xi,Yi]=CreateNURBSBasis(mesh_file,p,'pixels');
                phi0=sparse(size(phi1,1),size(phi1,2));
                phix=[phi1,phi0];
                phiy=[phi0,phi1];
                
            else
                
                load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
                mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod2,iscale-1));
                load(mesh_file,'xo','yo','ng')
                load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
                if doms
                    phi1=CreateFiniteElementBasis(mesh_file,sizeim,pscale,unmasked_nodes);
                    phi0=sparse(size(phi1,1),size(phi1,2));
                    wdetJ=GetWeigthDetJ(mesh_file,sizeim,pscale);
                    [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                    Xi=Xi(:);
                    Yi=Yi(:);
                    phix=[phi1,phi0];
                    on=sum(abs(phi1),2)>0;
                    save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','on','-v7.3');
                    phiy=[phi0,phi1];
                    save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
                end
                [wdetJ,inde]=GetWeigthDetJ(mesh_file,sizeim);
                phi1=CreateFiniteElementBasis(mesh_file,sizeim,1,unmasked_nodes);
                [Yi Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
                
                Xi=Xi(:);
                Yi=Yi(:);
                phi0=sparse(size(phi1,1),size(phi1,2));
                phix=[phi1,phi0];
                phiy=[phi0,phi1];
                
                
                
                
            end
            
        case 'btri'
            if iscale==1
                load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
                mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod2,iscale-1));
                wdetJ=1;
                [phi1,Xi,Yi,wdetJ]=CreateBezierTriangleBasis(mesh_file,'pixels');
                phi0=sparse(size(phi1,1),size(phi1,2));
                phix=[phi1,phi0];
                phiy=[phi0,phi1];
                
            else
                
                load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'sizeim');
                mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod2,iscale-1));
                load(mesh_file,'xo','yo','ng')
                load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
                if doms
                    phi1=CreateFiniteElementBasis(mesh_file,sizeim,pscale,unmasked_nodes);
                    phi0=sparse(size(phi1,1),size(phi1,2));
                    wdetJ=GetWeigthDetJ(mesh_file,sizeim,pscale);
                    [Yi Xi]=meshgrid(1:(sizeim(2)*pscale),1:(sizeim(1)*pscale));
                    Xi=Xi(:);
                    Yi=Yi(:);
                    phix=[phi1,phi0];
                    on=sum(abs(phi1),2)>0;
                    save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,10*(iscale-1))),'phix','Xi','Yi','wdetJ','sizeim','on','-v7.3');
                    phiy=[phi0,phi1];
                    save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,10*(iscale-1))),'phiy','sizeim','-v7.3');
                end
                [wdetJ,inde]=GetWeigthDetJ(mesh_file,sizeim);
                phi1=CreateFiniteElementBasis(mesh_file,sizeim,1,unmasked_nodes);
                [Yi Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
                
                Xi=Xi(:);
                Yi=Yi(:);
                phi0=sparse(size(phi1,1),size(phi1,2));
                phix=[phi1,phi0];
                phiy=[phi0,phi1];
                
                
                
                
            end
        otherwise
            error ('INVALID FUNCTIONAL BASIS');
            
    end
    
    
    if (iscale==1)&&(isfield(param,'enrichment'))
        enr=param.enrichment;
        h=1;
        if isfield(param,'mesh_size_ratio')
            h=param.mesh_size_ratio;
        end
        if isfield(param,'enrichment_type')
            strat=param.enrichment_type;
        else
            strat='node';
        end
        mesh_file=fullfile('TMP',sprintf('%d_meshx_%d',nmod2,iscale-1));
        if ~iscell(param0.levelset_file)
            nbfis=1;
        else
            nbfis=numel(param0.levelset_file);
        end
        for ic=1:nbfis
            load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'enriched_pixels','tip_nodes','face_nodes','enriched_nodes','face_elts');
            %         load('TMP/maskn_0','maskn');
            %         if numel(maskn)>1
            %             maskp=diag(sparse(phi1*maskn(:)));
            %         else
            %             maskp=maskn;
            %         end
            %         save('TMP/maskp_0','maskp');
            
            do_disc_enrich=true;
            if ~(strcmp(enr,'crack_disc'))
                ind=param.km_indices;
                kappa=param.kolosov;
                modes=param.modes;
                km=CreateKMBasis(modes,ind,kappa,enriched_pixels);
                
                if isempty(face_nodes)
                    do_disc_enrich=false;
                end
                
                if isempty(enriched_pixels)
                    if strcmp(strat,'zone')
                        phikm=1;
                    elseif strcmp(strat,'node')
                        if h==1
                            phikm=phi1;
                        else
                            phikm=CreateFiniteElementBasis(mesh_file,sizeim);
                        end
                    end
                else
                    phikm=CreateFiniteElementBasis(mesh_file,sizeim,pscale,tip_nodes);
                    if strcmp(strat,'zone')
                        phikm=sum(phikm,2);
                    end
                    %                 if strcmp(strat,'zone')
                    %                     load(fullfile('TMP',sprintf('%d_cutoff_%d',nmod,iscale-1)),'cutoff');
                    %
                    %                     phikm=cutoff;
                    %                     clear cutoff
                    %                 else
                    %                     phikm=CreateFiniteElementBasis(mesh_file,sizeim,pscale,tip_nodes);
                    %
                    %                 end
                end
                % figure
                % imagesc(reshape((abs(phikm)>0)&(~(abs(km(:,1))>0)),sizeim));colorbar;
                % figure
                % imagesc(reshape(phikm,sizeim));colorbar;
                % figure
                % imagesc(reshape(real(km(:,1)),sizeim));colorbar;
                % figure
                % imagesc(reshape(diag(sparse(phikm))*real(km(:,1)),sizeim));colorbar;
                
                if strcmp(strat,'zone')
                    phix=[phix,diag(sparse(phikm))*real(km)];
                    phiy=[phiy,diag(sparse(phikm))*imag(km)];
                elseif strcmp(strat,'node')
                    for kk=1:size(km,2)
                        phix= [phix,diag(sparse(real(km(:,kk))))*phikm];
                        phiy= [phiy,diag(sparse(imag(km(:,kk))))*phikm];
                    end
                end
                clear km phikm
                
            end
            
            if do_disc_enrich
                if strcmp(strat,'node')
                    phih=CreateFiniteElementBasis(mesh_file,sizeim,pscale,face_nodes);
                    phih0=sparse(size(phih,1),size(phih,2));
                    load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack');
                    if isfield(param,'shift_enrichment')
                        opti_shift=param.shift_enrichment;
                    else
                        opti_shift=0;
                    end
                    if ~opti_shift
                        crack=crack(enriched_pixels);
                        heaviside=zeros(prod(sizeim),1);
                        heaviside(enriched_pixels)=2*double(crack>=0)-1;
                        clear crack
                        heaviside=diag(sparse(heaviside));
                        phih=(heaviside*phih);
                        
                    else
                        
                        disp('WARNING SHIFTED ENRICHMENT !!!');
                        load(mesh_file,'xo','yo');
                        hn=interp2(crack,yo,xo,'*linear');
                        
                        hn=double(hn>=0);
                        hn=diag(sparse(hn(face_nodes)));
                        crack=crack(enriched_pixels);
                        heaviside=zeros(prod(sizeim),1);
                        heaviside(enriched_pixels)=double(crack>=0);
                        clear crack
                        heaviside=diag(sparse(heaviside));
                        phih=(heaviside*phih-phih*hn);
                    end
                    
                    phix= [phix,phih,phih0];
                    phiy= [phiy,phih0,phih];
                    clear phih phih0 heaviside
                elseif strcmp(strat,'element')
                    
                    phih=CreateFiniteElementBasis(mesh_file,sizeim,pscale,enriched_nodes);
                    phie=CreateElementwiseBasis(mesh_file,sizeim,pscale,face_elts);
                    load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack');
                    load(mesh_file,'xo','yo');
                    hn=interp2(crack,yo,xo,'*linear');
                    
                    hn=double(hn>=0);
                    hn=hn(enriched_nodes(:));
                    crack=crack(enriched_pixels);
                    heaviside=zeros(prod(sizeim),1);
                    heaviside(enriched_pixels)=double(crack>=0);
                    clear crack
                    heaviside=(heaviside-phih*hn);
                    heaviside=diag(sparse(heaviside));
                    
                    phie=heaviside*phie;
                    phie0=sparse(size(phie,1),size(phie,2));
                    phix= [phix,phie,phie0];
                    phiy= [phiy,phie0,phie];
                    clear phie phie0 heaviside
                end
                
                
                
                
                
            end
        end
        
    end
    Xi=Xi(:);
    Yi=Yi(:);
    
    on=sum(abs(phix+1i*phiy),2)>0;
    Nddl_tot=size(phix,2);
if ~onflight
    save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,iscale-1)),'phiy','sizeim','Nddl_tot','-v7.3');
    save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,iscale-1)),'phix','Xi','Yi','wdetJ','inde','sizeim','Nddl_tot','on','-v7.3');
else
    save(fullfile('TMP',sprintf('%d_phiy_%d',nmod2,iscale-1)),'sizeim','Nddl_tot','-v7.3');
    save(fullfile('TMP',sprintf('%d_phix_%d',nmod2,iscale-1)),'sizeim','Nddl_tot','-v7.3');    
end
% if( (iscale==1) && (strcmp(sbasis,'nurbs')) && ~strcmp(type_nurbs,'c0'))
    %     save(fullfile('TMP',sprintf('%d_dphiy_%d',nmod2,iscale-1)),'dphiy');
    %     save(fullfile('TMP',sprintf('%d_dphix_%d',nmod2,iscale-1)),'dphix');
    % end
end
disp(sprintf('Creating basis function for model %d...%6.2f s',nmod,toc));

end

