function [Xf]=ProjectDisplacement(Xg,iscale,nmod,im1,icamd)

load(fullfile('TMP','params'));
param0=param;
thermo=(param0.thermo==1);
load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nnodes');

if length(Xg)==2
    init='bilin';
    if isfield(param0,'initialization');
        init=param0.initialization;
    end
    if isfield(param0,'sampling_factor')
        psample=param0.sampling_factor;
    else
        psample=1;
    end
    if isempty(im1)
        init='rbt';
    end
    switch init
        %     Nddl=prod(Nnodes);
        %     Xf=zeros(Nddl,1);
        %     Xf(1:Nddl)=Xg(1);
        %     Xf(Nddl+(1:Nddl))=Xg(2);
        
        %% for using RBT
        case 'rbt'
            Uo=1*(Xg(1));
            Vo=1*(Xg(2));
            Un=zeros(Nnodes)+Uo;
            Vn=zeros(Nnodes)+Vo;
            
            %% for using elementwise RBT
        case 'rbt_elem'
            
            Uo=1*round(Xg(1));
            Vo=1*round(Xg(2));
            if isfield(param0,'sampling_factor')
                psample=param0.sampling_factor;
            else
                psample=1;
            end
            pscale=2^(iscale-1);
            
            load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
            wo=256;
            wo=2*pscale*mean(param.mesh_size);
            load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamd-1),iscale-1)),'xo','yo','Nelems','conn');
            xo=(xo-0.5)*pscale*psample+0.5;
            yo=(yo-0.5)*pscale*psample+0.5;
            Nddl=prod(Nnodes);
            incn(1)=1;incn(2)=incn(1)*Nnodes(1);
            Un=zeros(Nnodes);
            Vn=zeros(Nnodes);
            Nn=zeros(Nnodes);
            %        load(fullfile('TMP',sprintf('sample%d_0',icamd-1)),'im0','sizeim','roi');
            load(fullfile('TMP',sprintf('sample%d',icamd-1)),'im0','sizeim','roi');
            xmin=1;ymin=1;
            xmax=sizeim(1);ymax=sizeim(2);
            ndim=length(sizeim);
            dec=roi(1:2:(2*ndim-1))-1;
            for i1=1:prod(Nelems)
                found=find(conn(i1,:)>0);
                inods=conn(i1,found);
                xn=xo(inods);yn=yo(inods);
                xg=round(mean(xn));yg=round(mean(yn));
                wmin0=min(min(min(abs(xmin-xg)),min(abs(xmax-xg))),min(min(abs(ymin-yg)),min(abs(ymax-yg))));
                wmin1=min(min(min(abs(xmin-xg-Uo)),min(abs(xmax-xg-Uo))),min(min(abs(ymin-yg-Vo)),min(abs(ymax-yg-Vo))));
                w=min(wmin1,min(wmin0,wo));
                hh=ceil(-w/2):floor(w/2);
                xp=min(max(dec(1)+xg+hh,1),sizeim(1));
                yp=min(max(dec(2)+yg+hh,1),sizeim(2));
                im00=im0(xp,yp);
                im11=im1(xp+Uo,yp+Vo);
                %              yp=ceil(min(yo(inods))):floor(max(yo(inods)));
                %              xp=ceil(min(xo(inods))):floor(max(xo(inods)));
                %             im00=im0(dec(1)+xp,dec(2)+yp);
                %             im11=im1(dec(1)+xp+Uo,dec(2)+yp+Vo);
                [U,V]=rbt(im00,im11);
                %            [U,V]=rbt2(im00,im11);
                Un(inods)=Un(inods)+U+Uo;
                Vn(inods)=Vn(inods)+V+Vo;
                Nn(inods)=Nn(inods)+1;
            end
            Un=Un./Nn;
            Vn=Vn./Nn;
            %               figure
            %               imagesc(Un)
            %               colorbar
            
            %% for using bilinear displacement / element
            % check=0;
            % Uo=1*round(Xg(1));
            % Vo=1*round(Xg(2));
            % if isfield(param0,'sampling_factor')
            %     psample=param0.sampling_factor;
            % else
            %     psample=1;
            % end
            %     wo=256;
            %     load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','Nelems','conn');
            %     pscale=2^(iscale-1);
            %     xo=(xo-0.5)*pscale*psample+0.5;
            %     yo=(yo-0.5)*pscale*psample+0.5;
            %     Nddl=prod(Nnodes);
            %     incn(1)=1;incn(2)=incn(1)*Nnodes(1);
            %     Un=zeros(Nnodes);
            %     Vn=zeros(Nnodes);
            %     Nn=zeros(Nnodes);
            % %        load(fullfile('TMP',sprintf('sample%d_0',icamd-1)),'im0','sizeim','roi');
            %         load(fullfile('TMP',sprintf('sample%d',icamd-1)),'im0','sizeim','roi');
            % dynamic=max(im0(:))-min(im0(:));
            %             sizeim1=size(im1);
            % xmin=1;ymin=1;
            %     xmax=sizeim(1);ymax=sizeim(2);
            %         ndim=length(sizeim);
            %         dec=roi(1:2:(2*ndim-1))-1;
            %         for i1=1:prod(Nelems)
            %             found=find(conn(i1,:)>0);
            %             inods=conn(i1,found);
            %             xn=xo(inods);yn=yo(inods);
            %             xg=round(mean(xn));yg=round(mean(yn));
            %         wmin0=min(min(min(abs(xmin-xg)),min(abs(xmax-xg))),min(min(abs(ymin-yg)),min(abs(ymax-yg))));
            %         wmin1=min(min(min(abs(xmin-xg-Uo)),min(abs(xmax-xg-Uo))),min(min(abs(ymin-yg-Vo)),min(abs(ymax-yg-Vo))));
            %         w=min(wmin1,min(wmin0,wo));
            %     hh=ceil(-w/2):floor(w/2);
            %              xp=min(max(dec(1)+xg+hh,1),sizeim(1));
            %              yp=min(max(dec(2)+yg+hh,1),sizeim(2));
            %             [Yp,Xp]=meshgrid(yp,xp);
            %             im00=im0(xp,yp);
            %             [gy,gx]=gradient(im00);
            %             Xpo=Xp(:)-mean(Xp(:));Ypo=Yp(:)-mean(Yp(:));
            %             Xpo=Xpo/max(abs(Xpo));Ypo=Ypo/max(abs(Ypo));
            %
            %             N=[0.25*(1-Xpo).*(1-Ypo),0.25*(1+Xpo).*(1-Ypo),0.25*(1+Xpo).*(1+Ypo),0.25*(1-Xpo).*(1+Ypo)];
            %             Nx=[N,0*N];
            %             Ny=[0*N,N];
            %             phidf=diag(gx(:))*Nx+diag(gy(:))*Ny;
            %             M=phidf'*phidf;
            %             X=[repmat(Uo,4,1);repmat(Vo,4,1)];
            %             res=Inf;iter=1;itermax=20;conv=1.e-3;
            %             while res>conv&&iter<itermax
            %                 Ux=Nx*X;Uy=Ny*X;
            %          disc=mexInterpLinear(min(sizeim1(1)-2,max(3,Xp(:)+Ux)),min(sizeim1(2)-2,max(3,Yp(:)+Uy)),im1);
            %                mean1=mean(disc(:));
            %         std1=std(disc(:));
            %         std0=std(im00(:));
            %             disc=disc-mean1;
            %         st=std0/std1;
            %
            %         disc=(im00(:)-st*disc(:));
            %                 F=phidf'*disc;
            % merror=mean(abs(disc(:)))/dynamic;
            %
            %             dX=M\F;
            %             res=norm(dX)/norm(X);
            %             if check
            %                 disp(sprintf('At iteration # %d',iter));
            %         disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
            %
            %             disp(sprintf('|dU|=%f',res));
            %             end
            %         X=X+dX;
            %         iter=iter+1;
            %             end
            %             Un(inods)=Un(inods)+mean(X(1:4));
            %             Vn(inods)=Vn(inods)+mean(X(5:8));
            %             Nn(inods)=Nn(inods)+1;
            % %         if check
            % %             pause
            % %         end
            %         end
            %          Un=Un./Nn;
            %         Vn=Vn./Nn;
            %   figure
            %   imagesc(Un)
            %   colorbar
            %% for using bilinear displacement
            % tic
            % check=0;
            % Uo=1*round(Xg(1));
            % Vo=1*round(Xg(2));
            % if isfield(param0,'sampling_factor')
            %     psample=param0.sampling_factor;
            % else
            %     psample=1;
            % end
            %     load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamd-1),iscale-1)),'xo','yo');
            %         load(fullfile('TMP',sprintf('sample%d_0',icamd-1)),'im0','sizeim','roi');
            % dynamic=max(im0(:))-min(im0(:));
            % im0=im0-mean(im0(:));
            %     pscale=2^(iscale-1);
            %     xo=max(1,min(sizeim(1),(xo-0.5)*pscale*psample+0.5));
            %     yo=max(1,min(sizeim(2),(yo-0.5)*pscale*psample+0.5));
            %
            %     p=1;
            %     [Nx]=NURBSBasisFunc(1+p,p,1:sizeim(1),sizeim(1)*(-p:p+1));
            %     [Ny]=NURBSBasisFunc(1+p,p,1:sizeim(2),sizeim(2)*(-p:p+1));
            % %     figure
            % %      plot(Nx)
            % %      hold on
            % %      plot(sum(Nx,2),'x')
            %
            %     N=Nx(:)*(Ny(:)');
            % [indx,indy,val]=find(N);
            % indxm=mod(indx,sizeim(1));
            % found=find(indxm==0);
            % indxm(found)=sizeim(1);
            % indym=mod(indy,sizeim(2));
            % found=find(indym==0);
            % indym(found)=sizeim(2);
            % indp=indxm+sizeim(1)*(indym-1);
            % indn=ceil(indx/(sizeim(1)))+(2+p-1)*(ceil(indy/(sizeim(2)))-1);
            % N=sparse(indp,indn,val,prod(sizeim),(2+p-1)^2);
            %             sizeim1=size(im1);
            %             gx=mexFDGradient(im0);
            %             im0=im0';
            %             gy=mexFDGradient(im0);
            %             im0=im0';
            %                 gy=gy';
            %                 [Yp,Xp]=meshgrid(1:sizeim(2),1:sizeim(1));
            % %             Xpo=Xp(:)-mean(Xp(:));Ypo=Yp(:)-mean(Yp(:));
            % %             Xpo=Xpo/max(abs(Xpo));Ypo=Ypo/max(abs(Ypo));
            % %
            % %             N=sparse([0.25*(1-Xpo).*(1-Ypo),0.25*(1+Xpo).*(1-Ypo),0.25*(1+Xpo).*(1+Ypo),0.25*(1-Xpo).*(1+Ypo)]);
            %             Nx=[N,0*N];
            %             Ny=[0*N,N];
            %             phidf=diag(sparse(gx(:)))*Nx+diag(sparse(gy(:)))*Ny;
            %             M=phidf'*phidf;
            %             X=[repmat(Uo,size(N,2),1);repmat(Vo,size(N,2),1)];
            %             res=Inf;iter=1;itermax=100;conv=1.e-3;
            %         std0=std(im0(:));
            %             while res>conv&&iter<itermax
            %                 Ux=Nx*X;Uy=Ny*X;
            %          disc=mexInterpLinear(min(sizeim1(1)-2,max(3,(Xp(:)-1)*psample+roi(1)+Ux)),min(sizeim1(2)-2,max(3,(Yp(:)-1)*psample+Uy+roi(3))),im1);
            %                mean1=mean(disc(:));
            %         std1=std(disc(:));
            %             disc=disc-mean1;
            %         st=std0/std1;
            %
            %         disc=(im0(:)-st*disc(:));
            %                 F=phidf'*disc;
            % merror=mean(abs(disc(:)))/dynamic;
            % %         figure
            % %         imagesc(reshape(abs(disc),sizeim))
            % %         colorbar
            % %         pause
            %             dX=M\F;
            %             res=norm(dX)/norm(X);
            %             if check
            %                 disp(sprintf('At iteration # %d',iter));
            %         disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
            %
            %             disp(sprintf('|dU|=%f',res));
            %             end
            %         X=X+dX;
            %         iter=iter+1;
            %             end
            %
            %                             Ux=Nx*X;Uy=Ny*X;
            %
            %          Un=reshape(interp2(reshape(Ux,sizeim),yo,xo,'*linear'),Nnodes);
            %          Vn=reshape(interp2(reshape(Uy,sizeim),yo,xo,'*linear'),Nnodes);
            % %   figure
            % %   imagesc(Un)
            % %   colorbar
            %    toc
            %% for using bilinear displacement on a coarsened image
        case {'bilin','bilin+'}
            check=0;
            Uo=1*round(Xg(1));
            Vo=1*round(Xg(2));
            load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
            nscale=param.nscale;
            load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamd-1),iscale-1)),'xo','yo');
            load(fullfile('TMP',sprintf('sample%d_%d',icamd-1,nscale-1)),'im0','sizeim','roi');
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod*10^(1-1),nscale-1)),'mask');
            dynamic=max(im0(:))-min(im0(:));
            %            im0=im0-mean(im0(:));
            pscalei=2^(iscale-1);
            pscalen=2^(nscale-1);
            roi=(roi-1)/pscalen+1;
            xo=((xo-0.5)*pscalei)/pscalen+0.5;
            yo=((yo-0.5)*pscalei)/pscalen+0.5;
            xo=max(1,min(sizeim(1),(xo-0.5)*1*psample+0.5));
            yo=max(1,min(sizeim(2),(yo-0.5)*1*psample+0.5));
            [Yo,Xo]=meshgrid([0,sizeim(2)],[0,sizeim(1)]);
            %             p=1;
            %             [Nx]=NURBSBasisFunc(1+p,p,1:sizeim(1),sizeim(1)*(-p:p+1));
            %             [Ny]=NURBSBasisFunc(1+p,p,1:sizeim(2),sizeim(2)*(-p:p+1));
            %     figure
            %      plot(Nx)
            %      hold on
            %      plot(sum(Nx,2),'x')
            %
            %             N=Nx(:)*(Ny(:)');
            %             [indx,indy,val]=find(N);
            %             indxm=mod(indx,sizeim(1));
            %             found=find(indxm==0);
            %             indxm(found)=sizeim(1);
            %             indym=mod(indy,sizeim(2));
            %             found=find(indym==0);
            %             indym(found)=sizeim(2);
            %             indp=indxm+sizeim(1)*(indym-1);
            %             indn=ceil(indx/(sizeim(1)))+(2+p-1)*(ceil(indy/(sizeim(2)))-1);
            %             N=sparse(indp,indn,val,prod(sizeim),(2+p-1)^2);
            
            im10=im1;
            NestedCoarseImage(nscale-1);
            sizeim1=size(im1);
            gx=mexFDGradient(im0);
            im0=im0';
            gy=mexFDGradient(im0);
            im0=im0';
            gy=gy';
            [Yp,Xp]=meshgrid(1:sizeim(2),1:sizeim(1));
            %             Xpo=Xp(:)-mean(Xp(:));Ypo=Yp(:)-mean(Yp(:));
            %             Xpo=Xpo/max(abs(Xpo));Ypo=Ypo/max(abs(Ypo));
            %
            %             N=sparse([0.25*(1-Xpo).*(1-Ypo),0.25*(1+Xpo).*(1-Ypo),0.25*(1+Xpo).*(1+Ypo),0.25*(1-Xpo).*(1+Ypo)]);
            %             Nx=[N,0*N];
            %             Ny=[0*N,N];
            Nx=[0*Xp(:)+1,0*Xp(:),Yp(:),Xp(:),0*Xp(:),0.5*Yp(:)];
            Ny=[0*Xp(:),0*Xp(:)+1,-Xp(:),0*Xp(:),Yp(:),0.5*Xp(:)];
            % Nx=[0*Xp(:)+1,0*Xp(:),Yp(:),Xp(:),0*Xp(:),0.5*Yp(:),...
            %     Xp(:).^2,Yp(:).^2,Xp(:).*Yp(:),...
            %     Xp(:).^3,Yp(:).^3,Xp(:).^2.*Yp(:),Xp(:).*Yp(:).^2,zeros(numel(Xp),7)];
            % Ny=[0*Xp(:),0*Xp(:)+1,-Xp(:),0*Xp(:),Yp(:),0.5*Xp(:),zeros(numel(Xp),7),...
            %         Xp(:).^2,Yp(:).^2,Xp(:).*Yp(:),...
            %     Xp(:).^3,Yp(:).^3,Xp(:).^2.*Yp(:),Xp(:).*Yp(:).^2];
            
            phidf=diag(sparse(gx(:)))*Nx+diag(sparse(gy(:)))*Ny;
            M=phidf'*mask*phidf;
            %             Uxo=N\(repmat(Uo/pscalen,size(N,1),1));
            %             Uyo=N\(repmat(Vo/pscalen,size(N,1),1));
            %             X=[Uxo;Uyo];
            X=[Uo/pscalen;Vo/pscalen;0;0;0;0];
            %            X=[repmat(Uo/pscalen,size(N,2),1);repmat(Vo/pscalen,size(N,2),1)];
            res=1*(norm(X)>eps);erprev=0;
            iter=1;itermax=500;conv=1.e-2;
            std0=std(im0(:));
            if check
                figure
                imagesc(im1)
                colormap('gray')
                hold on
                on=zeros(sizeim);
                nstep=5;
                xgrid=floor(1:(sizeim(1)-1)/nstep:sizeim(1));
                ygrid=floor(1:(sizeim(2)-1)/nstep:sizeim(2));
                
                %                on(1,:)=1;on(:,1)=1;on(sizeim(1),:)=1;on(:,sizeim(2))=1;
                on(xgrid,:)=1;on(:,ygrid)=1;
                on=on==1;
            end
            if res
                while (res>conv&&iter<itermax)||(iter<10)
                    Ux=Nx*X;Uy=Ny*X;
                    Xi=(Xp(:)-1)*psample+roi(1)+Ux;
                    Yi=(Yp(:)-1)*psample+Uy+roi(3);
                    if check
                        %                     xc=Xo(:)+X((1:length(Xo(:))))+roi(1)-1;
                        %                     yc=Yo(:)+X(length(Xo(:))+(1:length(Xo(:))))+roi(3)-1;
                        %                     line(yc([1,2,4,3,1]),xc([1,2,4,3,1]),'LineWidth',2,'Color',[0,0,1]);
                        %                     axis image
                        %                     axis ij
                        xc=Xi(on);
                        yc=Yi(on);
                        plot(yc,xc,'.','LineWidth',1,'Color',[0,0,1]);
                        axis image
                        axis ij
                        pause(0.1)
                    end
                    disc=mexInterpLinear(Xi,Yi,im1);
                    maski=disc(:)<0;%double((Xi>1)&(Yi>1)&(Xi<sizeim1(1))&(Yi<sizeim1(2)));
                    %                maski=diag(sparse(maski(:)));
                    %                maski=mask*maski;
                    mean1=mean(disc(:));
                    std1=std(disc(:));
                    %disc=disc-mean1;
                    st=1;%std0/std1;
                    disc=(im0(:)-st*disc(:));
                    disc(maski)=0;
                    F=phidf'*(mask*disc);
                    merror=mean(abs(disc(:)))/dynamic;
                    %                 if check
                    %                     figure
                    %                     imagesc(reshape(abs(maski*disc),sizeim))
                    %                     colorbar
                    %                     pause
                    %                 end
                    dX=M\F;
                    resu=norm(dX)/norm(X);
                    if iter==1;
                        ero=merror;
                    end
                    resi=abs((merror-erprev)/ero);
                    res=resu;
                    if check
                        disp(sprintf('At iteration # %d',iter));
                        disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
                        
                        disp(sprintf('|dU/U|=%f   |phi/phio|=%f',resu,resi));
                    end
                    X=X+dX;
                    iter=iter+1;
                    erprev=merror;
                end
                if check
                    disp(sprintf('End of intialization'))
                    pause
                end
                if res>conv
                    disp(sprintf('CONVERGENCE FAILUIRE IN THE INITIALIZATION FOR CAMERA %d',icamd))
                    disp(sprintf('At iteration # %d',iter));
                    disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
                    disp(sprintf('|dU/U|=%f   |phi/phio|=%f',resu,resi));
                    disp(sprintf('THE RIGID BODY TRANSLATION IS USED INSTEAD'))
                    
                    figure
                    imagesc(im1)
                    colormap('gray')
                    hold on
                    %                 xc=Xo(:)+X((1:length(Xo(:))))+roi(1)-1;
                    %                 yc=Yo(:)+X(length(Xo(:))+(1:length(Xo(:))))+roi(3)-1;
                    %
                    %                 line(yc([1,2,4,3,1]),xc([1,2,4,3,1]),'LineWidth',2,'Color',[0,0,1]);
                    on=zeros(sizeim);
                    nstep=5;
                    xgrid=floor(1:(sizeim(1)-1)/nstep:sizeim(1));
                    ygrid=floor(1:(sizeim(2)-1)/nstep:sizeim(2));
                    
                    %                on(1,:)=1;on(:,1)=1;on(sizeim(1),:)=1;on(:,sizeim(2))=1;
                    on(xgrid,:)=1;on(:,ygrid)=1;
                    on=on==1;
                    xc=Xi(on);
                    yc=Yi(on);
                    plot(yc,xc,'.','LineWidth',1,'Color',[0,0,1]);
                    
                    X=[Uo/pscalen;Vo/pscalen;0;0;0;0]*0;
                    Ux=Nx*X;Uy=Ny*X;
                    Xi=(Xp(:)-1)*psample+roi(1)+Ux;
                    Yi=(Yp(:)-1)*psample+Uy+roi(3);
                    
                    on=zeros(sizeim);
                    nstep=5;
                    xgrid=floor(1:(sizeim(1)-1)/nstep:sizeim(1));
                    ygrid=floor(1:(sizeim(2)-1)/nstep:sizeim(2));
                    
                    %                on(1,:)=1;on(:,1)=1;on(sizeim(1),:)=1;on(:,sizeim(2))=1;
                    on(xgrid,:)=1;on(:,ygrid)=1;
                    on=on==1;
                    xc=Xi(on);
                    yc=Yi(on);
                    plot(yc,xc,'.','LineWidth',1,'Color',[1,0,0]);
                    
                    
                    
                    axis image
                    axis ij
                end
            end
            %                     figure
            %                     imagesc(im1)
            %              colormap('gray')
            %             hold on
            %                     xc=Xo(:)+X((1:length(Xo(:))))+roi(1)-1;
            %                     yc=Yo(:)+X(length(Xo(:))+(1:length(Xo(:))))+roi(3)-1;
            %
            %                 line(yc([1,2,4,3,1]),xc([1,2,4,3,1]),'LineWidth',2,'Color',[0,0,1]);
            % axis image
            % axis ij
            
            Ux=Nx*X*pscalen;Uy=Ny*X*pscalen;
            %             if check
            %                 figure
            %                 imagesc(reshape(Ux,sizeim))
            %                 colorbar
            %             end
            Un=reshape(interp2(reshape(Ux,sizeim),yo,xo,'*linear'),Nnodes);
            Vn=reshape(interp2(reshape(Uy,sizeim),yo,xo,'*linear'),Nnodes);
            %              toc
            %
            im1=im10;
            clear im10
            if strcmp(init,'bilin+')
                
                
                Uo=(Un);
                Vo=(Vn);
                if isfield(param0,'sampling_factor')
                    psample=param0.sampling_factor;
                else
                    psample=1;
                end
                pscale=2^(iscale-1);
                
                load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
                wo=pscale*mean(param.mesh_size);
                load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamd-1),iscale-1)),'xo','yo','Nelems','conn');
                xo=(xo-0.5)*pscale*psample+0.5;
                yo=(yo-0.5)*pscale*psample+0.5;
                Nddl=prod(Nnodes);
                incn(1)=1;incn(2)=incn(1)*Nnodes(1);
                Un=zeros(Nnodes);
                Vn=zeros(Nnodes);
                Nn=zeros(Nnodes);
                %        load(fullfile('TMP',sprintf('sample%d_0',icamd-1)),'im0','sizeim','roi');
                load(fullfile('TMP',sprintf('sample%d',icamd-1)),'im0','sizeim','roi');
                xmin=1;ymin=1;
                xmax=sizeim(1);ymax=sizeim(2);
                ndim=length(sizeim);
                dec=roi(1:2:(2*ndim-1))-1;
                for i1=1:prod(Nelems)
                    found=find(conn(i1,:)>0);
                    inods=conn(i1,found);
                    xn=xo(inods);yn=yo(inods);
                    xg=round(mean(xn));yg=round(mean(yn));
                    Ue=Uo(inods);Ve=Vo(inods);
                    Ue=round(mean(Ue));Ve=round(mean(Ve));
                    wmin0=min(min(min(abs(xmin-xg)),min(abs(xmax-xg))),min(min(abs(ymin-yg)),min(abs(ymax-yg))));
                    wmin1=min(min(min(abs(xmin-xg-Uo)),min(abs(xmax-xg-Uo))),min(min(abs(ymin-yg-Vo)),min(abs(ymax-yg-Vo))));
                    w=min(wmin1,min(wmin0,wo));
                    hh=ceil(-w/2):floor(w/2);
                    xp=min(max(dec(1)+xg+hh,1),sizeim(1));
                    yp=min(max(dec(2)+yg+hh,1),sizeim(2));
                    im00=im0(xp,yp);
                    im11=im1(xp+Ue,yp+Ve);
                    %              yp=ceil(min(yo(inods))):floor(max(yo(inods)));
                    %              xp=ceil(min(xo(inods))):floor(max(xo(inods)));
                    %             im00=im0(dec(1)+xp,dec(2)+yp);
                    %             im11=im1(dec(1)+xp+Uo,dec(2)+yp+Vo);
                    [U,V]=rbt(im00,im11);
                    %            [U,V]=rbt2(im00,im11);
                    Un(inods)=Un(inods)+U+Ue;
                    Vn(inods)=Vn(inods)+V+Ve;
                    Nn(inods)=Nn(inods)+1;
                end
                Un=Un./Nn;
                Vn=Vn./Nn;
                
                
            end
            
            
            
            
            
            
            
            
            
            
            
            %     load(fullfile('TMP',sprintf('sample%d,icamd-1)),'im0','sizeim');
            %         ndim=length(sizeim);
            %         dec=roi(1:2:(2*ndim-1))-1;
            %     wsize=50;
            %     for i1=1:prod(Nelems)
            %         found=find(conn(i1,:)>0);
            %         inods=conn(i1,found);
            %         xn=round(mean(xo(inods)));yn=round(mean(yo(inods)));
            %         yp=(yn-wsize/2):(yn+wsize/2);
            %         xp=(xn-wsize/2):(xn+wsize/2);
            %         im00=im0(dec(1)+xp,dec(2)+yp);
            %         im11=im1(dec(1)+xp,dec(2)+yp);
            %         [U,V]=rbt2(im00,im11);
            %         Un(inods)=Un(inods)+U;
            %         Vn(inods)=Vn(inods)+V;
            %         Nn(inods)=Nn(inods)+1;
            %     end
            %      Un=Un./Nn;
            %      Vn=Vn./Nn;
            %      figure
            %      scatter(yo(:),xo(:),20+0*yo(:),Un(:))
            
            
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
            
            if ~isempty(unmasked_nodes)
                Un=Un(unmasked_nodes);
                Vn=Vn(unmasked_nodes);
            end
    end
    Xf=[Un(:);Vn(:)];
    if thermo
        Xf=[Xf;0*Vn(:)];
    end
    %            figure
    %            imagesc(Un)
    %            colorbar;
    
else
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale)),'unmasked_nodes');
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale)),'xo','yo','Nnodes','conn','elt');
    xg=2*(xo-0.5)+0.5;
    yg=2*(yo-0.5)+0.5;
    
    meshg.conn=conn;
    meshg.elt=elt;
    meshg.xo=xg;
    meshg.yo=yg;
    Tg=[];
    xg=reshape(xg,Nnodes(1:2));
    yg=reshape(yg,Nnodes(1:2));
    if isempty(unmasked_nodes)
        Nddl=prod(Nnodes);
        Ug=reshape(Xg((1:Nddl)),Nnodes(1:2));
        Vg=reshape(Xg(Nddl+(1:Nddl)),Nnodes(1:2));
        if thermo
            Tg=reshape(Xg(2*Nddl+(1:Nddl)),Nnodes(1:2));
        end
    else
        Nddl=length(unmasked_nodes);
        Ug=repmat(0,Nnodes(1:2));
        Vg=repmat(0,Nnodes(1:2));
        Ug(unmasked_nodes)= Xg((1:Nddl));
        Vg(unmasked_nodes)= Xg(Nddl+(1:Nddl));
        if thermo
            Tg=repmat(0,Nnodes(1:2));
            Tg(unmasked_nodes)= Tg((1:Nddl));
        end
    end
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','Nnodes');
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
    
    xo=min(xo,max(xg(:)));xo=max(xo,min(xg(:)));
    yo=min(yo,max(yg(:)));yo=max(yo,min(yg(:)));
    Uf= griddata(xg,yg,Ug,xo(:),yo(:));
    Vf= griddata(xg,yg,Vg,xo(:),yo(:));
    %     Uf=interp2(yg,xg,Ug,yo(:),xo(:),'linear');
    %     Vf=interp2(yg,xg,Vg,yo(:),xo(:),'linear');
    
    %     Ug=Ug(:);
    %     Vg=Vg(:);
    %     Tg=Tg(:);
    %     coords.xi=xo;
    %     coords.yi=yo;
    %     [UVf]=interpMesh(meshg,[Ug,Vg,Tg],coords);
    %     Uf=UVf(:,1);
    %     Vf=UVf(:,2);
    
    if isempty(unmasked_nodes)
        Xf=[Uf(:);Vf(:)];
    else
        Xf=[Uf(unmasked_nodes);Vf(unmasked_nodes)];
    end
    if thermo
%        Tf=UVf(:,3);
    Tf= griddata(xg,yg,Tg,xo(:),yo(:));
        if isempty(unmasked_nodes)
            Xf=[Xf;Tf(:)];
        else
            Xf=[Xf;Tf(unmasked_nodes)];
        end
    end
end
    function NestedCoarseImage(scales)
        
        scale=2;
        for ip=1:scales
            imsiz0=size(im1);
            imsiz1=floor(imsiz0/2);
            nn=2*imsiz1;
            im1=im1(1:nn(1),1:nn(2));
            
            im1=reshape(im1,scale,prod(nn)/scale);
            im1=mean(im1,1);
            nn(1)=nn(1)/scale;
            im1=reshape(im1,nn);
            
            im1=im1';
            im1=reshape(im1,scale,prod(nn)/scale);
            im1=mean(im1,1);
            nn(2)=nn(2)/scale;
            im1=reshape(im1,nn([2,1]));
            im1=im1';
        end
        %toc();
        
    end


end