function [Uini]=InitializeSolution3D(v1,nmod,v3)
check=1;
if nargin<3
    h=v1;
    iscale=1;
else
    U=v1;
    iscale=v3;
end
load(fullfile('TMP','params'));
param0=param;
if iscell(param0.deformed_image)
    nim=length(param0.deformed_image);
else
    nim=1;
end
inorm=false;
if isfield(param0,'normalize_grey_level')
    inorm=param0.normalize_grey_level;
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
roi=param0.roi;
restart=1;
if isfield(param0,'restart')
    restart=param0.restart;
end
if restart
    nmax=nim;
else
    if nim<2
        nmax=nim;
    else
        nmax=2;
    end
end

if nargin==2
    load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','Nelems','Nnodes','conn');
    incn(1)=1;incn(2)=incn(1)*Nnodes(1);incn(3)=incn(2)*Nnodes(2);
    load(fullfile('TMP','sample0'),'sizeim');
    switch param0.stack_format
        case 'bin'
            fid=fopen(fullfile('TMP',sprintf('dsample0_%d',iscale-1)));
            im0=fread(fid,prod(sizeim));
            fclose(fid);
            im0=reshape(im0,sizeim);
        case 'mat'
            load(fullfile('TMP','sample0'),'im0');
    end
    dec=roi(1:2:(2*3-1))-1;
    hh=ceil(-h/2):floor(h/2);
    sizeim=size(im0);
    
    Uini=zeros(3*prod(Nnodes),nim);
    for iim=1:nim
        if nim==1
            fildef=param0.deformed_image;
        else
            fildef=param0.deformed_image{iim};
        end
        
        if isfield(param0,'stack_size')
            imsiz0=param0.stack_size;
            fid=fopen(fildef,'r');
            jm3=fread(fid,prod(imsiz0));
            jm3=reshape(jm3,imsiz0);
            fclose(fid);
        else
            [~, ~, ext] = fileparts(fildef);
            switch ext
                case '.mat'
                    load(fildef,'jm3');
                case {'.tif','.tiff'}
                    jm3=readTIFFasRAW(fildef);
            end
        end
        
        
        [Urbto]=rbt3(im0(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)),jm3(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)),inorm);
        Uel=round(Urbto(1));Vel=round(Urbto(2));Wel=round(Urbto(3));
        Un=zeros(Nnodes);
        Vn=zeros(Nnodes);
        Wn=zeros(Nnodes);
        Nn=zeros(Nnodes);
        for i1=1:prod(Nelems)
            found=find(conn(i1,:)>0);
            inods=conn(i1,found);
            xn=xo(inods);yn=yo(inods);zn=zo(inods);
            xg=round(mean(xn));yg=round(mean(yn));zg=round(mean(zn));
            xp=dec(1)+xg+hh;
            found=find((xp>=1)&(xp<=sizeim(1))&(xp+Uel>=1)&(xp+Uel<=sizeim(1)));
            xp=xp(found);
            yp=dec(2)+yg+hh;
            found=find((yp>=1)&(yp<=sizeim(2))&(yp+Vel>=1)&(yp+Vel<=sizeim(2)));
            yp=yp(found);
            zp=dec(3)+zg+hh;
            found=find((zp>=1)&(zp<=sizeim(3))&(zp+Wel>=1)&(zp+Wel<=sizeim(3)));
            zp=zp(found);
            im00=im0(xp,yp,zp);
            im11=jm3(xp+Uel,yp+Vel,zp+Wel);
            [Urbt]=rbt3(im00,im11);
            Un(inods)=Un(inods)+Uel+Urbt(1);
            Vn(inods)=Vn(inods)+Vel+Urbt(2);
            Wn(inods)=Wn(inods)+Wel+Urbt(3);
            Nn(inods)=Nn(inods)+1;
        end
        
        Ut=Un./Nn;
        Vt=Vn./Nn;
        Wt=Wn./Nn;
        
        
        Uini(:,iim)=[Ut(:);Vt(:);Wt(:)];
    end
else
    
    if isempty(U)||size(U,1)==3
        dorbm=0;
        if isfield(param0,'dorbm')
            dorbm=param0.dorbm;
        end
        pscale=2^(iscale-1);
        if dorbm
            if isempty(U)
                dorbt=1;
            else
                dorbt=0;
                Urbto=U;
            end
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            switch param0.stack_format
                case 'bin'
                    fid=fopen(fullfile('TMP',sprintf('dsample0_%d',iscale-1)));
                    im0=fread(fid,prod(sizeim));
                    fclose(fid);
                    im0=reshape(im0,sizeim);
                case 'mat'
                    load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
            end
            
            mesh_file=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1));
            load(mesh_file,'xo','yo','zo','Nnodes');
            mim0=mean(im0(:));
            sim0=max(1,std(im0(:)));
            Uini=zeros(3*prod(Nnodes),nim);
            [Yi,Xi,Zi]=meshgrid(1:sizeim(2),1:sizeim(1),1:sizeim(3));
            xc=0.5*(1+sizeim(1));yc=0.5*(1+sizeim(2));zc=0.5*(1+sizeim(3));
            phix=[1+0*Xi(:),0*Xi(:),0*Xi(:),0*Xi(:),-Zi(:)+zc,Yi(:)-yc];
            phiy=[0*Xi(:),1+0*Xi(:),0*Xi(:),Zi(:)-zc,0*Zi(:),-Xi(:)+xc];
            phiz=[0*Xi(:),0*Xi(:),1+0*Xi(:),-Yi(:)+yc,Xi(:)-xc,0*Xi(:)];
            phixn=[1+0*xo(:),0*xo(:),0*xo(:),0*xo(:),-zo(:)+zc,yo(:)-yc];
            phiyn=[0*xo(:),1+0*xo(:),0*xo(:),zo(:)-zc,0*zo(:),-xo(:)+xc];
            phizn=[0*xo(:),0*xo(:),1+0*xo(:),-yo(:)+yc,xo(:)-xc,0*xo(:)];
            [gy,gx,gz]=gradient(im0);
            phidf=diag(sparse(gx(:)))*phix+diag(sparse(gy(:)))*phiy+diag(sparse(gz(:)))*phiz;
            clear gx gy gz
            roip=(roi-1)/pscale+1;
            dynamic=max(im0(:))-min(im0(:));
            scale=2;
            for iim=1:nmax
                if nim==1
                    fildef=param0.deformed_image;
                else
                    fildef=param0.deformed_image{iim};
                end
                if isfield(param0,'stack_size')
                    imsiz0=param0.stack_size;
                    fid=fopen(fildef,'r');
                    jm3=fread(fid,prod(imsiz0));
                    jm3=reshape(jm3,imsiz0);
                    fclose(fid);
                else
                    [~, ~, ext] = fileparts(fildef);
                    switch ext
                        case '.mat'
                            load(fildef,'jm3');
                        case {'.tif','.tiff'}
                            jm3=readTIFFasRAW(fildef);
                    end
                end
                if dorbt
                    jm31=jm3((roi(1):roi(2)),(roi(3):roi(4)),(roi(5):roi(6)));
                    
                    for ip=1:(iscale-1)
                        imsiz0=size(jm31);
                        imsiz1=floor(imsiz0/2);
                        nn=2*imsiz1;
                        jm31=jm31(1:nn(1),1:nn(2),1:nn(3));
                        
                        jm31=reshape(jm31,scale,prod(nn)/scale);
                        jm31=mean(jm31,1);
                        nn(1)=nn(1)/scale;
                        jm31=reshape(jm31,nn);
                        
                        jm31=permute(jm31,[2,3,1]);
                        jm31=reshape(jm31,scale,prod(nn)/scale);
                        jm31=mean(jm31,1);
                        nn(2)=nn(2)/scale;
                        jm31=reshape(jm31,nn([2,3,1]));
                        jm31=permute(jm31,[3,1,2]);
                        
                        jm31=permute(jm31,[3,1,2]);
                        jm31=reshape(jm31,scale,prod(nn)/scale);
                        jm31=mean(jm31,1);
                        nn(3)=nn(3)/scale;
                        jm31=reshape(jm31,nn([3,1,2]));
                        jm31=permute(jm31,[2,3,1]);
                    end
                    [Urbt]=rbt3(im0,jm31,inorm);
                    clear jm31
                else
                    Urbt=Urbto(:,iim)/pscale;
                end
                
                for ip=1:(iscale-1)
                    imsiz0=size(jm3);
                    imsiz1=floor(imsiz0/2);
                    nn=2*imsiz1;
                    jm3=jm3(1:nn(1),1:nn(2),1:nn(3));
                    
                    jm3=reshape(jm3,scale,prod(nn)/scale);
                    jm3=mean(jm3,1);
                    nn(1)=nn(1)/scale;
                    jm3=reshape(jm3,nn);
                    
                    jm3=permute(jm3,[2,3,1]);
                    jm3=reshape(jm3,scale,prod(nn)/scale);
                    jm3=mean(jm3,1);
                    nn(2)=nn(2)/scale;
                    jm3=reshape(jm3,nn([2,3,1]));
                    jm3=permute(jm3,[3,1,2]);
                    
                    jm3=permute(jm3,[3,1,2]);
                    jm3=reshape(jm3,scale,prod(nn)/scale);
                    jm3=mean(jm3,1);
                    nn(3)=nn(3)/scale;
                    jm3=reshape(jm3,nn([3,1,2]));
                    jm3=permute(jm3,[2,3,1]);
                end
                sizeim1=size(jm3);
                iter=1;itermax=100;conv=1.e-3;
                res=Inf;
                X=[Urbt;zeros(3,1)];
                while (res>conv&&iter<itermax)
                    Ux=phix*X;Uy=phiy*X;Uz=phiz*X;
                    Xp=(Xi(:)-1)+roip(1)+Ux;
                    Yp=(Yi(:)-1)+Uy+roip(3);
                    Zp=(Zi(:)-1)+Uz+roip(5);
                    maski=double((Xp>1)&(Yp>1)&(Zp>1)&(Xp<sizeim1(1))&(Yp<sizeim1(2))&(Zp<sizeim1(3)));
                    maski=diag(sparse(maski(:)));
                    disc=mexInterpLinear3D(Xp,Yp,Zp,jm3);
                    if inorm
                        disc=(im0(:)-mim0)-sim0/max(1,std(disc(:)))*(disc-mean(disc(:)));
                    else
                        disc=(im0(:)-disc(:));
                    end
                    F=phidf'*(maski*disc);
                    M=phidf'*(maski*phidf);
                    merror=mean(abs(disc(:)))/dynamic;
                    dX=M\F;
                    res=norm(dX)/norm(X);
                    if check
                        disp(sprintf('At iteration # %d',iter));
                        disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
                        disp(sprintf('|dU/U|=%f',res));
                    end
                    X=X+dX;
                    iter=iter+1;
                end
                if res>conv
                    warning('Convergence not reached after %d iterations',iter);
                    disp(sprintf('Discrepancy wrt dyn. =%6.2f %%',merror*100));
                    disp(sprintf('|dU/U|=%f',res));
                end
                Uxn=phixn*X;
                Uyn=phiyn*X;
                Uzn=phizn*X;
                Uini(:,iim)=[Uxn(:);Uyn(:);Uzn(:)]*pscale;
                
            end
        else
            switch param.basis
                case 'fem'
                    load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'Nnodes');
                case 'nurbs'
                    load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'Nbsnodes');
                    Nnodes=Nbsnodes;
            end
            if isempty(U)
                Urbto=zeros(3,nim);
            else
                Urbto=round(U);
            end
            dorbt=1;
            if isfield(param0,'dorbt')
                dorbt=param0.dorbt;
            end
            if dorbt
                Urbt=zeros(3,nim);
                load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
                switch param0.stack_format
                    case 'bin'
                        fid=fopen(fullfile('TMP',sprintf('dsample0_%d',iscale-1)));
                        im0=fread(fid,prod(sizeim));
                        fclose(fid);
                        im0=reshape(im0,sizeim);
                    case 'mat'
                        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
                end
                scale=2;
                
                for iim=1:nmax
                    if nim==1
                        fildef=param0.deformed_image;
                    else
                        fildef=param0.deformed_image{iim};
                    end
                    if isfield(param0,'stack_size')
                        imsiz0=param0.stack_size;
                        fid=fopen(fildef,'r');
                        jm3=fread(fid,prod(imsiz0));
                        jm3=reshape(jm3,imsiz0);
                        fclose(fid);
                    else
                        [~, ~, ext] = fileparts(fildef);
                        switch ext
                            case '.mat'
                                load(fildef,'jm3');
                            case {'.tif','.tiff'}
                                jm3=readTIFFasRAW(fildef);
                        end
                    end
                    jm3=jm3(Urbto(1,iim)+(roi(1):roi(2)),Urbto(2,iim)+(roi(3):roi(4)),Urbto(3,iim)+(roi(5):roi(6)));
                    for ip=1:(iscale-1)
                        imsiz0=size(jm3);
                        imsiz1=floor(imsiz0/2);
                        nn=2*imsiz1;
                        jm3=jm3(1:nn(1),1:nn(2),1:nn(3));
                        
                        jm3=reshape(jm3,scale,prod(nn)/scale);
                        jm3=mean(jm3,1);
                        nn(1)=nn(1)/scale;
                        jm3=reshape(jm3,nn);
                        
                        jm3=permute(jm3,[2,3,1]);
                        jm3=reshape(jm3,scale,prod(nn)/scale);
                        jm3=mean(jm3,1);
                        nn(2)=nn(2)/scale;
                        jm3=reshape(jm3,nn([2,3,1]));
                        jm3=permute(jm3,[3,1,2]);
                        
                        jm3=permute(jm3,[3,1,2]);
                        jm3=reshape(jm3,scale,prod(nn)/scale);
                        jm3=mean(jm3,1);
                        nn(3)=nn(3)/scale;
                        jm3=reshape(jm3,nn([3,1,2]));
                        jm3=permute(jm3,[2,3,1]);
                    end
                    
                    [Urbti]=rbt3(im0,jm3,inorm);
                    Urbt(:,iim)=Urbti*2^(iscale-1)+Urbto(:,iim);
                end
            else
                Urbt=Urbto;
            end
            Ut=repmat(Urbt(1,:),prod(Nnodes),1);
            Vt=repmat(Urbt(2,:),prod(Nnodes),1);
            Wt=repmat(Urbt(3,:),prod(Nnodes),1);
            Uini=[Ut;Vt;Wt];
            
        end
        
    else
        
        
        load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale)),'xo','yo','zo','Nnodes','conn','elt');
        xg=2*(xo-0.5)+0.5;
        yg=2*(yo-0.5)+0.5;
        zg=2*(zo-0.5)+0.5;
        
        meshg.conn=conn;
        meshg.elt=elt;
        meshg.xo=xg;
        meshg.yo=yg;
        meshg.zo=zg;
        
        xg=reshape(xg,Nnodes);
        yg=reshape(yg,Nnodes);
        zg=reshape(zg,Nnodes);
        Nddl=prod(Nnodes);
        Ng=Nnodes;
        load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','Nnodes');
        Uini=zeros(3*prod(Nnodes),nim);
        xo=min(xo,max(xg(:)));xo=max(xo,min(xg(:)));
        yo=min(yo,max(yg(:)));yo=max(yo,min(yg(:)));
        zo=min(zo,max(zg(:)));zo=max(zo,min(zg(:)));
        
        coords.xi=xo;
        coords.yi=yo;
        coords.zi=zo;
        
        UVWg=[U((1:Nddl),:),U(Nddl+(1:Nddl),:),U(2*Nddl+(1:Nddl),:)];
        [UVWf]=interpMesh(meshg,UVWg,coords);
        Uini=[UVWf(:,(1:nmax));UVWf(:,nim+(1:nmax));UVWf(:,2*nim+(1:nmax))];
        
        %         for iim=1:nim
        %             Ug=reshape(U((1:Nddl),iim),Ng);
        %             Vg=reshape(U(Nddl+(1:Nddl),iim),Ng);
        %             Wg=reshape(U(2*Nddl+(1:Nddl),iim),Ng);
        %
        %             Uf=interp3(yg,xg,zg,Ug,yo(:),xo(:),zo(:),'linear');
        %             Vf=interp3(yg,xg,zg,Vg,yo(:),xo(:),zo(:),'linear');
        %             Wf=interp3(yg,xg,zg,Wg,yo(:),xo(:),zo(:),'linear');
        %
        %             Uini(:,iim)=[Uf(:);Vf(:);Wf(:)];
        %
        %         end
        
        
        
        
        if iscale==1&&~strcmp(param.basis,'fem')
            load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'Nbsnodes');
            load(fullfile('TMP',sprintf('%d_3d_phio_%d',nmod,iscale-1)),'phio');
            Uinio=Uini;
            Uini=zeros(3*prod(Nbsnodes),nim);
            for iim=1:nmax
                Xg=Uinio(:,iim);
                Ng=size(Xg,1)/3;
                Ug=Xg((1:Ng),:);
                Vg=Xg(Ng+(1:Ng),:);
                Wg=Xg(2*Ng+(1:Ng),:);
                M=phio'*phio;
                F=phio'*Ug;
                Uf=M\F;
                F=phio'*Vg;
                Vf=M\F;
                F=phio'*Wg;
                Wf=M\F;
                Uini(:,iim)=[Uf;Vf;Wf];
            end
            
        end
        if nim>nmax
            Uini=[Uini,zeros(size(Uini,1),nim-nmax)];
        end
    end
    
    
    
end

