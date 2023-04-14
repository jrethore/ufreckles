function RegularMeshes(nmod,check)

if nargin<2,check=0;end
rfac=1.e-3;
rflag=true;
load(fullfile('TMP','params'),'param');
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
rint=false;
if isfield(param,'reduced_integration')
    rint=param.reduced_integration;
end
if isfield(param,'element_type')
    eelt=param.element_type;
else
    eelt=4;
end
if strcmp(param.basis,'btri')
    eelt=33;
end
nscale=param.nscale;
if isfield(param,'nb_gauss_points')
    ngi=param.nb_gauss_points;
else
    ngi=0;
end
nsi=ones(1,3);
if isfield(param,'nb_sub_cells')
    nsi=param.nb_sub_cells;
end
for iscale=1:nscale
    mesh_size=(param.mesh_size)/psample;
    load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
    
    if iscale>1
        ng=0;ns=ones(1,3);
    else
        ng=ngi;ns=nsi;
    end
    xo=(1:mesh_size(1):sizeim(1)-1)+1;
    yo=(1:mesh_size(2):sizeim(2)-1)+1;
    
    if length(sizeim)==3
        zo=(1:mesh_size(3):sizeim(3)-1)+1;
        if (length(xo)==1)
            mesh_size(1)=sizeim(1)-1;
            xo=round(0:mesh_size(1):sizeim(1)-1)+1;
        end
        if (length(yo)==1)
            mesh_size(2)=sizeim(2)-1;
            yo=round(0:mesh_size(2):sizeim(2)-1)+1;
        end
        if (length(zo)==1)
            mesh_size(3)=sizeim(3)-1;
            zo=round(0:mesh_size(3):sizeim(3)-1)+1;
        end
        zo=zo-floor(mean(zo)-sizeim(3)/2);
        
    else
        zo=1;
        if (length(xo)==1)||(length(yo)==1)
            if (length(xo)==1)&&(length(yo)==1)
                mesh_size(1)=sizeim(1)-1;
                mesh_size(2)=sizeim(2)-1;
            elseif length(xo)==1
                mesh_size(1)=sizeim(1)-1;
            elseif length(yo)==1
                mesh_size(2)=sizeim(2)-1;
            end
            xo=round(1:mesh_size(1):sizeim(1)-1)+1;
            yo=round(1:mesh_size(2):sizeim(2)-1)+1;
            
        end
    end
    if min(mesh_size)>2
    xo=xo-floor(mean(xo)-sizeim(1)/2);
    yo=yo-floor(mean(yo)-sizeim(2)/2);
    end
    
    
    Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
    Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
    Nnodes=[length(xo),length(yo),length(zo)];
    Nelems=max(Nnodes-1,1);
    if length(sizeim)==3
        [yo,xo,zo]=meshgrid(yo,xo,zo);
        yo=yo(:)-0.5;
        xo=xo(:)-0.5;
        zo=zo(:)-0.5;
        
    if eelt==4||iscale>1
            incn=[0,1,Nnodes(1)+1,Nnodes(1)];
            incn=[incn,incn+prod(Nnodes(1:2))];
            nroot=reshape(1:prod(Nnodes),Nnodes);
            nroot=nroot(1:max(1,Nnodes(1)-1),:,:);
            nroot=nroot(:,1:max(1,Nnodes(2)-1),:);
            nroot=nroot(:,:,1:max(1,Nnodes(3)-1));
            conn=repmat(nroot(:),1,8)+repmat(incn,prod(Nelems),1);
            elt=repmat(8,prod(Nelems),1);
        else
            conn=delaunay(xo,yo,zo);
            conn=[conn,zeros(size(conn,1),4)];
            Nelems=[size(conn,1),1,1];
            elt=repmat(4,prod(Nelems),1);
            
        end
        selected=ones(Nnodes);
        selected(1,:,:)=0;
        selected(:,1,:)=0;
        selected(:,:,1)=0;
        selected(Nnodes(1),:,:)=0;
        selected(:,Nnodes(2),:)=0;
        selected(:,:,Nnodes(3))=0;
        save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'rflag','rint','selected','ns','ng','xo','yo','zo','Nnodes','Nelems','Smesh','Vmesh','mesh_size','conn','elt');
        
        
    else
        [yo,xo]=meshgrid(yo,xo);
        yo=yo(:)-0.5;
        xo=xo(:)-0.5;
        incn=[0,1,Nnodes(1)+1,Nnodes(1)];
        nroot=reshape(1:prod(Nnodes),Nnodes);
        nroot=nroot(1:max(1,Nnodes(1)-1),:);
        nroot=nroot(:,1:max(1,Nnodes(2)-1));
        conn=repmat(nroot(:),1,4)+repmat(incn,prod(Nelems),1);
        elt=repmat(4,prod(Nelems),1);
        selected=ones(Nnodes);
        selected(1,:)=0;
        selected(:,1)=0;
        selected(Nnodes(1),:)=0;
        selected(:,Nnodes(2))=0;
        
        if (eelt==3||eelt==33)&&iscale==1
            %figure
            %plot(xo,yo,'bo')
            %hold on
            %            xm=mean(xo(conn),2);
            %            ym=mean(yo(conn),2);
            %             plot(xm,ym,'rx')
            %             eroot=prod(Nnodes)+reshape(1:prod(Nelems),Nelems);
            %              nroot=repmat(nroot(:)',4,1);
            %              eroot=repmat(eroot(:)',4,1);
            %             incn=[0,1;...
            %                 1,Nnodes(1)+1;...
            %                 Nnodes(1)+1,Nnodes(1);...
            %                 Nnodes(1),0];
            %             conn=[eroot(:),repmat(nroot(:),1,2)+repmat(incn,prod(Nelems),1),0*eroot(:)];
            %  elt=repmat(3,4*prod(Nelems),1);
            %  xo=[xo;xm];
            %  yo=[yo;ym];
            %     Nnodes=[length(xo),1,1];
            %     Nelems=[length(elt),1,1];
            % selected=[selected(:);repmat(0,length(xm),1)];
            
            
            incn=[1,Nnodes(1)+1,0;...
                Nnodes(1),0,Nnodes(1)+1;...
                0,1,Nnodes(1);...
                Nnodes(1)+1,Nnodes(1),1];
            Nelems2=ceil((Nelems)/2);
            incni=repmat(incn,Nelems2(1),1);
            incni=incni(1:(2*Nelems(1)),:);
            incni=[incni;circshift(incni,-2)];
            incnj=repmat(incni,Nelems2(2),1);
            incnj=incnj(1:2*prod(Nelems),:);
            nroot=reshape(1:prod(Nnodes),Nnodes);
            nroot=nroot(1:Nnodes(1)-1,:);
            nroot=nroot(:,1:Nnodes(2)-1);
            nroot=repmat(nroot(:)',2,1);
            conn=[repmat(nroot(:),1,3)+incnj,zeros(numel(nroot),1)];
            elt=repmat(3,2*prod(Nelems),1);
            Nelems=[length(elt),1];
            %triplot(conn,xo,yo)
            if strcmp(param0.analysis,'correlation')
                xo=xo+rfac*(2*rand(size(xo))-1);
                yo=yo+rfac*(2*rand(size(xo))-1);
            end
        end
        
        if check
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
            figure
            imagesc(im0')
            colormap('gray')
            hold on
            roi=param0.roi;
            plot(xo,yo,'bo')
            axis xy
            axis image
        end
        
        save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rflag','rint','selected','ns','ng','xo','yo','zo','Nnodes','Nelems','Smesh','Vmesh','mesh_size','conn','elt');
    end
    if iscale==1
        mesh_size_1=mesh_size;
    end
end
if isfield(param,'enrichment')&&(length(sizeim)==2)
    iscale=1;
    mesh_size=mesh_size_1;
    ng=ngi;ns=nsi;
    h=1;
    if isfield(param,'mesh_size_ratio')
        h=param.mesh_size_ratio;
    end
    if ~(h==1)
        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
        mesh_size_h=h*mesh_size;
        xh=round(0:mesh_size_h(1):sizeim(1)-1)+1;
        yh=round(0:mesh_size_h(2):sizeim(2)-1)+1;
        if length(sizeim)==3
            zh=round(0:mesh_size_h(3):sizeim(3)-1)+1;
        else
            zh=1;
        end
        
        
        xo=round(0:mesh_size(1):sizeim(1)-1)+1;
        yo=round(0:mesh_size(2):sizeim(2)-1)+1;
        if length(sizeim)==3
            zo=round(0:mesh_size(3):sizeim(3)-1)+1;
        else
            zo=1;
        end
        if xo(length(xo))>xh(length(xh));
            nn=find(xo==max(xh));
            xo=xo(1:nn);
        else
            nn=find(xh==max(xo));
            xh=xh(1:nn);
        end
        if yo(length(yo))>yh(length(yh));
            nn=find(yo==max(yh));
            yo=yo(1:nn);
        else
            nn=find(yh==max(yo));
            yh=yh(1:nn);
        end
        if zo(length(zo))>zh(length(zh));
            nn=find(zo==max(zh));
            zo=zo(1:nn);
        else
            nn=find(zh==max(zo));
            zh=zh(1:nn);
        end
        
        xh=xh-floor(mean(xo)-sizeim(1)/2);
        yh=yh-floor(mean(yo)-sizeim(2)/2);
        
        xo=xo-floor(mean(xo)-sizeim(1)/2);
        yo=yo-floor(mean(yo)-sizeim(2)/2);
        
        Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
        Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
        [yo,xo]=meshgrid(yo-0.5,xo-0.5);
        
        if isfield(param,'adapt_mesh')
            if param.adapt_mesh>0
                ic=param.adapt_mesh;
                load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'theta','dist');
                [px,py]=find(dist==min(dist(:)));
                px=mean(px);
                py=mean(py);
                if abs(abs(theta)-pi/2)<pi/4
                    theta=sign(theta)*abs(abs(theta)-pi/2);
                end
                theta=rem(theta,pi);
                uo=xo-px;vo=yo-py;
                xo=xo-vo*theta;
                yo=yo+uo*theta;
                xout=(xo<1.5)|(xo>sizeim(1)-0.5);
                found=find(sum(xout,2)>0);
                xo(found,:)=[];yo(found,:)=[];
                yout=(yo<1.5)|(yo>sizeim(2)-0.5);
                found=find(sum(yout,1)>0);
                xo(:,found)=[];yo(:,found)=[];
                %                 xo=max(1,min(sizeim(1),xo));
                %                 yo=max(1,min(sizeim(2),yo));
            end
            rflag=false;
        end
        Nnodes=[size(xo),1];
        Nelems=max(Nnodes-1,1);
        yo=yo(:);
        xo=xo(:);
        incn=[0,1,Nnodes(1)+1,Nnodes(1)];
        nroot=reshape(1:prod(Nnodes),Nnodes);
        nroot=nroot(1:Nnodes(1)-1,:);
        nroot=nroot(:,1:Nnodes(2)-1);
        conn=repmat(nroot(:),1,4)+repmat(incn,prod(Nelems),1);
        elt=repmat(4,prod(Nelems),1);
        if eelt==3&&iscale==1
            %            xm=mean(xo(conn),2);
            %            ym=mean(yo(conn),2);
            %             plot(xm,ym,'rx')
            %             eroot=prod(Nnodes)+reshape(1:prod(Nelems),Nelems);
            %              nroot=repmat(nroot(:)',4,1);
            %              eroot=repmat(eroot(:)',4,1);
            %             incn=[0,1;...
            %                 1,Nnodes(1)+1;...
            %                 Nnodes(1)+1,Nnodes(1);...
            %                 Nnodes(1),0];
            %             conn=[eroot(:),repmat(nroot(:),1,2)+repmat(incn,prod(Nelems),1),0*eroot(:)];
            %  elt=repmat(3,4*prod(Nelems),1);
            %  xo=[xo;xm];
            %  yo=[yo;ym];
            %     Nnodes=[length(xo),1,1];
            %     Nelems=[length(elt),1,1];
            % selected=[selected(:);repmat(0,length(xm),1)];
            
            
            incn=[1,Nnodes(1)+1,0;...
                Nnodes(1),0,Nnodes(1)+1;...
                0,1,Nnodes(1);...
                Nnodes(1)+1,Nnodes(1),1];
            Nelems2=ceil((Nelems)/2);
            incni=repmat(incn,Nelems2(1),1);
            incni=incni(1:(2*Nelems(1)),:);
            incni=[incni;circshift(incni,-2)];
            incnj=repmat(incni,Nelems2(2),1);
            incnj=incnj(1:2*prod(Nelems),:);
            nroot=reshape(1:prod(Nnodes),Nnodes);
            nroot=nroot(1:Nnodes(1)-1,:);
            nroot=nroot(:,1:Nnodes(2)-1);
            nroot=repmat(nroot(:)',2,1);
            conn=[repmat(nroot(:),1,3)+incnj,zeros(numel(nroot),1)];
            elt=repmat(3,2*prod(Nelems),1);
            Nelems=length(elt);
            xo=xo+rfac*(2*rand(size(xo))-1);
            yo=yo+rfac*(2*rand(size(xo))-1);
            figure
            plot(xo,yo,'bo')
            hold on
            triplot(conn,xo,yo)
        end
        
        
        
        
        
        selected=ones(Nnodes);
        selected(1,:)=0;
        selected(:,1)=0;
        selected(Nnodes(1),:)=0;
        selected(:,Nnodes(2))=0;
        
        %         if isfield(param,'adapt_mesh')
        %             if param.adapt_mesh>0
        %                 found=find((xo>=1)&(xo<=sizeim(1))&(yo>=1)&(yo<=sizeim(2)));
        %                 nconn=zeros(size(conn));
        %                 nelt=zeros(size(elt));
        %                 nselected=ones(length(found),1);
        %                 ielt=0;
        %                 for ie=1:length(elt)
        %                     inods=conn(ie,1:elt(ie));
        %                     xn=xo(inods);
        %                     yn=yo(inods);
        %                     founde=find((xn>=1)&(xn<=sizeim(1))&(yn>=1)&(yn<=sizeim(2)));
        %                     if length(founde)==elt(ie)
        %                         for in=1:elt(ie)
        %                             nconn(ielt+1,in)=find(found==inods(in));
        %                         end
        %                         nelt(ielt+1)=elt(ie);
        %                         ielt=ielt+1;
        %                     elseif any(~selected(inods))
        %                         for in=1:elt(ie)
        %                             if any(founde==in)
        %                                 nselected(find(found==inods(in)))=0;
        %                             end
        %                         end
        %
        %                     end
        %                 end
        %                 xo=xo(found);
        %                 yo=yo(found);
        %                 conn=nconn(1:ielt,:);
        %                 elt=nelt(1:ielt);
        %                 selected=nselected;
        %                 Nnodes=[length(xo),1,1];
        %                 Nelems=[size(conn,1),1,1];
        %             end
        %         end
        save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rflag','rint','selected','ns','ng','xo','yo','zo','Nnodes','Nelems','Smesh','Vmesh','mesh_size','conn','elt');
        
        xo=xh;
        yo=yh;
        zo=zh;
        
        
        
        
        Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
        Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
        [yo,xo]=meshgrid(yo-0.5,xo-0.5);
        if isfield(param,'adapt_mesh')
            if param.adapt_mesh>0
                ic=param.adapt_mesh;
                load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'theta','dist');
                [px,py]=find(dist==min(dist(:)));
                px=mean(px);
                py=mean(py);
                if abs(abs(theta)-pi/2)<pi/4
                    theta=sign(theta)*abs(abs(theta)-pi/2);
                end
                theta=rem(theta,pi);
                uo=xo-px;vo=yo-py;
                xo=xo-vo*theta;
                yo=yo+uo*theta;
                xout=(xo<1.5)|(xo>sizeim(1)-0.5);
                found=find(sum(xout,2)>0);
                xo(found,:)=[];yo(found,:)=[];
                yout=(yo<1.5)|(yo>sizeim(2)-0.5);
                found=find(sum(yout,1)>0);
                xo(:,found)=[];yo(:,found)=[];
                
                %                 xo=max(1,min(sizeim(1),xo));
                %                 yo=max(1,min(sizeim(2),yo));
            end
            rflag=false;
        end
        
        Nnodes=[size(xo),1];
        Nelems=max(Nnodes-1,1);
        yo=yo(:);
        xo=xo(:);
        
        incn=[0,1,Nnodes(1)+1,Nnodes(1)];
        nroot=reshape(1:prod(Nnodes),Nnodes);
        nroot=nroot(1:Nnodes(1)-1,:);
        nroot=nroot(:,1:Nnodes(2)-1);
        conn=repmat(nroot(:),1,4)+repmat(incn,prod(Nelems),1);
        elt=repmat(4,prod(Nelems),1);
        
        if eelt==3&&iscale==1
            %            xm=mean(xo(conn),2);
            %            ym=mean(yo(conn),2);
            %             plot(xm,ym,'rx')
            %             eroot=prod(Nnodes)+reshape(1:prod(Nelems),Nelems);
            %              nroot=repmat(nroot(:)',4,1);
            %              eroot=repmat(eroot(:)',4,1);
            %             incn=[0,1;...
            %                 1,Nnodes(1)+1;...
            %                 Nnodes(1)+1,Nnodes(1);...
            %                 Nnodes(1),0];
            %             conn=[eroot(:),repmat(nroot(:),1,2)+repmat(incn,prod(Nelems),1),0*eroot(:)];
            %  elt=repmat(3,4*prod(Nelems),1);
            %  xo=[xo;xm];
            %  yo=[yo;ym];
            %     Nnodes=[length(xo),1,1];
            %     Nelems=[length(elt),1,1];
            % selected=[selected(:);repmat(0,length(xm),1)];
            
            
            incn=[1,Nnodes(1)+1,0;...
                Nnodes(1),0,Nnodes(1)+1;...
                0,1,Nnodes(1);...
                Nnodes(1)+1,Nnodes(1),1];
            Nelems2=ceil((Nelems)/2);
            incni=repmat(incn,Nelems2(1),1);
            incni=incni(1:(2*Nelems(1)),:);
            incni=[incni;circshift(incni,-2)];
            incnj=repmat(incni,Nelems2(2),1);
            incnj=incnj(1:2*prod(Nelems),:);
            nroot=reshape(1:prod(Nnodes),Nnodes);
            nroot=nroot(1:Nnodes(1)-1,:);
            nroot=nroot(:,1:Nnodes(2)-1);
            nroot=repmat(nroot(:)',2,1);
            conn=[repmat(nroot(:),1,3)+incnj,zeros(numel(nroot),1)];
            elt=repmat(3,2*prod(Nelems),1);
            Nelems=length(elt);
            xo=xo+rfac*(2*rand(size(xo))-1);
            yo=yo+rfac*(2*rand(size(xo))-1);
        end
        
        %         if isfield(param,'adapt_mesh')
        %             if param.adapt_mesh>0
        %                 found=find((xo>=1)&(xo<=sizeim(1))&(yo>=1)&(yo<=sizeim(2)));
        %                 nconn=zeros(size(conn));
        %                 nelt=zeros(size(elt));
        %                 ielt=0;
        %                 for ie=1:length(elt)
        %                     inods=conn(ie,1:elt(ie));
        %                     xn=xo(inods);
        %                     yn=yo(inods);
        %                     founde=find((xn>=1)&(xn<=sizeim(1))&(yn>=1)&(yn<=sizeim(2)));
        %                     if length(founde)==elt(ie)
        %                         for in=1:elt(ie)
        %                             nconn(ielt+1,in)=find(found==inods(in));
        %                         end
        %                         nelt(ielt+1)=elt(ie);
        %                         ielt=ielt+1;
        %                     end
        %                 end
        %                 xo=xo(found);
        %                 yo=yo(found);
        %                 conn=nconn(1:ielt,:);
        %                 elt=nelt(1:ielt);
        %                 Nnodes=[length(xo),1,1];
        %                 Nelems=[size(conn,1),1,1];
        %             end
        %         end
        save(fullfile('TMP',sprintf('%d_meshx_%d',nmod,iscale-1)),'rflag','rint','ns','ng','xo','yo','zo','Nnodes','Nelems','Smesh','Vmesh','mesh_size','conn','elt');
    else
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rflag','rint','ns','ng','xo','yo','zo','Nnodes','Nelems','Smesh','Vmesh','mesh_size','conn','elt');
        save(fullfile('TMP',sprintf('%d_meshx_%d',nmod,iscale-1)),'rflag','rint','ns','ng','xo','yo','zo','Nnodes','Nelems','Smesh','Vmesh','mesh_size','conn','elt');
        
    end
end

end

