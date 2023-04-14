function NURBSMeshes(nmod)
rint=false;
rflag=0;
load(fullfile('TMP','params'),'param');
param0=param;
tic;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if strcmp(param.basis,'nurbs-beam')
    BEAMMesh(nmod);
    return
end
nscale=param.nscale;
if isfield(param,'continuity')
    type_nurbs=param.continuity;
else
    type_nurbs='open';
end

for iscale=1:nscale
    load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
    
    if ~((iscale==1)&&isfield(param,'nurbs_mesh'))
        if iscale==1
            degree=param.degree;
            mesh_size=param.mesh_size;
            rflag=false;
        else
            mesh_size=round(param.mesh_size(:)./param.degree(:));
            degree=ones(1,3);
            rflag=true;
        end
        xo=(1:mesh_size(1):sizeim(1)-1)+1;
        yo=(1:mesh_size(2):sizeim(2)-1)+1;
        if (length(xo)==1)
            mesh_size(1)=sizeim(1)-1;
            xo=(0:mesh_size(1):sizeim(1)-1)+1;
        end
        if (length(yo)==1)
            mesh_size(2)=sizeim(2)-1;
            yo=(0:mesh_size(2):sizeim(2)-1)+1;
        end
        
        xo=xo-floor(mean(xo)-sizeim(1)/2);
        yo=yo-floor(mean(yo)-sizeim(2)/2);
        if length(sizeim)==3
            zo=(0:mesh_size(3):sizeim(3)-1)+1;
            if (length(zo)==1)
                mesh_size(3)=sizeim(3)-1;
                zo=(0:mesh_size(3):sizeim(3)-1)+1;
            end
        else
            zo=1;
        end
        Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
        Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
        Nnodes=[length(xo),length(yo),length(zo)];
        Nelems=max(Nnodes-1,1);
        Nbselems=Nelems;
        yo=yo-0.5;
        xo=xo-0.5;
        uo=xo-xo(1)-0.5;
        vo=yo-yo(1)-0.5;
        
        Px=[min(xo),max(xo)];
        Py=[min(yo),max(yo)];
        ud=uo([1,1,length(uo),length(uo)]);
        vd=vo([1,1,length(vo),length(vo)]);
        if degree(1)>1
            [Px,ud] = bspdegelev(1,Px,ud,degree(1)-1);
        end
        if degree(2)>1
            [Py,vd] = bspdegelev(1,Py,vd,degree(2)-1);
        end
        if length(uo)>2
            [Px,unew] = bspkntins(degree(1),Px,ud,uo(2:length(uo)-1));
        end
        if length(vo)>2
            [Py,vnew] = bspkntins(degree(2),Py,vd,vo(2:length(vo)-1));
        end
        if length(sizeim)==2
            wo=zo;
            Nbsnodes=[length(Px),length(Py),1];
            
            [Py,Px]=meshgrid(Py,Px);
            if degree(1)>1
                for ir=1:degree(1)-1
                    xo=sort([xo,0.5*diff(xo)+xo(1:length(xo)-1)]);
                end
            end
            if degree(2)>1
                for ir=1:degree(2)-1
                    yo=sort([yo,0.5*diff(yo)+yo(1:length(yo)-1)]);
                end
            end
            Nnodes=[length(xo),length(yo),length(zo)];
            Nelems=max(Nnodes-1,1);
            [yo,xo]=meshgrid(yo,xo);
            yo=yo(:);
            xo=xo(:);
            incn=[0,1,Nnodes(1)+1,Nnodes(1)];
            nroot=reshape(1:prod(Nnodes),Nnodes);
            nroot=nroot(1:max(1,Nnodes(1)-1),:);
            nroot=nroot(:,1:max(1,Nnodes(2)-1));
            conn=repmat(nroot(:),1,4)+repmat(incn,prod(Nelems),1);
            elt=repmat(4,prod(Nelems),1);
            ng=0;
            ns=ones(1,2);
            elt=repmat(4,prod(Nelems),1);
            selected=ones(Nnodes);
            selected(1,:)=0;
            selected(:,1)=0;
            selected(Nnodes(1),:)=0;
            selected(:,Nnodes(2))=0;
            %             if max(degree)>1
            %                 [elt,conn,xo,yo,Nnodes,Nelems,selected]=BuiltRefinedConnectivity(elt,conn,xo,yo,selected(:),max(degree));
            %             end
            ui=xo-xo(1)-0.5;
            vi=yo-yo(1)-0.5;
        else
            zo=zo-0.5;
            wo=zo-zo(1)-0.5;
            Pz=[min(zo),max(zo)];
            wd=wo([1,1,length(wo),length(wo)]);
            if degree(3)>1
                [Pz,wd] = bspdegelev(1,Pz,wd,degree(3)-1);
            end
            if length(wo)>2
                [Pz,wnew] = bspkntins(degree(3),Pz,wd,wo(2:length(wo)-1));
            end
            Nbsnodes=[length(Px),length(Py),length(Pz)];
            [Py,Px,Pz]=meshgrid(Py,Px,Pz);
            
            if degree(1)>1
                for ir=1:degree(1)-1
                    xo=sort([xo,0.5*diff(xo)+xo(1:length(xo)-1)]);
                end
            end
            if degree(2)>1
                for ir=1:degree(2)-1
                    yo=sort([yo,0.5*diff(yo)+yo(1:length(yo)-1)]);
                end
            end
            if degree(3)>1
                for ir=1:degree(3)-1
                    zo=sort([zo,0.5*diff(zo)+zo(1:length(zo)-1)]);
                end
            end
            Nnodes=[length(xo),length(yo),length(zo)];
            Nelems=max(Nnodes-1,1);
            
            
            [yo,xo,zo]=meshgrid(yo,xo,zo);
            yo=yo(:);
            xo=xo(:);
            zo=zo(:);
            
            
            incn=[0,1,Nnodes(1)+1,Nnodes(1)];
            incn=[incn,incn+prod(Nnodes(1:2))];
            nroot=reshape(1:prod(Nnodes),Nnodes);
            nroot=nroot(1:max(1,Nnodes(1)-1),:,:);
            nroot=nroot(:,1:max(1,Nnodes(2)-1),:);
            nroot=nroot(:,:,1:max(1,Nnodes(3)-1));
            conn=repmat(nroot(:),1,8)+repmat(incn,prod(Nelems),1);
            ng=0;
            ns=ones(1,3);
            elt=repmat(8,prod(Nelems),1);
            selected=ones(Nnodes);
            selected(1,:,:)=0;
            selected(:,1,:)=0;
            selected(:,:,1)=0;
            selected(Nnodes(1),:,:)=0;
            selected(:,Nnodes(2),:)=0;
            selected(:,:,Nnodes(3))=0;
            ui=xo-xo(1)-0.5;
            vi=yo-yo(1)-0.5;
            wi=zo-zo(1)-0.5;
        end
        
        switch type_nurbs
            % uniform open knot vector
            case 'open'
                for i=1:degree(1)
                    uo=[uo(1),uo,uo(length(uo))];
                end
                for i=1:degree(2)
                    vo=[vo(1),vo,vo(length(vo))];
                end
                if length(sizeim)==3
                    for i=1:degree(3)
                        wo=[wo(1),wo,wo(length(wo))];
                    end
                end
                % periodic knot vector
            case 'periodic'
                for i=1:degree(1)
                    uo=[uo(1)-mesh_size(1),uo,uo(length(uo))+mesh_size(1)];
                end
                for i=1:degree(2)
                    vo=[vo(1)-mesh_size(2),vo,vo(length(vo))+mesh_size(2)];
                end
                if length(sizeim)==3
                    for i=1:degree(3)
                        wo=[wo(1)-mesh_size(3),wo,wo(length(wo))+mesh_size(3)];
                    end
                end
            case 'c0'
                j=1;
                for i=1:length(uo)
                    uuo(j)=uo(i);
                    for k=1:(degree(1)-1)
                        uuo(j+k)=uo(i);
                    end
                    j=j+degree(1);
                end
                uuo=[uo(1),uuo,uo(length(uo))];
                uo=uuo;
                j=1;
                for i=1:length(vo)
                    vvo(j)=vo(i);
                    for k=1:(degree(2)-1)
                        vvo(j+k)=vo(i);
                    end
                    j=j+degree(2);
                end
                vvo=[vo(1),vvo,vo(length(vo))];
                vo=vvo;
                
            otherwise
                error('Invalid type of nurbs basis')
        end
        
        
    else
        roi=param0.roi;
        if length(roi)>4
            mesh_size=param.mesh_size;
            load(param.nurbs_mesh)
            Px=Px-roi(1)+1;
            Py=Py-roi(3)+1;
            Pz=Pz-roi(5)+1;
        else
            load(param.nurbs_mesh,'Px','Py','uo','vo','degree')
            pdegree=param.degree;
            
            
            if isfield(param,'gluing_parameters')
                [Px,Py]=GlueMesh(param.gluing_parameters,Px,Py);
            end
            uo=uo-uo(1);
            vo=vo-vo(1);
            du=diff(uo);
            ne=sum(du>0);
            du=sum(du)/ne;
            uo=uo/du*floor(sizeim(1)/ne)-0.5;
            du=diff(vo);
            ne=sum(du>0);
            du=sum(du)/ne;
            vo=vo/du*floor(sizeim(2)/ne)-0.5;
            
            Px=Px-roi(1)+1;
            Py=Py-roi(3)+1;
            
            if pdegree(1)>degree(1)
                Pxn=[];Pyn=[];
                for jj=1:size(Px,2)
                    [Pxi,ud] = bspdegelev(degree(1),Px(:,jj)',uo,pdegree(1)-degree(1));
                    Pxn=[Pxn,Pxi'];
                    [Pyi,ud] = bspdegelev(degree(1),Py(:,jj)',uo,pdegree(1)-degree(1));
                    Pyn=[Pyn,Pyi'];
                end
                Px=Pxn;
                Py=Pyn;
            else
                ud=uo;
            end
            if pdegree(2)>degree(2)
                Pxn=[];Pyn=[];
                for ii=1:size(Px,1)
                    [Pxi,vd] = bspdegelev(degree(2),Px(ii,:),vo,pdegree(2)-degree(2));
                    Pxn=[Pxn;Pxi];
                    [Pyi,vd] = bspdegelev(degree(2),Py(ii,:),vo,pdegree(2)-degree(2));
                    Pyn=[Pyn;Pyi];
                end
                Px=Pxn;
                Py=Pyn;
            else
                vd=vo;
            end
            
            if isfield(param,'refinement')
                
                if param.refinement(1)
                    urnew=[];
                    for ie=1:(length(ud)-2*pdegree(1)-1)
                        
                        un=ud(ie+pdegree(1)+(0:1));
                        unnew=[];
                        for ir=1:param.refinement(1)
                            unew=[];
                            for in=1:length(un)-1
                                unew=[unew,mean(un(in+(0:1)))];
                            end
                            un=sort([un,unew]);
                            unnew=sort([unnew,unew]);
                        end
                        urnew=[urnew,unnew];
                        
                    end
                    unew=sort(urnew);
                    Pxn=[];Pyn=[];
                    for jj=1:size(Px,2)
                        [Pxi,utmp] = bspkntins(pdegree(1),Px(:,jj)',ud,unew);
                        Pxn=[Pxn,Pxi'];
                        [Pyi,utmp] = bspkntins(pdegree(1),Py(:,jj)',ud,unew);
                        Pyn=[Pyn,Pyi'];
                    end
                    Px=Pxn;
                    Py=Pyn;
                    ud=utmp;
                end
                
                if param.refinement(2)
                    urnew=[];
                    for ie=1:(length(vd)-2*pdegree(2)-1)
                        
                        un=vd(ie+pdegree(2)+(0:1));
                        unnew=[];
                        for ir=1:param.refinement(2)
                            unew=[];
                            for in=1:length(un)-1
                                unew=[unew,mean(un(in+(0:1)))];
                            end
                            un=sort([un,unew]);
                            unnew=sort([unnew,unew]);
                        end
                        urnew=[urnew,unnew];
                        
                    end
                    unew=sort(urnew);
                    Pxn=[];Pyn=[];
                    for ii=1:size(Px,1)
                        [Pxi,vtmp] = bspkntins(pdegree(2),Px(ii,:),vd,unew);
                        Pxn=[Pxn;Pxi];
                        [Pyi,vtmp] = bspkntins(pdegree(2),Py(ii,:),vd,unew);
                        Pyn=[Pyn;Pyi];
                    end
                    Px=Pxn;
                    Py=Pyn;
                    vd=vtmp;
                end
            end
            uo=ud;
            vo=vd;
            Nbsnodes=[size(Px),1];
            xo=uo((pdegree(1)+1):(length(uo)-pdegree(1)))+1;
            yo=vo((pdegree(2)+1):(length(vo)-pdegree(2)))+1;
            
            Nnodes=[length(xo),length(yo),1];
            Nbselems=max(Nnodes-1,1);
            
            if pdegree(1)>1
                for ir=1:pdegree(1)
                    xo=sort([xo,0.5*diff(xo)+xo(1:length(xo)-1)]);
                end
            end
            if pdegree(2)>1
                for ir=1:pdegree(2)
                    yo=sort([yo,0.5*diff(yo)+yo(1:length(yo)-1)]);
                end
            end
            
            [yo,xo]=meshgrid(yo,xo);
            Nnodes=[size(xo),1];
            yo=yo(:);
            xo=xo(:);
            Nelems=max(Nnodes-1,1);
            incn=[0,1,Nnodes(1)+1,Nnodes(1)];
            nroot=reshape(1:prod(Nnodes),Nnodes);
            nroot=nroot(1:max(1,Nnodes(1)-1),:);
            nroot=nroot(:,1:max(1,Nnodes(2)-1));
            conn=repmat(nroot(:),1,4)+repmat(incn,prod(Nelems),1);
            elt=repmat(4,prod(Nelems),1);
            ng=0;
            ns=ones(1,2);
            elt=repmat(4,prod(Nelems),1);
            selected=ones(Nnodes);
            selected(1,:)=0;
            selected(:,1)=0;
            selected(Nnodes(1),:)=0;
            selected(:,Nnodes(2))=0;
            degree=pdegree;
            %         if max(degree)>1
            %             [elt,conn,xo,yo,Nnodes,Nelems,selected]=BuiltRefinedConnectivity(elt,conn,xo,yo,selected(:),max(degree));
            %         end
             ui=xo-xo(1)-0.5;
             vi=yo-yo(1)-0.5;
            zo=1;
            wo=zo;
            Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
            Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
        end
        
        
    end
    
    if length(sizeim)==2
        save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rflag','rint','ns','ng','conn','elt','ui','vi','zo','uo','vo','wo','Px','Py','Nbsnodes','Nnodes','Nbselems','Nelems','Smesh','Vmesh','selected','degree','-v7.3');
        [phio,xo,yo]=CreateNURBSBasis(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),degree,'nodes');
        save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','-append');
        save(fullfile('TMP',sprintf('%d_phio_%d',nmod,iscale-1)),'phio','-v7.3');
        
        if (iscale==1)&&isfield(param,'nurbs_mesh')
            
            load(fullfile('TMP','sample0'),'im0');
%             figure
%             
%             colormap(gray)
%             imagesc(im0')
%             hold on;
%             %             rect1=rectangle;
%             %             set(rect1,'position',[roi(1),roi(3),roi(2)-roi(1),roi(4)-roi(3)],...
%             %                 'EdgeColor','yellow',...
%             %                 'LineStyle','--',...
%             %                 'LineWidth',2)
%             %            plot((xo-1)+roi(1),(yo-1)+roi(3),'b+','LineWidth',0.5,'MarkerSize',10);
%             plot((xo-1)+roi(1),(yo-1)+roi(3),'b.','LineWidth',1.5,'MarkerSize',10);
%             axis xy
%             axis image
%             
%             title(param.nurbs_mesh,'FontSize',24);
            
            
            figure
            
            colormap(gray)
            imagesc(im0')
            hold on;
            
            vtmp=vo;
            utmp=uo;
            vo=utmp;
            yo=vo((degree(1)+1):(length(vo)-degree(1)))+1;
            yi=(min(yo)):(max(yo));
            mesh_file=fullfile('TMP',sprintf('export_mesh'));
            save(mesh_file,'yo','vo');
            [phi]=CreateNURBSBasis0D(mesh_file,degree(1),yi);
            for iy=1:size(Px,2)
                Pxi=Px(:,iy);
                Pyi=Py(:,iy);
                xi=phi*Pxi;
                yi=phi*Pyi;
                plot((xi-1)+roi(1),(yi-1)+roi(3),'-','LineWidth',1.,'MarkerSize',10);
                
            end
            vo=vtmp;
            yo=vo((degree(2)+1):(length(vo)-degree(2)))+1;
            yi=(min(yo)):(max(yo));
            mesh_file=fullfile('TMP',sprintf('export_mesh'));
            save(mesh_file,'yo','vo');
            [phi]=CreateNURBSBasis0D(mesh_file,degree(2),yi);
            for ix=1:size(Px,1)
                Pxi=Px(ix,:)';
                Pyi=Py(ix,:)';
                xi=phi*Pxi;
                yi=phi*Pyi;
                plot((xi-1)+roi(1),(yi-1)+roi(3),'-','LineWidth',1.,'MarkerSize',10);
                
            end
            axis xy
            axis image
            axis off;
            
            title(param.nurbs_mesh,'FontSize',24);
            
        end
        
    else
        save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'rflag','rint','ns','ng','conn','elt','ui','vi','wi','uo','vo','wo','Px','Py','Pz','Nbsnodes','Nnodes','Nbselems','Nelems','Smesh','Vmesh','mesh_size','selected','degree','-v7.3');
        [phio,xo,yo,zo]=CreateNURBSBasis3D(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),degree,'nodes');
        save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','-append');
        save(fullfile('TMP',sprintf('%d_3d_phio_%d',nmod,iscale-1)),'phio','-v7.3');
    end
    
    
    if isfield(param,'topography')
        LoadTopography(nmod);
    end
    
    
end
    disp(sprintf('Constructing meshes for model %d...%6.2fs',nmod,toc));

