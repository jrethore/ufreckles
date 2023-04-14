function BEAMMesh(nmod)
rint=false;

load(fullfile('TMP','params'),'param');
param0=param;
tic;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
assert(strcmp(param.basis,'nurbs-beam'))
nscale=param.nscale;
if isfield(param,'continuity')
    type_nurbs=param.continuity;
else
    type_nurbs='open';
end

for iscale=1:nscale
    load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
    
        if iscale==1
            degree=param.degree;
            mesh_size=param.mesh_size;
            degree(2)=degree(1);
        else
            mesh_size=round(param.mesh_size./param.degree);
            degree=ones(1,3);
        end
        xo=(0:mesh_size(1):sizeim(1)-1)+1;
        yo=(0:mesh_size(2):sizeim(2)-1)+1;
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
        uo=xo;
        vo=yo;
        
        if length(sizeim)==2
            wo=zo;
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
        else
            zo=zo-0.5;
            wo=zo;
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
        
    
    if length(sizeim)==2
        save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','rint','ns','ng','conn','elt','zo','uo','vo','wo','Nnodes','Nelems','Smesh','Vmesh','selected');        
    else
        save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','rint','ns','ng','conn','elt','uo','vo','wo','Nnodes','Nelems','Smesh','Vmesh','mesh_size','selected','degree');
    end
    
end


disp(sprintf('Constructing meshes for model %d...%6.2fs',nmod,toc));
end

