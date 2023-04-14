function VICMeshes(nmod)

load(fullfile('TMP','params'),'param');
tic;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
nscale=param.nscale;
type_nurbs=param.continuity;
load(fullfile('TMP',sprintf('sample0_%d',1-1)),'ls1','sizeim');
Nel=param.Nelems;
degree=param.degree;
for iscale=1:nscale
    load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'nband');
    if numel(ls1)>numel(nband)
        ls1b=ls1(nband);
    else
        ls1b=ls1(:);
    end
    spany=-min(ls1b)+max(ls1b);
    %    mesh_size=spany/max(1,Nel(1)/iscale);
    mesh_size=spany/ceil(max(1,Nel(1)/2^(iscale-1)));
    
    yo=(min(ls1b):mesh_size:max(ls1b));
    Nnodes=[length(yo),1,1];
    Nelems=max(Nnodes-1,1);
    vo=yo;
    
    
    switch type_nurbs
        % uniform open knot vector
        case 'open'
            for i=1:degree(1)
                vo=[vo(1),vo,vo(length(vo))];
            end
            % periodic knot vector
        case 'periodic'
            for i=1:degree(1)
                vo=[vo(1)-mesh_size(1),vo,vo(length(vo))+mesh_size(1)];
            end
            
        otherwise
            error('Invalid type of nurbs basis')
    end
    save(fullfile('TMP',sprintf('%d_vicmesh_%d',nmod,iscale-1)),'yo','vo','Nnodes','Nelems','mesh_size');
    vo=yo;
    for i=1:degree(1)
        vo=[vo(1),vo,vo(length(vo))];
    end
    save(fullfile('TMP',sprintf('%d_vicmeshs_%d',nmod,iscale-1)),'yo','vo','Nnodes','Nelems','mesh_size');
    if  length(sizeim)==3
        if iscale==1
            load(fullfile('TMP',sprintf('sample0_%d',1-1)),'ls2');
        end
        spanx=-min(ls2)+max(ls2);
        %        mesh_size=spanx/max(1,Nel(2)/iscale);
        mesh_size=spanx/max(1,Nel(2)/2^(iscale-1));
        
        xo=(min(ls2):mesh_size:max(ls2));
        Nnodes=[length(yo),length(xo),1];
        Nelems=max(Nnodes-1,1);
        uo=xo;
        
        switch type_nurbs
            % uniform open knot vector
            case 'open'
                for i=1:degree(2)
                    uo=[uo(1),uo,uo(length(uo))];
                end
                % periodic knot vector
            case 'periodic'
                for i=1:degree(2)
                    uo=[uo(1)-mesh_size,uo,uo(length(uo))+mesh_size];
                end
                
            otherwise
                error('Invalid type of nurbs basis')
        end
        save(fullfile('TMP',sprintf('%d_vicmesh_%d',nmod,iscale-1)),'xo','uo','Nnodes','Nelems','-append');
        
        
    end
    %    Nel=max(4,ceil(0.5*Nel));
end
disp(sprintf('Constructing vic mesh for model %d...%6.2fs',nmod,toc));
end

