function LoadMeshes(nmod,check)
if nargin<2,check=0;end
tic;
load(fullfile('TMP','params'),'param');
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
rflag=false;
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
load(fullfile('TMP',sprintf('sample0')),'sizeim');
if isfield(param,'mesh_file')&&(~isfield(param,'topography'))
    
    if isfield(param,'mesh_size')
        RegularMeshes(nmod,check);
        if length(sizeim)==3
            UnstructuredMesh3D(nmod,1);
        else
            UnstructuredMesh(nmod,1);
        end
        AdaptRegularMeshes(nmod);
    else
        if iscell(param.mesh_file)
            for iscale=1:nscale
                if length(sizeim)==3
                    UnstructuredMesh3D(nmod,iscale);
                else
                    UnstructuredMesh(nmod,iscale);
                end
            end
        else
            
            if length(sizeim)==3
                UnstructuredMesh3D(nmod,1);
                load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,1-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
                for iscale=2:nscale
                    xo=(xo-0.5)*0.5+0.5;
                    yo=(yo-0.5)*0.5+0.5;
                    zo=(zo-0.5)*0.5+0.5;
                    save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
                end
            else
                UnstructuredMesh(nmod,1);
                load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,1-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
                for iscale=2:nscale
                    xo=(xo-0.5)*0.5+0.5;
                    yo=(yo-0.5)*0.5+0.5;
                    save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
                end
            end
        end
    end
else
    RegularMeshes(nmod,check);
end

if isfield(param,'topography')
    LoadTopography(nmod);
end
if eelt==33
    GenerateBezierTriangles(nmod);
end
disp(sprintf('Constructing meshes for model %d...%6.2fs',nmod,toc));

end

