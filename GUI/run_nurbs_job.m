function run_nurbs_job(param,model,Up)
nmod=0;
if nargin<3,Up=[];end
clear functions
param.onflight=0;
model.basis='nurbs';
model.continuity='open';
LoadParameters(param);
LoadParameters(model,nmod);
ReferenceImage(nmod);
LoadMask(nmod);
NURBSMeshes(nmod);
LoadMat(nmod);

nscale=model.nscale;
filres=param.result_file;
%%
load(fullfile('TMP','sample0_0.mat'),'sizeim')
if prod(sizeim)>4e6
    param.initialization='bilin+';
%    param.initialization='rbt_elem';
    LoadParameters(param);
end
load(fullfile('TMP','sample0_0.mat'),'sizeim')
%%
U1=[];
if prod(sizeim)<4e6&&~param.normalize_grey_level
    
    for iscale=nscale:-1:1
        disp(sprintf('Pre-processing scale %d...',iscale));
        CreateBasisFunction(iscale,nmod);
        ComputeGradFPhi(iscale,nmod);
        CreateGradBasisFunction(iscale,nmod);
        AssembleCorrelationOperator(iscale,nmod);

            Uini=InitializeSolution(U1,iscale,nmod);
            [U1]=Solve(Uini,iscale,nmod);        
    end
    
else
    for iscale=nscale:-1:1
            Uini=InitializeSolution(U1,iscale,nmod);
            [U1]=Solve2D(Uini,iscale,nmod);
    end
end
%%
filreso=strrep(filres,'.res','');
UP=U1;
U=CPToNodes(U1,nmod,1);
copyfile(fullfile('TMP','0_error_0.mat'), [filreso,'-error.res']);
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
model=param;
load(fullfile('TMP','params'),'param');
load(fullfile('TMP','0_mesh_0'),'Px','Py','ui','vi','uo','vo','Nbsnodes','Nbselems','Nnodes','Nelems','xo','yo','conn','elt','ng','rflag','rint');
%ExportImageToVTK(filres)
delete([fullfile('VTK',[filreso,'-error']),'-0*.vtk']);
delete(fullfile('VTK',['camr*','-error.vtk']));

mesh.xo=xo;
mesh.yo=yo;
mesh.conn=conn;
mesh.elt=elt;
LoadMeshes(nmod)
model.basis='fem';
LoadParameters(model,nmod);
RefineMesh(nmod,ceil(mean(model.degree-1)));
load(fullfile('TMP','0_mesh_0'),'xo','yo');
coords.xi=xo+0.5;
coords.yi=yo+0.5;
[U1]=interpMesh(mesh,reshape(U,prod(Nnodes),2*size(U,2)*size(U,3)),coords);
load(fullfile('TMP','0_mesh_0'),'Nnodes','Nelems','xo','yo','conn','elt','ng','rflag','rint');
U1=reshape(U1,[2*prod(Nnodes),size(U,2),size(U,3)]);
U=U1;
model.mesh_file=strrep(filres,'.res','-res.vtk');
save(filres,'U','UP','Px','Py','ui','vi','uo','vo','Nbsnodes','Nbselems','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
if isfield(param,'calibration_data')
    [U,Xo,Yo,Zo]=Extract3DFields(U1,nmod);
    CreateGradBasisFunction25D(iscale,nmod);
    
    save(filres,'U','U1','Xo','Yo','Zo','-append');
    writeVTKmesh25D(filres);
    postproVTK25D(filres,1);
else
postproVTK([filres,''],0,0);
end
param.result_file=strrep(filres,'.res','-res.res');
copyfile(filres,strrep(filres,'.res','-res.res'));
xo=xo+param.roi(1)-1;
yo=yo+param.roi(3)-1;
save(strrep(filres,'.res','-res.res'),'xo','yo','-append');
writeVTKmesh(strrep(filres,'.res','-res.res'));
delete(strrep(filres,'.res','-res.res'));
end