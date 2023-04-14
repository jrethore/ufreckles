function run_fem_job_3D(param,model)

nmod=0;
LoadParameters(param);
LoadParameters(model,nmod);
ReferenceImage(nmod);
LoadMeshes(nmod);
LoadMask(nmod);
nscale=model.nscale;
U=[];
for iscale=nscale:-1:1
    Uini=InitializeSolution3D(U,nmod,iscale);
    [U]=Solve3D(Uini,nmod,iscale);
end
unix(['cp ' ,fullfile('TMP','0_error_0.mat'),' ', strrep(param.result_file,'.res','-error.res')]);
load(fullfile('TMP','0_3d_mesh_0'),'Nnodes','Nelems','xo','yo','zo','conn','elt','ng','rint','Smesh','ns');
save(param.result_file,'U','Nnodes','Nelems','xo','yo','zo','param','model','nmod','conn','elt','ng','rint','Smesh','ns');
postproVTK3D(param.result_file,0,1);
end