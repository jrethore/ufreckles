function postproFEMVIC3D(filres)
iscale=1;
load(filres);
Ut=U;
nn=0;
if ~exist('models','var')
    models=model;
end
nmods=param.model_id;
for imod=1:length(nmods)
    nmod=nmods(imod);
    if iscell(models)
        model=models{imod};
    else
        model=models;
    end
    filresi=[filres,sprintf('-%02d',nmod)];
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nnodes','n')
    Un=Ut(nn+(1:prod(Nnodes)),:)+0.*model.transition_length;
    nn=nn+prod(Nnodes);
    LoadMat(nmod);
mu=0.5;lambda=0;
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambda');
        CreateGradBasisFunction3D(iscale,nmod);
    K=AssembleStiffnessMatrix(iscale,nmod);
     load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'rint','Nnodes','Nelems','xo','yo','zo','conn','elt','ng','selected')
     indio=find(~selected);
     indi=[indio(:);indio(:)+prod(Nnodes);indio(:)+2*prod(Nnodes)];
     Up=[n(:,1).*Un;n(:,2).*Un;n(:,3).*Un];
     indj=1:length(indi);
     C=sparse(indi,indj,1,size(K,1),length(indj));
     K=[K,C;C',sparse(size(C,2),size(C,2))];
     U=K\([zeros(3*prod(Nnodes),size(Up,2));Up]);
     U=U(1:3*prod(Nnodes),:);
     save(filresi,'param','model','nmod','U','rint','Nnodes','Nelems','xo','yo','zo','conn','elt','ng','selected')
     postproVTK3D(filresi,0,0)
     
end




end