function AdaptRegularMeshes(nmod)
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
nscale=param.nscale;
load(fullfile('TMP',sprintf('sample0')),'sizeim');
dim=numel(sizeim);
switch dim
    case 2
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,1-1)),'xo','yo');
        zi=0;
    case 3
        load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,1-1)),'xo','yo','zo');
        zi=zo;
end
xi=xo;yi=yo;
for iscale=2:nscale
    pscale=2^(iscale-1);
    switch dim
        case 2
            load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','conn','elt','Nnodes','Nelems');
            zo=0*xo;
        case 3
            load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','conn','elt','Nnodes','Nelems');
    end
    xo=(xo-0.5)*pscale+0.5;
    yo=(yo-0.5)*pscale+0.5;
    zo=(zo-0.5)*pscale+0.5;
    
    keep=zeros(prod(Nnodes),1);
    for ie=1:prod(Nelems)
        inods=conn(ie,1:elt(ie));
        switch dim
            case 2
                keepe=any(inpolygon(xi,yi,xo(inods([1:numel(inods),1])),yo(inods([1:numel(inods),1]))));
            case 3
                xn=xo(inods);yn=yo(inods);zn=zo(inods);
                inbox=(xi>min(xn))&(xi<max(xn))&(yi>min(yn))&(yi<max(yn))&(zi>min(zn))&(zi<max(zn));
                keepe=any(inbox);
        end
        if keepe,keep(inods)=1;end
    end
    keep=find(keep);
    ielt=GetEltsFromNodes(conn,elt,keep,1);
    if length(ielt)>0
    elt=elt(ielt);
    conn=conn(ielt,:);
    xo=xo(keep);
    yo=yo(keep);
    zo=zo(keep);
    newids=zeros(1,prod(Nnodes)+1);
    newids(keep)=1:length(keep);
    conn(conn==0)=prod(Nnodes)+1;
    conn=newids(conn);
    Nnodes=[numel(xo),1,1];
    Nelems=[numel(elt),1,1];
    
    xo=(xo-0.5)/pscale+0.5;
    yo=(yo-0.5)/pscale+0.5;
    zo=(zo-0.5)/pscale+0.5;
    switch dim
        case 2
            save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','conn','elt','Nnodes','Nelems','-append');
        case 3
            save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','conn','elt','Nnodes','Nelems','-append');
    end
    end
end
end