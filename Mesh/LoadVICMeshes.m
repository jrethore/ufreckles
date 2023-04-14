function LoadVICMeshes(nmods,check)
if nargin<2,check=0;end
iscale=1;
for im=1:length(nmods)
    nmod=nmods(im);
    disp(sprintf('Reading meshes for solid %d...',nmod));
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
    lc=param.transition_length;
    rint=false;
    if isfield(param,'reduced_integration')
        rint=param.reduced_integration;
    end
    rflag=false;
    meshfile=param.mesh_file;
    [xo,yo,zo,conn,elt,selected]=ReadVTK3D(meshfile);
    
    if isfield(param,'gluing_parameters')
        if ~isempty(param.gluing_parameters)
            [xo,yo,zo]=GlueMesh(nmod,xo,yo,zo);
        end
    end
    Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
    Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
    Nnodes=[length(xo),1,1];
    Nelems=[size(conn,1),1,1];
    ng=0;
    if isfield(param,'nb_gauss_points')
        ng=param.nb_gauss_points;
    end
    ns=ones(1,3);
    if isfield(param,'nb_sub_cells')
        ns=param.nb_sub_cells;
    end
    save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
    getSurf(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),xo,yo,zo,conn,elt,find(~selected),0)
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','Nnodes','Nelems','conn','elt');
    n=zeros(prod(Nnodes),3);
    for in=1:prod(Nnodes)
        nconn=find(conn==in);
        nn=0;
        for ic=1:length(nconn)
            [ie,ien]=ind2sub(size(conn),nconn(ic));
            inods=conn(ie,1:elt(ie));
            xn=xo(inods(1+mod(ien+[-1,0,1]-1,elt(ie))));
            yn=yo(inods(1+mod(ien+[-1,0,1]-1,elt(ie))));
            zn=zo(inods(1+mod(ien+[-1,0,1]-1,elt(ie))));
            nx=(yn(1)-yn(2))*(zn(3)-zn(2))-(zn(1)-zn(2))*(yn(3)-yn(2));
            ny=(zn(1)-zn(2))*(xn(3)-xn(2))-(xn(1)-xn(2))*(zn(3)-zn(2));
            nz=(xn(1)-xn(2))*(yn(3)-yn(2))-(yn(1)-yn(2))*(xn(3)-xn(2));
            nn=nn+[nx,ny,nz];
        end
        nn=nn/length(nconn);
        nn=nn/norm(nn);
        n(in,:)=-nn;
    end
    selected=[zeros(length(xo),1);ones(length(xo),1)];
    xo=[xo-0.5*lc*n(:,1);xo+0.5*lc*n(:,1)];
    yo=[yo-0.5*lc*n(:,2);yo+0.5*lc*n(:,2)];
    zo=[zo-0.5*lc*n(:,3);zo+0.5*lc*n(:,3)];
    conn=[conn,zeros(length(elt),8-size(conn,2))];
    for ie=1:length(elt)
        conn(ie,elt(ie)+(1:elt(ie)))=conn(ie,(1:elt(ie)))+prod(Nnodes);
        elt(ie)=2*elt(ie);
    end
    Nnodes=[length(xo),1,1];
    save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'n','-append');
    save(fullfile('TMP',sprintf('%d_3d_vicmesh_%d',nmod,iscale-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','ng','ns','selected');
    if check
        writeVTKmesh25D(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)));
        writeVTKmesh3D(fullfile('TMP',sprintf('%d_3d_vicmesh_%d',nmod,iscale-1)));
    end
end
end