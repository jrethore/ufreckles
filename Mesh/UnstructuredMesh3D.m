function UnstructuredMesh3D(nmod,iscale)
check=0;

load(fullfile('TMP','params'),'param');
roi=param.roi;
    tic;
load(fullfile('TMP','sample0'),'sizeim');
sizeim0=sizeim;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
rint=false;
if isfield(param,'reduced_integration')
    rint=param.reduced_integration;
end
meshfile=param.mesh_file;
if iscale>1
    assert(iscell(meshfile));
end
rflag=false;
if iscell(meshfile)
    meshfile=meshfile{iscale};
end 
pscale=2^(iscale-1);
[xo,yo,zo,conn,elt,selected]=ReadVTK3D(meshfile);
if isfield(param,'gluing_parameters')
    if ~isempty(param.gluing_parameters)
   [xo,yo,zo]=GlueMesh(param.gluing_parameters,xo,yo,zo);
    end
end
    if check
       figure 
       plot3(xo,yo,zo,'x')
       hold on
%       line(roi([1,2,1,2]),roi([4,4,3,3]),roi([5,5,5,5]),'r-')
%       line(roi([1,2,1,2]),roi([4,4,3,3]),roi([6,6,6,6]),'r-')
    end
xo=xo-roi(1)+1;
yo=yo-roi(3)+1;
zo=zo-roi(5)+1;

load(fullfile('TMP',sprintf('sample0_%d',1-1)),'sizeim');

if any(xo<1)||any(xo>sizeim(1))||any(yo<1)||any(yo>sizeim(2))||any(zo<1)||any(zo>sizeim(3))
    display('WARNING, NODES OUT OF THE ROI')
    display(sprintf('IMAGE SIZE: %d %d %d',sizeim0));
    display(sprintf('ROI: %d %d %d %d %d %d',roi));
    display(sprintf('MESH AREA: %g %g %g %g %g %g',min(xo)+roi(1)-1,max(xo)+roi(1)-1,min(yo)+roi(3)-1,max(yo)+roi(3)-1,min(zo)+roi(5)-1,max(zo)+roi(5)-1));
    display('THE MESH IS SHRINKED INSIDE THE ROI')
    conng=conn;
    conng(conn==0)=numel(xo)+1;
    xg=[xo;0];
    xg=xg(conng);
    keep=sum((xg<1)|(xg>sizeim(1)),2);
    xg=[yo;0];
    xg=xg(conng);
    keep=keep+sum((xg<1)|(xg>sizeim(2)),2);
    xg=[zo;0];
    xg=xg(conng);
    keep=keep+sum((xg<1)|(xg>sizeim(3)),2);
    conn(keep>0,:)=[];
    elt(keep>0,:)=[];
    
    [pind,~,j1]=unique(conn);
    conn=reshape(j1,size(conn));
    xo=xo(pind);
    yo=yo(pind);
    zo=zo(pind);
    
end
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
if iscale>1
xo=round((xo-0.5)/pscale)+0.5;
yo=round((yo-0.5)/pscale)+0.5;
zo=round((zo-0.5)/pscale)+0.5;
end
    xo=max(xo,1);
xo=min(xo,sizeim(1));
yo=max(yo,1);
yo=min(yo,sizeim(2));
zo=max(zo,1);
zo=min(zo,sizeim(3));


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
    disp(sprintf('Loading mesh for model %d...%6.2fs',nmod,toc));
    
    

end

