function UnstructuredMesh(nmod,iscale,dec)
if nargin<3,dec=zeros(2,1);end
check=0;
tic;
rflag=false;
load(fullfile('TMP','params'),'param');
roi=param.roi;
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
meshfile=param.mesh_file;
if iscale>1
    assert(iscell(meshfile));
end
if iscell(meshfile)
    meshfile=meshfile{iscale};
end
pscale=2^(iscale-1);
rint=false;
if isfield(param,'reduced_integration')
    rint=param.reduced_integration;
end
[xi,yi,zi,conn,elt,selected]=ReadVTK(meshfile);
if isfield(param,'gluing_parameters')
    if isempty(param.gluing_parameters)
        gluings=GetGluingParameters(nmod,xi,yi,zi,conn,elt);
        param.gluing_parameters=gluings;
        save(fullfile('TMP',sprintf('%d_params',nmod)),'param','-append');
        check=0;
    end
   [xi,yi,zi]=GlueMesh(param.gluing_parameters,xi,yi,zi);
 xo=(xi+1-roi(1))/psample+dec(1);
yo=(yi+1-roi(3))/psample+dec(2);
zo=0;
else
load(fullfile('TMP','sample0'),'sizeim');
sizeim0=sizeim;
xo=(-(yi-sizeim0(1))+1-roi(1))/psample-1+dec(1);
yo=(xi-roi(3))/psample+1+dec(2);
zo=0;
end
 load(fullfile('TMP','sample0'),'sizeim');
sizeim0=sizeim;
load(fullfile('TMP',sprintf('sample0_%d',1-1)),'sizeim');

if any(xo<.5)||any(xo>sizeim(1)+.5)||any(yo<.5)||any(yo>sizeim(2)+0.5)
   display('WARNING, NODES OUT OF THE ROI')
    display(sprintf('IMAGE SIZE: %d %d',sizeim0/psample));
    display(sprintf('ROI: %d %d %d %d',roi));
    display(sprintf('MESH AREA: %g %g %g %g',(min(xo)-1)*psample+roi(1),(max(xo)-1)*psample+roi(1),(min(yo)-1)*psample+roi(3),(max(yo)-1)*psample+roi(3)))
    %   display('THE MESH IS SHRINKED INSIDE THE ROI')
    display('NODES OUT OF THE ROI ARE REMOVED')
    load(fullfile('TMP','sample0'),'im0');
    
        seg3=reshape(conn(elt==3,[1,2,1,3,2,3])',2,3*sum(elt==3))';
    seg4=reshape(conn(elt==4,[1,2,2,3,3,4,4,1])',2,4*sum(elt==4))';
    
    seg=unique(sort([seg3;seg4],2),'rows');
    
    
    figure

    colormap(gray)
    imagesc(im0')
    hold on;
    rect1=rectangle;
    set(rect1,'position',[roi(1),roi(3),roi(2)-roi(1),roi(4)-roi(3)],...
        'EdgeColor','yellow',...
        'LineStyle','--',...
        'LineWidth',2)
    plot((xo(seg)'-1)*psample+roi(1),(yo(seg)'-1)*psample+roi(3),'b+','LineWidth',0.5);
    axis xy
    axis image
    %            axis off;

    title(meshfile,'FontSize',24);
    
    
    
    found=find((xo>=0)&(xo<=sizeim(1)+1)&(yo>=0)&(yo<=sizeim(2)+1));
    nconn=zeros(size(conn));
    nelt=zeros(size(elt));
    nselected=ones(length(found),1);
    ielt=0;
    for ie=1:length(elt)
        inods=conn(ie,1:elt(ie));
        xn=xo(inods);
        yn=yo(inods);
        founde=find((xn>=0)&(xn<=sizeim(1)+1)&(yn>=0)&(yn<=sizeim(2)+1));
        if length(founde)==elt(ie)
            for in=1:elt(ie)
                nconn(ielt+1,in)=find(found==inods(in));
            end
            nelt(ielt+1)=elt(ie);
            ielt=ielt+1;
            if any(~selected(inods))
                founds=find(~selected(inods));
                for in=1:length(founds)
                    nselected(find(found==inods(founds(in))))=0;
                end
            end
       elseif any(~selected(inods))
            for in=1:elt(ie)
                if any(founde==in)
                    nselected(find(found==inods(in)))=0;
                end
            end

        end
    end
    xo=xo(found);
    yo=yo(found);
    conn=nconn(1:ielt,:);
    elt=nelt(1:ielt);
    selected=nselected;
end
%     xo=max(xo,1);
% xo=min(xo,sizeim(1));
% yo=max(yo,1);
% yo=min(yo,sizeim(2));
    if check
        
            
        seg3=reshape(conn(elt==3,[1,2,1,3,2,3])',2,3*sum(elt==3))';
    seg4=reshape(conn(elt==4,[1,2,2,3,3,4,4,1])',2,4*sum(elt==4))';
    
    seg=unique(sort([seg3;seg4],2),'rows');

        
    figure

    if exist(fullfile('TMP','sample0.mat'),'file')
        load(fullfile('TMP','sample0'));
    if exist('im0','var')
        colormap(gray)
    imagesc(im0')
    end
    end
    hold on;
    rect1=rectangle;
    set(rect1,'position',[roi(1),roi(3),roi(2)-roi(1),roi(4)-roi(3)],...
        'EdgeColor','yellow',...
        'LineStyle','--',...
        'LineWidth',2)
        plot((xo(seg)'-1)*psample+roi(1),(yo(seg)'-1)*psample+roi(3),'b','LineWidth',0.5);
        found=find(~selected);
        plot((xo(found)-1)*psample+roi(1),(yo(found)-1)*psample+roi(3),'ro','LineWidth',0.5,'MarkerSize',10);

    axis xy
    axis image
    %            axis off;
    title(meshfile,'FontSize',24);


end

if iscale>1
xo=((xo-0.5)/pscale)+0.5;
yo=((yo-0.5)/pscale)+0.5;
end

Smesh=floor([max(xo)-min(xo),max(yo)-min(yo)]);
Vmesh=floor([max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)]);
Nnodes=[length(xo),1,1];
Nelems=[size(conn,1),1,1];

ng=0;
if isfield(param,'nb_gauss_points')
    ng=param.nb_gauss_points;
end
ns=ones(1,2);
if isfield(param,'nb_sub_cells')
    ns=param.nb_sub_cells;
end

save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
disp(sprintf('Loading mesh for model %d...%6.2fs',nmod,toc));
if isfield(param,'enrichment')
    save(fullfile('TMP',sprintf('%d_meshx_%d',nmod,iscale-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');

end




end

