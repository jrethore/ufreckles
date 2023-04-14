function [face_nodes,face_elts,crackn]=TreatmentofEnrichment3D(ic,nmod,crack,front,zone)

iscale=1;
pscale=1;
load(fullfile('TMP',sprintf('sample0_%d',0)),'sizeim');
load(fullfile('TMP','params'),'param');
roi=param.roi;
check=1;
       load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','Nnodes','Nelems','conn','elt','Smesh');
       xo=xo+roi(1)-zone(1);
       yo=yo+roi(3)-zone(3);
       zo=zo+roi(5)-zone(5);
 
    sizec=size(crack);
    inbox=~((xo<1)|(xo>sizec(1))|(yo<1)|(yo>sizec(2))|(zo<1)|(zo>sizec(3)));
 inside_nodes=find(inbox);
 ind_inside=zeros(Nnodes);
 ind_inside(inside_nodes)=1:length(inside_nodes);
 
 xn=xo(inside_nodes);yn=yo(inside_nodes);zn=zo(inside_nodes);
 crackn=mexInterpLinear3D(xn,yn,zn,crack);
 frontn=mexInterpLinear3D(xn,yn,zn,front);
    hnodes=2*double(crackn>=0)-1;
    face_elts=[];
    selected_nodes=zeros(Nnodes);
   for ie=1:prod(Nelems)
        inods=conn(ie,1:elt(ie));
        if ~any(~inbox(inods))
        xn=xo(inods);yn=yo(inods);zn=zo(inods);
        iin=ind_inside(inods);
            above_below=abs(mean(hnodes(iin)));
            if (above_below<1)
                infront_behind=frontn(iin)<0;
                if any(infront_behind)
                     face_elts=[face_elts,ie];
                     selected_nodes(inods(infront_behind))=1;
                end
            end
        end
    end
    face_nodes=find(selected_nodes);
 xn=xo(face_nodes);yn=yo(face_nodes);zn=zo(face_nodes);
  crackn=mexInterpLinear3D(xn,yn,zn,crack);
 frontn=mexInterpLinear3D(xn,yn,zn,front);
mesh_file=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1));
enriched_face_nodes=AddOneNodeLayer(conn,face_nodes);
phi=CreateFiniteElementBasis3D(mesh_file,sizeim,pscale,enriched_face_nodes,'sub_cells',false,face_elts);
xp=phi*xo(enriched_face_nodes);yp=phi*yo(enriched_face_nodes);zp=phi*zo(enriched_face_nodes);
 crackp=mexInterpLinear3D(xp,yp,zp,crack);
 frontp=mexInterpLinear3D(xp,yp,zp,front);
if check

           figure
           plot3(xo,yo,zo,'k+','MarkerSize',7)
           hold on
%           plot3(xo(inside_nodes),yo(inside_nodes),zo(inside_nodes),'rx','MarkerSize',7)
           plot3(xo(face_nodes),yo(face_nodes),zo(face_nodes),'bo','MarkerSize',7)
          xlabel('X');ylabel('Y');zlabel('Z');
          axis equal
          print('-djpeg',[param.result_file,'-enrichment.jpg']);
          print('-depsc2',[param.result_file,'-enrichment.eps']);
           pause(0.1)
end
        enriched_nodes=face_nodes;
    save(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'enriched_nodes','face_elts','face_nodes','frontn','crackn','zone','crackp','frontp');
end