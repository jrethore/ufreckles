function CreateVICBasisFunction(iscale,nmod)
%persistent ny nx
load(fullfile('TMP','params'));
param0=param;
tic;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
pscale=2^(iscale-1);
inde=[];
sbasis=param.basis;
%if isempty(nx)
%end
switch sbasis
    case 'vic-beam'
        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
        [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
        if iscale==1
            p=1;
        else
            p=param.degree;
        end
        mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
        type_nurbs=~strcmp(param.continuity,'c0');
        [phi]=CreateNURBSBasis0D(mesh_file,sizeim,p,pscale,type_nurbs);
        wdetJ=1;
    case {'vic-nurbs'}
        load(fullfile('TMP',sprintf('sample0_%d',1-1)),'ls1');
        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim','nband','on');
%        mesh_file=fullfile('TMP',sprintf('%d_vicmesh_%d',nmod,iscale-1));
        mesh_file=fullfile('TMP',sprintf('%d_vicmesh_%d',nmod,1-1))
        type_nurbs=~strcmp(param.continuity,'c0');
        p=param.degree;
        %         if iscale==1
        %             p=param.degree;
        %         else
        %             p=1;
        %         end
        if length(sizeim)==2
            load(fullfile('TMP',sprintf('sample0_%d',1-1)),'nx','ny');
            nx=-diag(sparse(nx(:)));
            ny=-diag(sparse(ny(:)));

            [phi,dphi,ddphi]=CreateNURBSBasis0D(mesh_file,p,ls1,type_nurbs,param.closed,nband,sizeim);
            ls1i=ls1(nband(on));
            [phii,dphii,ddphii]=CreateNURBSBasis0D(mesh_file,p,ls1i,type_nurbs,param.closed);
            % phi1=phi;
            %         [phi,dphi]=CreateNURBSBasis0D(mesh_file,p,ls1,type_nurbs,param.closed);

            wdetJ=1;
            [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
            phix=nx*phi;
            phiy=ny*phi;


            Xi=Xi(:);
            Yi=Yi(:);

            Nddl_tot=size(phi,2);
            save(fullfile('TMP',sprintf('%d_phi_%d',nmod,iscale-1)),'phi','phii','dphii','ddphii','dphi','ddphi','sizeim','ls1i','Nddl_tot');
            save(fullfile('TMP',sprintf('%d_phix_%d',nmod,iscale-1)),'phix','Xi','Yi','wdetJ','inde','sizeim','Nddl_tot');
            save(fullfile('TMP',sprintf('%d_phiy_%d',nmod,iscale-1)),'phiy','sizeim','Nddl_tot');
        else
            load(fullfile('TMP',sprintf('sample0_%d',1-1)),'ls2');
            ls1=ls1(on);ls2=ls2(on);
            load(fullfile('TMP',sprintf('sample0_%d',1-1)),'nx','ny','nz');
            nx=-diag(sparse(nx(on)));
            ny=-diag(sparse(ny(on)));
            nz=-diag(sparse(nz(on)));
            [phi,dphix,dphiy]=CreateNURBSBasis2D(mesh_file,p,ls1,ls2,type_nurbs,param.closed);
            phix=nx*phi;
            phiy=ny*phi;
            phiz=nz*phi;
            wdetJ=1;
            [Xi,Yi,Zi]=ind2sub(sizeim,nband(on));

            Nddl_tot=size(phi,2);
            save(fullfile('TMP',sprintf('%d_3d_phi_%d',nmod,iscale-1)),'phi','dphix','dphiy','sizeim','Nddl_tot');
            save(fullfile('TMP',sprintf('%d_3d_phix_%d',nmod,iscale-1)),'phix','Xi','Yi','Zi','wdetJ','inde','sizeim','Nddl_tot');
            save(fullfile('TMP',sprintf('%d_3d_phiy_%d',nmod,iscale-1)),'phiy','sizeim','Nddl_tot');
            save(fullfile('TMP',sprintf('%d_3d_phiz_%d',nmod,iscale-1)),'phiz','sizeim','Nddl_tot');


        end

    otherwise
        error ('INVALID FUNCTIONAL BASIS');

end

disp(sprintf('Creating basis function for model %d...%6.2f s',nmod,toc));

end

