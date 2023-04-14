function CreateGradBasisFunction25D(iscale,nmod)
load(fullfile('TMP','params'),'param');
param0=param;
roi=param0.roi;
tic;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
pscale=2^(iscale-1);
sbasis=param.basis;
load(fullfile('TMP','sample0'),'sizeim');

switch sbasis

    case 'fem'
        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
        mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
        [dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis25D(mesh_file,sizeim,pscale,[],'Gauss_points');
        wdetJ=GetWeigthDetJ25D(mesh_file,sizeim,1,'Gauss_points');
        Nddl_tot=3*size(dphidx,2);
        phi0=sparse(size(dphidx,1),size(dphidx,2));
        epsxx=[dphidx,phi0,phi0];
        save(fullfile('TMP',sprintf('%d_25d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ','sizeim');
        epsyy=[phi0,dphidy,phi0];
        save(fullfile('TMP',sprintf('%d_25d_epsyy_%d',nmod,10*(iscale-1))),'epsyy','sizeim');
        epsxy=[dphidy,dphidx,phi0];
            Uxy=[dphidy,phi0,phi0];
            Uyx=[phi0,dphidx,phi0];
        save(fullfile('TMP',sprintf('%d_25d_epsxy_%d',nmod,10*(iscale-1))),'Uxy','Uyx','epsxy','sizeim');
        epsxz=[dphidz,phi0,dphidx];
        Uxz=[dphidz,phi0,phi0];
        Uzx=[phi0,phi0,dphidx];
        save(fullfile('TMP',sprintf('%d_25d_epsxz_%d',nmod,10*(iscale-1))),'Uxz','Uzx','epsxz','sizeim');
        epsyz=[phi0,dphidz,dphidy];
        Uyz=[phi0,dphidz,phi0];
        Uzy=[phi0,phi0,dphidy];
        save(fullfile('TMP',sprintf('%d_25d_epsyz_%d',nmod,10*(iscale-1))),'Uyz','Uzy','epsyz','sizeim');
        epszz=[phi0,phi0,dphidz];
        save(fullfile('TMP',sprintf('%d_25d_epszz_%d',nmod,10*(iscale-1))),'epszz','sizeim');
    case 'nurbs'

        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            p=param.degree;
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
            [dphidx,dphidy,dphidz,Xi,Yi,Zi,wdetJ]=CreateGradNURBSBasis25D(mesh_file,p,'Gauss_points',1,'physical');
        Nddl_tot=3*size(dphidx,2);
        phi0=sparse(size(dphidx,1),size(dphidx,2));
        epsxx=[dphidx,phi0,phi0];
        save(fullfile('TMP',sprintf('%d_25d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ','sizeim');
        epsyy=[phi0,dphidy,phi0];
        save(fullfile('TMP',sprintf('%d_25d_epsyy_%d',nmod,10*(iscale-1))),'epsyy','sizeim');
        epsxy=[dphidy,dphidx,phi0];
            Uxy=[dphidy,phi0,phi0];
            Uyx=[phi0,dphidx,phi0];
        save(fullfile('TMP',sprintf('%d_25d_epsxy_%d',nmod,10*(iscale-1))),'Uxy','Uyx','epsxy','sizeim');
        epsxz=[dphidz,phi0,dphidx];
        Uxz=[dphidz,phi0,phi0];
        Uzx=[phi0,phi0,dphidx];
        save(fullfile('TMP',sprintf('%d_25d_epsxz_%d',nmod,10*(iscale-1))),'Uxz','Uzx','epsxz','sizeim');
        epsyz=[phi0,dphidz,dphidy];
        Uyz=[phi0,dphidz,phi0];
        Uzy=[phi0,phi0,dphidy];
        save(fullfile('TMP',sprintf('%d_25d_epsyz_%d',nmod,10*(iscale-1))),'Uyz','Uzy','epsyz','sizeim');
        epszz=[phi0,phi0,dphidz];
        save(fullfile('TMP',sprintf('%d_25d_epszz_%d',nmod,10*(iscale-1))),'epszz','sizeim');

    otherwise
        error ('INVALID FUNCTIONAL BASIS');

end
disp(sprintf('Creating basis function for model %d...%6.2f s',nmod,toc));
end
