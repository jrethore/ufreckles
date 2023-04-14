function CreateGradBasisFunction3D(iscale,nmod)
%sbasis,frestart
load(fullfile('TMP','params'));
param0=param;
tic;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
pscale=2^(iscale-1);
sbasis=param.basis;
sizeim=ones(1,3);
switch sbasis
    case 'fem'

        mesh_file=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1));
        [dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis3D(mesh_file,sizeim,pscale,[],'Gauss_points');
        wdetJ=GetWeigthDetJ(mesh_file,sizeim,1,'Gauss_points');

    case 'nurbs'
        mesh_file=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1));
        load(mesh_file,'degree');
            [dphidx,dphidy,dphidz,Xi,Yi,Zi,wdetJ]=CreateGradNURBSBasis3D(mesh_file,degree,'Gauss_points',1,'physical');

    otherwise
        error ('INVALID FUNCTIONAL BASIS');

end
        Nddl_tot=3*size(dphidx,2);
        phi0=sparse(size(dphidx,1),size(dphidx,2));
        epsxx=[dphidx,phi0,phi0];
        save(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,1*(iscale-1))),'epsxx','wdetJ','sizeim');
        epsyy=[phi0,dphidy,phi0];
        save(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,1*(iscale-1))),'epsyy','sizeim');
        epsxy=[dphidy,dphidx,phi0];
            Uxy=[dphidy,phi0,phi0];
            Uyx=[phi0,dphidx,phi0];
        save(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,1*(iscale-1))),'Uxy','Uyx','epsxy','sizeim');
        epsxz=[dphidz,phi0,dphidx];
        Uxz=[dphidz,phi0,phi0];
        Uzx=[phi0,phi0,dphidx];
        save(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,1*(iscale-1))),'Uxz','Uzx','epsxz','sizeim');
        epsyz=[phi0,dphidz,dphidy];
        Uyz=[phi0,dphidz,phi0];
        Uzy=[phi0,phi0,dphidy];
        save(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,1*(iscale-1))),'Uyz','Uzy','epsyz','sizeim');
        epszz=[phi0,phi0,dphidz];
        save(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,1*(iscale-1))),'epszz','sizeim');

        if (iscale==1)&&(isfield(param,'enrichment'))

    mesh_file=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1));
    load(mesh_file,'ng','Nelems');
    wdetJo=GetWeigthDetJ(mesh_file,sizeim,1,'sub_cells',1);
    npt=size(wdetJo,1);
    epsxx=sparse(npt,size(epsxx,2));
    epsyy=sparse(npt,size(epsxx,2));
    epsxy=sparse(npt,size(epsxx,2));
    epsxz=sparse(npt,size(epsxx,2));
    epsyz=sparse(npt,size(epsxx,2));
    epszz=sparse(npt,size(epsxx,2));
    epscxx=sparse(npt,size(epsxx,2));
    epscyy=sparse(npt,size(epsxx,2));
    epscxy=sparse(npt,size(epsxx,2));
    epscxz=sparse(npt,size(epsxx,2));
    epscyz=sparse(npt,size(epsxx,2));
    epsczz=sparse(npt,size(epsxx,2));
        Uxy=sparse(npt,size(epsxx,2));
        Uyx=sparse(npt,size(epsxx,2));
        Uxz=sparse(npt,size(epsxx,2));
        Uzx=sparse(npt,size(epsxx,2));
        Uzy=sparse(npt,size(epsxx,2));
        Uyz=sparse(npt,size(epsxx,2));
        Ucxy=sparse(npt,size(epsxx,2));
        Ucyx=sparse(npt,size(epsxx,2));
        Ucxz=sparse(npt,size(epsxx,2));
        Uczx=sparse(npt,size(epsxx,2));
        Uczy=sparse(npt,size(epsxx,2));
        Ucyz=sparse(npt,size(epsxx,2));
    wdetJ=sparse(npt,1);
    if ~iscell(param.levelset_file)
        nbfis=1;
    else
        nbfis=numel(param.levelset_file);
    end
    roi=param0.roi;
    for ic=1:nbfis
        load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'zone','face_nodes','face_elts','crackn','crackp');
        wdetJo=GetWeigthDetJ(mesh_file,sizeim,1,'sub_cells',face_elts);

        [dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis3D(mesh_file,sizeim,pscale,face_nodes,'sub_cells',true,face_elts);
        phi0=sparse(size(dphidx,1),size(dphidx,2));
        epscxx=epscxx+[dphidx,phi0,phi0];
        epscyy=epscyy+[phi0,dphidy,phi0];
        epscxy=epscxy+[dphidy,dphidx,phi0];
        epscxz=epscxz+[dphidz,phi0,dphidx];
        epscyz=epscyz+[phi0,dphidz,dphidy];
        epsczz=epsczz+[phi0,phi0,dphidz];
        Ucxy=Ucxy+[dphidy,phi0,phi0];
        Ucyx=Ucyx+[phi0,dphidx,phi0];
        Ucxz=Ucxz+[dphidz,phi0,phi0];
        Uczx=Uczx+[phi0,phi0,dphidx];
        Uczy=Uczy+[phi0,phi0,dphidy];
        Ucyz=Ucyz+[phi0,dphidz,phi0];
        hn=double(crackn>=0);
        hn=diag(sparse(hn));
        heaviside=double(crackp>=0);
        heaviside=diag(sparse(heaviside));
        dphihdx=heaviside*dphidx(:,face_nodes)-dphidx(:,face_nodes)*hn;
        dphihdy=heaviside*dphidy(:,face_nodes)-dphidy(:,face_nodes)*hn;
        dphihdz=heaviside*dphidz(:,face_nodes)-dphidz(:,face_nodes)*hn;
        phih0=sparse(size(dphihdx,1),size(dphihdx,2));
        epsxx=[epsxx,dphihdx,phih0,phih0];
        epsyy=[epsyy,phih0,dphihdy,phih0];
        epsxy=[epsxy,dphihdy,dphihdx,phih0];
        epsxz=[epsxz,dphihdz,phih0,dphihdx];
        epsyz=[epsyz,phih0,dphihdz,dphihdy];
        epszz=[epszz,phih0,phih0,dphihdz];
        Uxy=[Uxy,dphihdy,phih0,phih0];
        Uyx=[Uyx,phih0,dphihdx,phih0];
        Uxz=[Uxz,dphihdz,phih0,phih0];
        Uzx=[Uzx,phih0,phih0,dphihdx];
        Uzy=[Uzy,phih0,phih0,dphihdy];
        Uyz=[Uyz,phih0,dphihdz,phih0];
        wdetJ=wdetJ+diag(wdetJo);
    end
        eps0=sparse(npt,size(epsxx,2)-size(epscxx,2));        
                   epscxx=[epscxx,eps0];
                    epscyy=[epscyy,eps0];
                   epscxy=[epscxy,eps0];
                   epscxz=[epscxz,eps0];
                   epscyz=[epscyz,eps0];
                   epsczz=[epsczz,eps0];
                   Ucxy=[Ucxy,eps0];
                   Ucxz=[Ucxz,eps0];
                   Uczy=[Uczy,eps0];
                   Ucyx=[Ucyx,eps0];
                   Uczx=[Uczx,eps0];
                   Ucyz=[Ucyz,eps0];
    wdetJ=diag(wdetJ);
    epsxx=epscxx+epsxx;
    epsyy=epscyy+epsyy;
    epsxy=epscxy+epsxy;
    epsxz=epscxz+epsxz;
    epsyz=epscyz+epsyz;
    epszz=epsczz+epszz;
    Uxy=Ucxy+Uxy;
    Uyx=Ucyx+Uyx;
    Uxz=Ucxz+Uxz;
    Uzx=Uczx+Uzx;
    Uzy=Uczy+Uzy;
    Uyz=Ucyz+Uyz;
    Nddl_tot=size(epsxx,2);
    save(fullfile('TMP',sprintf('%d_3d_xepsxx_%d',nmod,1*(iscale-1))),'epsxx','wdetJ','sizeim','Nddl_tot');
    save(fullfile('TMP',sprintf('%d_3d_xepsyy_%d',nmod,1*(iscale-1))),'epsyy','sizeim');
    save(fullfile('TMP',sprintf('%d_3d_xepsxy_%d',nmod,1*(iscale-1))),'epsxy','sizeim','Uxy','Uyx');
    save(fullfile('TMP',sprintf('%d_3d_xepsxz_%d',nmod,1*(iscale-1))),'epsxz','sizeim','Uxz','Uzx');
    save(fullfile('TMP',sprintf('%d_3d_xepsyz_%d',nmod,1*(iscale-1))),'epsyz','sizeim','Uzy','Uyz');
    save(fullfile('TMP',sprintf('%d_3d_xepszz_%d',nmod,1*(iscale-1))),'epszz','sizeim');
end
save(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,1*(iscale-1))),'Nddl_tot','-append');
disp(sprintf('Creating gradient basis function 3D for model %d...%6.2f s',nmod,toc));


end