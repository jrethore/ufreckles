function CreateGradBasisFunction(iscale,nmod)
load(fullfile('TMP','params'));
param0=param;
inde=[];
    tic;
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
    pscale=2^(iscale-1);
    sbasis=param.basis;
if strcmp(sbasis,'nurbs')&&iscale>1
    sbasis='fem';
end   
    switch sbasis
        case 'KM'
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            wdetJ=1;
            ind=param.km_indices;
            kappa=param.kolosov;
            modes=param.modes;
            if isfield(param,'crack_id')
                ic=param.crack_id;
            else
                ic=1;
            end
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
            sizemask=numel(mask);
            if sizemask==1
                ipix=[];
            else
                ipix=find(diag(mask)>0);
            end
            [dkmdx,dkmdy]=CreateGradKMBasis(modes,ind,kappa,ipix,ic);
            epsxx=real(dkmdx);
            save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,(iscale-1))),'epsxx','wdetJ','sizeim');
            epsyy=imag(dkmdy);
            save(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,(iscale-1))),'epsyy','sizeim');
            epsxy=imag(dkmdx)+real(dkmdy);
            save(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,(iscale-1))),'epsxy','sizeim');


            %             [km]=CreateKMBasis(modes,ind,kappa,[]);
            %
            %
            %
            %
            %             [kmy,kmx]=gradient(reshape(real(km(:,1)),sizeim));
            %             kmx=kmx.*reshape(diag(mask),sizeim);
            %             kmy=kmy.*reshape(diag(mask),sizeim);
            %             figure
            %             toto=reshape(real(dkmdx(:,1)),sizeim);
            %             imagesc(toto(2:sizeim(1)-1,2:sizeim(2)-1))
            %             colorbar;
            %             figure
            %             imagesc(kmx(2:sizeim(1)-1,2:sizeim(2)-1))
            %             colorbar;
            %             figure
            %             imagesc(abs(kmx(2:sizeim(1)-1,2:sizeim(2)-1)-toto(2:sizeim(1)-1,2:sizeim(2)-1)))
            %             colorbar;


        case 'KM+'
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            wdetJ=1;
            ind=[0,2];
            kappa=param.kolosov;
            modes=param.modes;
            lmin=param.cz_length;
            if isfield(param,'crack_id')
                ic=param.crack_id;
            else
                ic=1;
            end
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
            sizemask=numel(mask);
            if sizemask==1
                ipix=[];
            else
                ipix=find(diag(mask)>0);
            end
            load(fullfile('TMP',sprintf('%d_iphi_%d',nmod,iscale-1)),'iphi','ztip','dz');

            epsxx=[];epsxy=[];epsyy=[];
            if any(modes==1)
                for iz=1:length(ztip)
                    [dkmdx,dkmdy]=CreateGradKMBasis(1,1,kappa,ipix,ic,ztip(iz));

                    epsxx=[epsxx,real(dkmdx)*dz(iz)];
                    epsyy=[epsyy,imag(dkmdy)*dz(iz)];
                    epsxy=[epsxy,imag(dkmdx)*dz(iz)+real(dkmdy)*dz(iz)];

                end
            end
            if any(modes==2)
                for iz=1:length(ztip)
                    [dkmdx,dkmdy]=CreateGradKMBasis(2,1,kappa,ipix,ic,ztip(iz));

                    epsxx=[epsxx,real(dkmdx)*dz(iz)];
                    epsyy=[epsyy,imag(dkmdy)*dz(iz)];
                    epsxy=[epsxy,imag(dkmdx)*dz(iz)+real(dkmdy)*dz(iz)];

                end
            end
            if length(modes)==2
                iphi=blkdiag(iphi,iphi);
            end
            epsxx=epsxx*iphi;
            epsyy=epsyy*iphi;
            epsxy=epsxy*iphi;

            [dkmdx,dkmdy]=CreateGradKMBasis([1,2],ind,kappa,ipix,ic);
            epsxx=[real(dkmdx),epsxx];
            epsyy=[imag(dkmdy),epsyy];
            epsxy=[imag(dkmdx)+real(dkmdy),epsxy];
%             if any(modes==2)
%                     [dkmdx,dkmdy]=CreateGradKMBasis(2,1,kappa,ipix,ic);
% 
%                     epsxx=[epsxx,real(dkmdx)];
%                     epsyy=[epsyy,imag(dkmdy)];
%                     epsxy=[epsxy,imag(dkmdx)+real(dkmdy)];
% 
%             end

            if isfield(param,'sub_indices')
                sind=param.sub_indices;
                if ~isempty(sind)
                    [dkmdx,dkmdy]=CreateGradKMBasis(modes,sind,kappa,ipix,ic,lmin);
                    epsxx=[epsxx,real(dkmdx)];
                    epsyy=[epsyy,imag(dkmdy)];
                    epsxy=[epsxy,imag(dkmdx)+real(dkmdy)];
                end
            end

            save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,(iscale-1))),'epsxx','wdetJ','sizeim');
            save(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,(iscale-1))),'epsyy','sizeim');
            save(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,(iscale-1))),'epsxy','sizeim');



%         case 'sfem'
% 
%             load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
%             mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
%             load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
%             [dphidx,dphidy]=CreateGradSmoothFiniteElementBasis(mesh_file,sizeim,pscale,unmasked_nodes,'Gauss_points');
%             phi0=sparse(size(dphidx,1),size(dphidx,2));
%             epsxx=[dphidx,phi0];
%             save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','sizeim');
%             epsyy=[phi0,dphidy];
%             save(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy','sizeim');
%             epsxy=[dphidy,dphidx];
%             save(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'epsxy','sizeim');
        case 'fem'

            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
            load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'cut_nodes','unmasked_nodes');
            [dphidx,dphidy]=CreateGradFiniteElementBasis(mesh_file,sizeim,pscale,unmasked_nodes,'Gauss_points');
            [wdetJ,inde]=GetWeigthDetJ(mesh_file,sizeim,1,'Gauss_points');
            if (~isempty(cut_nodes))
             [dphidxp,dphidyp]=CreateGradFiniteElementBasis(mesh_file,sizeim,pscale,cut_nodes,'pixels',1);   
            wdetJp=GetWeigthDetJ(mesh_file,sizeim,1,'pixels');
            wdetJ=[wdetJ*ones(size(dphidx,2),1);wdetJp*ones(size(dphidxp,2),1)];
            wdetJ=diag(sparse(wdetJ));
            dphidx=[dphidx;dphidxp(:,unmasked_nodes)];
            dphidy=[dphidy;dphidyp(:,unmasked_nodes)];
            end
            phi0=sparse(size(dphidx,1),size(dphidx,2));
            epsxx=[dphidx,phi0];
            Uxx=[dphidx,phi0];
            Uxy=[dphidy,phi0];
            Uyx=[phi0,dphidx];
            Uyy=[phi0,dphidy];
            save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'Uxx','epsxx','wdetJ','sizeim');
            epsyy=[phi0,dphidy];
            save(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'Uyy','epsyy','sizeim');
            epsxy=[dphidy,dphidx];
            save(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uyx','Uxy','epsxy','sizeim');
        case 'nurbs'
            
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            p=param.degree;
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
            [dphidx,dphidy,Xi,Yi,wdetJ]=CreateGradNURBSBasis(mesh_file,p,'Gauss_points',1,'physical');
            phi0=sparse(size(dphidx,1),size(dphidx,2));
            epsxx=[dphidx,phi0];
            Uxx=[dphidx,phi0];
            Uxy=[dphidy,phi0];
            Uyx=[phi0,dphidx];
            Uyy=[phi0,dphidy];
            save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'Uxx','epsxx','wdetJ','Xi','Yi','sizeim');
            epsyy=[phi0,dphidy];
            save(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'Uyy','epsyy','sizeim');
            epsxy=[dphidy,dphidx];
            save(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uyx','Uxy','epsxy','sizeim');
    [dphidxx,dphidyy,dphidxy,xg,yg,wdetJ]=CreateGradGradNURBSBasis(mesh_file,p,'Gauss_points',1,'physical');
    phi0=sparse(size(dphidxx,1),size(dphidxx,2));
    epsxxx=[dphidxx,phi0];
    epsxxy=[dphidxy,phi0];
    epsyyy=[phi0,dphidyy];
    epsyyx=[phi0,dphidxy];
    epsxyx=[dphidxy,dphidxx];
    epsxyy=[dphidyy,dphidxy];
            save(fullfile('TMP',sprintf('%d_epsijk_%d',nmod,10*(iscale-1))),'epsxxx','epsxxy','epsxyx','epsxyy','epsyyx','epsyyy');

        case 'btri'
            
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
            [phi,Xi,Yi,wdetJ]=CreateBezierTriangleBasis(mesh_file,'Gauss_points');
            [dphidx,dphidy,wdetJ]=CreateGradBezierTriangleBasis(mesh_file,'Gauss_points');
            phi0=sparse(size(dphidx,1),size(dphidx,2));
            epsxx=[dphidx,phi0];
            Uxx=[dphidx,phi0];
            Uxy=[dphidy,phi0];
            Uyx=[phi0,dphidx];
            Uyy=[phi0,dphidy];
            save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'Uxx','epsxx','wdetJ','sizeim');
            epsyy=[phi0,dphidy];
            save(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'Uyy','epsyy','sizeim');
            epsxy=[dphidy,dphidx];
            save(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uyx','Uxy','epsxy','sizeim');


        otherwise
            error ('INVALID FUNCTIONAL BASIS');

    end

    
    if (iscale==1)&&(isfield(param,'enrichment'))
        if isfield(param,'mesh_size_ratio')
            assert(param.mesh_size_ratio==1);
        end
        enr=param.enrichment;
        mesh_file=fullfile('TMP',sprintf('%d_meshx_%d',nmod,iscale-1));
        load(mesh_file,'ng','Nelems');
        wdetJ=GetWeigthDetJ(mesh_file,sizeim,1,'sub_cells',1);
        npt=size(wdetJ,1);
        if npt==1
            npt=prod(sizeim);
        end
        epsxx=sparse(npt,size(epsxx,2));
        epsyy=sparse(npt,size(epsxx,2));
        epsxy=sparse(npt,size(epsxx,2));
        Uxy=sparse(npt,size(epsxx,2));
        Uyx=sparse(npt,size(epsxx,2));
        epscxx=sparse(npt,size(epsxx,2));
        epscyy=sparse(npt,size(epsxx,2));
        epscxy=sparse(npt,size(epsxx,2));
        Ucxy=sparse(npt,size(epsxx,2));
        Ucyx=sparse(npt,size(epsxx,2));

        if ~iscell(param0.levelset_file)
            nbfis=1;
        else
            nbfis=numel(param0.levelset_file);
        end
        wdetJ=0*wdetJ;
        for ic=1:nbfis

            load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'enriched_nodes','enriched_elts','tip_nodes','face_nodes','enriched_pixels');
            wdetJi=GetWeigthDetJ(mesh_file,sizeim,1,'sub_cells',enriched_elts);
            wdetJ=wdetJ+wdetJi;
            phi=CreateFiniteElementBasis(mesh_file,sizeim,pscale,enriched_nodes,'sub_cells',true);
            load(mesh_file,'xo','yo');
            xp=phi*xo;yp=phi*yo;
            %             figure
            %             imagesc(reshape(sum(phi,2),sizeim))
            % figure
            % plot(sum(phi,2))
            % % any(sum(phi,2)>1)
            %              figure
            %              scatter(xp,yp,1+sum(phi,2))
% 
%             figure
%             plot(xo,yo,'b+')
%             hold on
%             plot(xp,yp,'ko')
%             axis xy
%             axis image
            %         pause
            [dphidx,dphidy]=CreateGradFiniteElementBasis(mesh_file,sizeim,pscale,enriched_nodes,'sub_cells',true);
            if isempty(unmasked_nodes)
                phi0=sparse(size(phi,1),size(phi,2));
                epscxx=epscxx+[dphidx,phi0];
                epscyy=epscyy+[phi0,dphidy];
                epscxy=epscxy+[dphidy,dphidx];
                Ucxy=Ucxy+[dphidy,phi0];
                Ucyx=Ucyx+[phi0,dphidx];

            else
                phi0=sparse(size(phi,1),length(unmasked_nodes));
                epscxx=epscxx+[dphidx(:,unmasked_nodes),phi0];
                epscyy=epscyy+[phi0,dphidy(:,unmasked_nodes)];
                epscxy=epscxy+[dphidy(:,unmasked_nodes),dphidx(:,unmasked_nodes)];
                Ucxy=Ucxy+[dphidy(:,unmasked_nodes),phi0];
                Ucyx=Ucyx+[phi0,dphidx(:,unmasked_nodes)];
            end
            do_disc_enrich=true;
            if ~(strcmp(enr,'crack_disc'))
                ind=param.km_indices;
                kappa=param.kolosov;
                modes=param.modes;
                [dkmdx,dkmdy]=CreateGradKMBasis(modes,ind,kappa,enriched_pixels,ic);
                [km]=CreateKMBasis(modes,ind,kappa,enriched_pixels,ic);

                strat=param.enrichment_type;

                if isempty(face_nodes)
                    do_disc_enrich=false;
                end

                if isempty(enriched_pixels)&&strcmp(strat,'zone')
                    phikm=1;
                    dphikmdx=0;
                    dphikmdy=0;
                else
                    phikm=phi(:,tip_nodes);
                    dphikmdx=dphidx(:,tip_nodes);
                    dphikmdy=dphidy(:,tip_nodes);
                    if strcmp(strat,'zone')
                        phikm=sum(phikm,2);
                        dphikmdx=sum(dphikmdx,2);
                        dphikmdy=sum(dphikmdy,2);
                    end
                    %                 if strcmp(strat,'zone')
                    %                        load(fullfile('TMP',sprintf('%d_cutoff_%d',nmod,iscale-1),'cutoff','dcdx','dcdy');
                    %                     phikm=cutoff;
                    %                     dphikmdx=dcdx;
                    %                     dphikmdy=dcdy;
                    %                     clear cutoff dcdx dcdy
                    %                 else
                    %                 phikm=phi(:,tip_nodes);
                    %                 dphikmdx=dphidx(:,tip_nodes);
                    %                 dphikmdy=dphidy(:,tip_nodes);
                    %                 end
                end
                if strcmp(strat,'zone')
                    epsxx=[epsxx,diag(sparse(phikm))*real(dkmdx)+diag(sparse(dphikmdx))*real(km)];
                    epsyy=[epsyy,diag(sparse(phikm))*imag(dkmdy)+diag(sparse(dphikmdy))*imag(km)];
                    epsxy=[epsxy,diag(sparse(phikm))*(real(dkmdy)+imag(dkmdx))+diag(sparse(dphikmdy))*real(km)+diag(sparse(dphikmdx))*imag(km)];
                    Uxy=[Uxy,diag(sparse(phikm))*(real(dkmdy))+diag(sparse(dphikmdy))*real(km)];
                    Uyx=[Uyx,diag(sparse(phikm))*(imag(dkmdx))+diag(sparse(dphikmdx))*imag(km)];
                elseif strcmp(strat,'node')
                    for kk=1:size(km,2)
                        epsxx= [epsxx,diag(sparse(real(dkmdx(:,kk))))*phikm+diag(sparse(real(km(:,kk))))*dphikmdx];
                        epsyy= [epsyy,diag(sparse(imag(dkmdy(:,kk))))*phikm+diag(sparse(imag(km(:,kk))))*dphikmdy];
                        epsxy= [epsxy,diag(sparse(real(dkmdy(:,kk))+imag(dkmdx(:,kk))))*phikm+diag(sparse(imag(km(:,kk))))*dphikmdx+diag(sparse(real(km(:,kk))))*dphikmdy];
                        Uxy=[Uxy,diag(sparse(phikm))*(real(dkmdy(:,kk)))+diag(sparse(dphikmdy))*real(km(:,kk))];
                        Uyx=[Uyx,diag(sparse(phikm))*(imag(dkmdx(:,kk)))+diag(sparse(dphikmdx))*imag(km(:,kk))];
                    end
                end
                clear km dkmdx dkmdy dphikmdx dphikmdy

            end

            if do_disc_enrich
                load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack');
                if isfield(param,'shift_enrichment')
                    opti_shift=param.shift_enrichment;
                else
                    opti_shift=0;
                end
                if ~opti_shift

                    crack=interp2(crack,yp,xp,'*linear');
                    heaviside=2*double(crack>=0)-1;
                    %                     crack=crack(enriched_pixels);
                    %                     heaviside=zeros(prod(sizeim),1);
                    %                     heaviside(enriched_pixels)=2*double(crack>=0)-1;
                    %                     figure
                    %                     imagesc(full(reshape(heaviside,sizeim)))

%                     figure
%                     scatter(xp,yp,10+0*yp,heaviside)
%                     hold on
%                     plot(xo,yo,'+','MarkerSize',10,'LineWidth',2)
                    % figure
                    % scatter(xp,yp,10+0*yp,sum(abs(dphidx),2))
                    % face_nodes
                    % figure
                    % scatter(xp,yp,10+0*yp,sum(abs(dphidx(:,face_nodes)),2))


                    clear crack
                    heaviside=diag(sparse(heaviside));
                    dphihdx=heaviside*dphidx(:,face_nodes);
                    dphihdy=heaviside*dphidy(:,face_nodes);

                else

                    disp('WARNING SHIFTED ENRICHMENT !!!');
                    load(mesh_file,'xo','yo');
                    hn=interp2(crack,yo,xo,'*linear');
                    hn=double(hn>=0);
                    hn=diag(sparse(hn(face_nodes)));
                    %                     crack=crack(enriched_pixels);
                    %                     heaviside=zeros(prod(sizeim),1);
                    %                     heaviside(enriched_pixels)=double(crack>=0);
                    crack=interp2(crack,yp,xp,'*linear');
                    heaviside=double(crack>=0);
                    %                     figure
                    %                     imagesc(reshape(heaviside,sizeim)')
                    %                     hold on
                    %                     plot(xo,yo,'ko')
                    %             axis xy
                    %         axis image
                    %         axis off;
                    clear crack
                    heaviside=diag(sparse(heaviside));
                    dphihdx=heaviside*dphidx(:,face_nodes)-dphidx(:,face_nodes)*hn;
                    dphihdy=heaviside*dphidy(:,face_nodes)-dphidy(:,face_nodes)*hn;
                end
                phih0=sparse(size(dphihdx,1),size(dphihdx,2));
                epsxx=[epsxx,dphihdx,phih0];
                epsyy=[epsyy,phih0,dphihdy];
                epsxy=[epsxy,dphihdy,dphihdx];
                Uxy=[Uxy,dphihdy,phih0];
                Uyx=[Uyx,phih0,dphihdx];
                clear dphidx dphidy phih0 heaviside
            end
        end
        eps0=sparse(npt,size(epsxx,2)-size(epscxx,2));
        epscxx=[epscxx,eps0];
        epscyy=[epscyy,eps0];
        epscxy=[epscxy,eps0];
        Ucxy=[Ucxy,eps0];
        Ucyx=[Ucyx,eps0];

        epsxx=epscxx+epsxx;
        epsyy=epscyy+epsyy;
        epsxy=epscxy+epsxy;
        Uxy=Ucxy+Uxy;
        Uyx=Ucyx+Uyx;
        save(fullfile('TMP',sprintf('%d_xepsxx_%d',nmod,(iscale-1))),'epsxx','wdetJ','inde','sizeim');
        save(fullfile('TMP',sprintf('%d_xepsyy_%d',nmod,(iscale-1))),'epsyy','sizeim');
        save(fullfile('TMP',sprintf('%d_xepsxy_%d',nmod,(iscale-1))),'Uxy','Uyx','epsxy','sizeim');


    end
    Nddl_tot=size(epsxx,2);
    %    disp(sprintf('Creating gradient basis function for model %d...%6.2f s',nmod,toc));
save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'Nddl_tot','-append');


end