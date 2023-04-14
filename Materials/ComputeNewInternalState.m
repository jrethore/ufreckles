function ComputeNewInternalState(nmod,U,keep)
if nargin <3, keep=1;end
persistent nmodo epsxx epsyy epszz epsxy dphincz dphitcz
persistent Uxz Uzx Uxy Uyx Uyz Uzy notonboundary
load(fullfile('TMP','params'),'param');
param0=param;
nb_sub_cycles=1;
if isfield(param0,'nb_sub_cycles')
    nb_sub_cycles=param0.nb_sub_cycles;
end
helper='none';
if isfield(param0,'helper')
    helper=param0.helper;
end
if strcmp(helper,'none')
    pix2m=param.pixel_size;
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
    dflag=isfield(param,'extrusion_parameters')||(length(param0.roi)==6);
    nlflag=false;
    if isfield(param,'nlgeom')
        nlflag=param.nlgeom;
    end
    iscale=1;
    matmod=param.material_model;
    restart=1;
    if ~isempty(nmodo)
        if (nmodo==nmod)
            restart=0;
        else
            restart=1;
        end
    end
    switch matmod
        case 'elastic_homogeneous_isotropic_damage'
            load(fullfile('TMP',sprintf('%d_matmod',nmod)),'model','eeqmaxo');
            if restart
                load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,iscale-1)),'epsxx');
                load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,iscale-1)),'epsyy');
                load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,iscale-1)),'epsxy');
                nmodo=nmod;
            end
            exxm=epsxx*U;
            eyym=epsyy*U;
            exym=0.5*epsxy*U;
            %        if any(isnan(U))
            %            found=find(isnan(exxm));
            %            exxm(found)=0;
            %            exym(found)=0;
            %            eyym(found)=0;
            %        end
            if isfield(model,'nuo')
                nuo=model.nuo;
            else
                nuo=0;
            end
            switch model.equivalent_strain
                case 'Mazars'
                    e1=0.5*(exxm+eyym)+sqrt((0.5*(exxm-eyym)).^2+exym.^2);
                    e2=0.5*(exxm+eyym)-sqrt((0.5*(exxm-eyym)).^2+exym.^2);
                    if isfield(model,'exponent')
                        m=model.exponent;
                    else
                        m=2;
                    end
                    eeq=(max(0,e1).^m+max(0,e2).^m+nuo*max(0,-e1-e2).^m/(1-nuo)).^(1/m);
                case 'Marigo'
                    e1=0.5*(exxm+eyym)+sqrt((0.5*(exxm-eyym)).^2+exym.^2);
                    e2=0.5*(exxm+eyym)-sqrt((0.5*(exxm-eyym)).^2+exym.^2);
                    m=2;
                    eeq=(max(0,e1).^m+max(0,e2).^m+2*nuo*max(0,e1).*max(0,e2)).^(1/m);
                case 'vonMises'
                    
                    if isfield(model,'gamma')
                        g=model.gamma;
                    else
                        g=1;
                    end
                    nu=model.nu;
                    J1=exxm+eyym;
                    J2=(exxm.^2+eyym.^2-exxm.*eyym+3*(exym.^2));
                    eeq=(g-1)*J1/(2*g*(1-2*nu))+sqrt(((g-1)*J1/(1-2*nu)).^2+12*g*J2/((1+nu)^2))/(2*g);
            end
            eeqmax=max(eeqmaxo,eeq);
            switch model.law
                case 'linear'
                    d=max(0,min(model.a*eeqmax,model.dmax));
                case 'exponential'
                    d=(1-exp(-eeqmax*(1./model.b)))*model.a;
                case 'ulb'
                    ko=model.ko;
                    a=model.a;
                    b=model.b;
                    d=(1-ko*(1-a+a*exp(-b*(eeqmax-ko)))./eeqmax).*double(eeqmax>ko);
                otherwise
                    error ('INVALID DAMAGE LAW');
                    
            end
            D=max(0,1-d);
            D=diag(sparse(D));
            save(fullfile('TMP',sprintf('%d_matmod',nmod)),'D','eeqmax','-append');
        case 'elastic_plastic_homogeneous_isotropic'
            load(fullfile('TMP',sprintf('%d_matmod',nmod)),'model','E','S','Sp','Ep','Eeqp','Ui');
            
            DU=sparse(U-Ui)/nb_sub_cycles;
            
            for ic=1:nb_sub_cycles
                
                if dflag
                    meshfile=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,0));
                else
                    meshfile=fullfile('TMP',sprintf('%d_mesh_%d',nmod,0));
                end
                
                if nlflag
                    load(fullfile('TMP',sprintf('%d_matmod',nmod)),'R','Rp');
                    %                if keep
                    load(meshfile,'xo','yo','zo','Nnodes');
                    %                 else
                    %                     if restart
                    %                         load(meshfile,'xo','yo','zo','Nnodes');
                    %                     else
                    %                         load([meshfile,'-tmp'],'xo','yo','zo','Nnodes');
                    %                     end
                    %                 end
                    nn=prod(Nnodes);
                    Ux=0.5*DU((1:nn));
                    Uy=0.5*DU(nn+(1:nn));
                    if dflag
                        Uz=0.5*DU(2*nn+(1:nn));
                    else
                        Uz=zeros(size(zi));
                    end
                    xo=xo+Ux;
                    yo=yo+Uy;
                    zo=zo+Uz;
                    save(meshfile,'xo','yo','zo','-append');
                    if restart
                        save([meshfile,'-tmp'],'xo','yo','zo','Nnodes');
                    end
                    if dflag
                        CreateGradBasisFunction3D(1,nmod);
                    else
                        CreateGradBasisFunction(1,nmod);
                    end
                end
                
                restart=restart||nlflag;
                if restart
                    load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,iscale-1)),'epsyy');
                    load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,iscale-1)),'Uxy','Uyx');
                    load(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,iscale-1)),'epszz');
                    load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,iscale-1)),'Uxz','Uzx');
                    load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,iscale-1)),'Uyz','Uzy');
                    nmodo=nmod;
                    if isfield(param,'coupling_parameter')
                        load(fullfile('TMP','sample0'),'sizeim');
                        phi=CreateFiniteElementBasis3D(meshfile,sizeim,1,[],'Gauss_points');
                        load(meshfile,'selected');
                        notonboundary=sparse(((phi*selected)>0.9));
                    else
                        notonboundary=spones(size(epsxx,1),1);
                    end
                    %                     load(meshfile,'xo','yo','zo','Nnodes');
                    % figure
                    % plot3(xo,yo,zo,'kx')
                    % hold on
                    % plot3(phi*xo,phi*yo,phi*zo,'b+')
                    % xg=phi*xo;
                    % yg=phi*yo;
                    % zg=phi*zo;
                    % plot3(xg(~notonboundary),yg(~notonboundary),zg(~notonboundary),'ro')
                    % pause
                end
                ngauss=size(epsxx,1);
                %             if ~keep
                %                 if numel(Sp)==1
                %                 S=sparse(ngauss,6);
                %                 else
                %                 S=Sp;
                %                 end
                %                 if nlflag
                %                     R=Rp;
                %                 end
                %             end
                if numel(E)==1
                    E=sparse(ngauss,6);
                    S=sparse(ngauss,6);
                    Ep=sparse(ngauss,6);
                    Eeqp=sparse(ngauss,1);
                    if nlflag
                        R=repmat(sparse([1,0,0,0,1,0,0,0,1]),ngauss,1);
                    end
                end
                mu=model.mu;
                lambda=model.lambda;
                Delas=sparse([2*mu+lambda,lambda,lambda,0,0,0;...
                    lambda,2*mu+lambda,lambda,0,0,0;...
                    lambda,lambda,2*mu+lambda,0,0,0;...
                    0,0,0,mu,0,0;...
                    0,0,0,0,mu,0;...
                    0,0,0,0,0,mu...
                    ]');
                halfshear=sparse(1:6,1:6,[1,1,1,0.5,0.5,0.5],6,6);
                doubleshear=sparse(1:6,1:6,[1,1,1,2,2,2],6,6);
                
                
                DE=[epsxx*DU,epsyy*DU,epszz*DU,(Uxy+Uyx)*DU,(Uyz+Uzy)*DU,(Uxz+Uzx)*DU];
                E=E+DE;
                if nlflag
                    DR=0.5*[(Uxy-Uyx)*DU,(Uxz-Uzx)*DU,(Uyz-Uzy)*DU];
                    [S,R]=GloToLoc(S,R,DR);
                    [DE]=GloToLoc(DE*halfshear,R);
                    DE=DE*doubleshear;
                    [Ep]=GloToLoc(Ep*halfshear,R);
                    Ep=Ep*doubleshear;
                end
                Seq=GetEquivalentStress(model,S);
                S=S+DE*Delas;
                nSeq=GetEquivalentStress(model,S);
                [Sy,H]=ComputePlasticFlow(model,Eeqp);
                %             figure
                %             plot(Eeqp,Seq,'bx')
                %             hold on
                %             plot(Eeqp,nSeq,'ro')
                %             title('Seq and nSeq before plastic flow')
                %              figure
                %             plot(Eeqp,Sy,'x')
                %             title('Sy before plastic flow')
                
                if keep
                    doplasticres=(Seq>=Sy)&(nSeq>Seq)&notonboundary;
                else
                    doplasticres=(nSeq>=Sy)&(nSeq>Seq)&notonboundary;
                end
                %             figure
                %             plot(Seq,doplasticres,'x')
                
                
                Dmatx=repmat(Delas(:)',ngauss,1);
                if any(doplasticres)
                    %                 fdoplasticres=find(~doplasticres);
                    %                             figure
                    %             plot(Eeqp(fdoplasticres),Seq(fdoplasticres),'bx')
                    %             hold on
                    %             plot(Eeqp(fdoplasticres),nSeq(fdoplasticres),'ro')
                    %             title('Seq and nSeq no plastic flow')
                    %              figure
                    %             plot(Eeqp(fdoplasticres),Sy(fdoplasticres),'x')
                    %             title('Sy before no flow')
                    
                    
                    
                    found=find(doplasticres);
                    no=sparse(1:length(found),found,1,length(found),ngauss);
                    Dmat=no*Dmatx;
                    hmat=sparse([1,1,1,2,2,2,3,3,3,4,5,6],...
                        [1,2,3,1,2,3,1,2,3,4,5,6],...
                        [1,-0.5,-0.5,-0.5,1,-0.5,-0.5,-0.5,1,3,3,3],6,6);
                    hmat=repmat(hmat(:)',length(found),1);
                    [indj,indi]=meshgrid(1:6,1:6);
                    ondiag=find(indi(:)==indj(:));
                    se=no*S;
                    s=se;
                    eeqp=no*Eeqp;
                    db=sparse(length(found),1);
                    res=Inf;iter=1;
                    
                    %                 seq=GetEquivalentStress(model,s);
                    %                     [sy,h]=ComputePlasticFlow(model,eeqp);
                    %                       figure
                    %             plot(eeqp,seq,'bx')
                    %             hold on
                    %             plot(eeqp,sy,'r+')
                    %             title('seq sy before plastic flow')
                    
                    while res>1e-4&&iter<20
                        seq=GetEquivalentStress(model,s);
                        [sy,h]=ComputePlasticFlow(model,eeqp);
                        if iter==1
                            [gradf,dgradf,beta]=YieldSurfaceGradient(model,s,seq,h);
                            s=s-diag((seq-sy).*beta)*dgradf;
                        else
                            [gradf,dgradf]=YieldSurfaceGradient(model,s,seq,h);
                            smat=hmat-gradf(:,indi(:)).*gradf(:,indj(:));
                            smat=diag(sparse(1./seq))*smat;
                            qmat=sparse(1,1);
                            for kk=1:6
                                qmat=qmat+Dmat(:,indi(:)+(kk-1)*6).*smat(:,kk+(indj(:)-1)*6);
                            end
                            qmat=diag(db)*qmat;
                            qmat(:,ondiag)=qmat(:,ondiag)+1.;
                            
                            if iter==2
                                iqmat=qmat;
                            else
                                iqmat=mexInv66(full(qmat)');
                                iqmat=sparse(iqmat)';
                            end
                            rvec=s-se+diag(db)*dgradf;
                            
                            qvec=sparse(1,1);
                            for k=1:6
                                qvec=qvec+diag(rvec(:,k))*iqmat(:,indi(:,k)+6*(indj(:,k)-1));
                            end
                            
                            svec=sparse(1,1);
                            for k=1:6
                                svec=svec+diag(dgradf(:,k))*iqmat(:,indi(:,k)+6*(indj(:,k)-1));
                            end
                            
                            dnum=(seq-sy)-sum(gradf.*qvec,2);
                            dnom=h+sum(gradf.*svec,2);
                            ddbar=max(0,dnum./dnom);
                            s=s-qvec-diag(ddbar)*svec;
                            db=db+ddbar;
                            eeqp=eeqp+ddbar;
                        end
                        res=max(abs(seq-sy)./sy);
                        iter=iter+1;
                    end
                    if res>1e-4
                        warning('CONVERGENCE FAILURE IN THE PLASTIC FLOW RES=%f',full(res));
                    end
                    seq=GetEquivalentStress(model,s);
                    [sy,h]=ComputePlasticFlow(model,eeqp);
                    [gradf,dgradf,beta]=YieldSurfaceGradient(model,s,seq,h);
                    smat=hmat-gradf(:,indi(:)).*gradf(:,indj(:));
                    smat=diag(sparse(1./seq))*smat;
                    qmat=sparse(1,1);
                    for k=1:6
                        qmat=qmat+Dmat(:,indi(:)+(k-1)*6).*smat(:,k+(indj(:)-1)*6);
                    end
                    qmat=diag(db)*qmat;
                    qmat(:,ondiag)=qmat(:,ondiag)+1.;
                    if iter==2
                        iqmat=qmat;
                    else
                        iqmat=mexInv66(full(qmat)');
                        iqmat=sparse(iqmat)';
                    end
                    rmat=sparse(1,1);
                    for k=1:6
                        rmat=rmat+iqmat(:,indi(:)+(k-1)*6).*Dmat(:,k+(indj(:)-1)*6);
                    end
                    dmatx=rmat-diag(beta)*(dgradf(:,indi(:)).*dgradf(:,indj(:)));
                    if nlflag
                        [dmatx]=GloToLoc(dmatx,no*R,1);
                    end
                    on=sparse(found,1:length(found),1,ngauss,length(found));
                    Dmatx=Dmatx+on*(dmatx-Dmat);
                    S=S+on*(s-se);
                    dEp=eeqp-no*Eeqp;
                    Eeqp=Eeqp+on*dEp;
                    Ep=Ep+on*(diag(dEp)*gradf);
                    %                 figure
                    %                 plot(eeqp,seq,'bx')
                    %                 hold on
                    %                  plot(eeqp,sy,'r+')
                    % title('after plastic flow')
                    %                 figure
                    %
                    %             nSeq=GetEquivalentStress(model,S);
                    %                  plot(Eeqp,nSeq,'k*')
                    % title('all after plastic flow')
                    %  [syo,ho]=ComputePlasticFlow(model,(0:0.005:0.2)');
                    %  eo=syo(1)/model.young;
                    %  figure
                    %  plot(eo+(0:0.005:0.2)',syo,'r-')
                    %  hold on
                    %  plot([0,eo],[0,syo(1)],'b-')
                    %                      plot(eeqp,sy,'bo')
                    %            pause
                end
                
                if nlflag
                    [S]=GloToLoc(S,R,1);
                    [Ep]=GloToLoc(Ep*halfshear,R,1);
                    Ep=Ep*doubleshear;
                    clear Uxz Uzx Uxy Uyx Uyz Uzy epsxx epsyy epszz
                    save(fullfile('TMP',sprintf('%d_matmod',nmod)),'R','-append');
                end
            end
            %            if keep
            Ui=U;
            save(fullfile('TMP',sprintf('%d_matmod',nmod)),'Dmatx','E','S','Ep','Eeqp','DU','Ui','-append');
            %             else
            %                 save(fullfile('TMP',sprintf('%d_matmod',nmod)),'Dmatx','S','-append');
            %             end
            
            %                 if dflag
            %                     meshfile=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,0));
            %                 else
            %                     meshfile=fullfile('TMP',sprintf('%d_mesh_%d',nmod,0));
            %                 end
            %                         load(meshfile,'xo','yo','zo','Nnodes');
            %             phi=CreateFiniteElementBasis3D(meshfile,[100,100,100],1,[],'Gauss_points');
            %             xg=phi*xo;
            %             yg=phi*yo;
            %             zg=phi*zo;
            %
            %                 Seq=GetEquivalentStress(model,S);
            %                 figure
            %             scatter3(xg,yg,zg,10+0*zg,Seq)
            %             pause
            
    end
    if isfield(param,'interface_model')
        matmod=param.interface_model;
        switch matmod
            case 'cohesive'
                load(fullfile('TMP',sprintf('%d_intmod',nmod)),'model');
                if restart
                    load(fullfile('TMP',sprintf('%d_dphin_%d',nmod,iscale-1)),'dphincz');
                    load(fullfile('TMP',sprintf('%d_dphit_%d',nmod,iscale-1)),'dphitcz');
                end
                un=dphincz*U;
                ut=dphitcz*U;
                
                ueq=sign(un).*abs(pix2m*un+i*model.alpha*model.beta*pix2m*ut);
                % figure
                % plot(ueq)
                % title('ComputeNewInternalState')
                save(fullfile('TMP',sprintf('%d_intmod',nmod)),'ueq','-append');
                
                
            otherwise
                error ('INVALID INTERFACE MODEL');
                
        end
    end
    
end
