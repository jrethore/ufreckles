function UpdateInternalState(nmod,U,step)
if nargin<3, step=[];end
persistent nmodo epsxx epsyy epsxy epsxz epsyz epszz dphincz dphitcz
load(fullfile('TMP','params'),'param');
param0=param;
helper='none';
if isfield(param,'helper')
    helper=param.helper;
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
    restart=restart||nlflag;
    switch matmod
        
        case {'elastic_homogeneous_isotropic','elastic_heterogeneous_anisotropic','elastic_heterogeneous_orthotropic','elastic_homogeneous_orthotropic'}
            if ~isempty(step)
                if restart
                    
                    
                    
                    if dflag
                        load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
                        load(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,iscale-1)),'epsyy');
                        load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,iscale-1)),'epsxy');
                        load(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,iscale-1)),'epszz');
                        load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,iscale-1)),'epsxz');
                        load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,iscale-1)),'epsyz');
                        
                    else
                        
                        load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,iscale-1)),'epsxx');
                        load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,iscale-1)),'epsyy');
                        load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,iscale-1)),'epsxy');
                        epsxz=sparse(size(epsxx,1),size(epsxx,2));epsyz=epsxz;epszz=epsxz;
                    end
                    if isfield(param,'enrichment')&&(iscale==1)
                        zeps=sparse(size(epsxx,1),size(U,1)-size(epsxx,2));
                        epsxx=[epsxx,zeps];
                        epsyy=[epsyy,zeps];
                        epsxy=[epsxy,zeps];
                        epszz=[epszz,zeps];
                        epsxz=[epsxz,zeps];
                        epsyz=[epsyz,zeps];
                        
                    end
                end
                [sxx,syy,sxy,sxz,syz,szz]=ComputeStress(epsxx*U,epsyy*U,epsxy*U,nmod,epsxz*U,epsyz*U,epszz*U);
                S=[sxx,syy,szz,sxy,syz,sxz];
                save(sprintf('%s_%04d',strrep(param0.result_file,'.res',''),step),'S','-v7.3');
            end
        case 'elastic_homogeneous_isotropic_damage'
            load(fullfile('TMP',sprintf('%d_matmod',nmod)),'model','eeqmaxo');
            if restart
                load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,iscale-1)),'epsxx');
                load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,iscale-1)),'epsyy');
                load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,iscale-1)),'epsxy');
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
            eeqmaxo=max(eeqmaxo,eeq);
            switch model.law
                case 'linear'
                    d=max(0,min(model.a*eeqmaxo,model.dmax));
                    dd=(model.a).*(eeqmaxo).*(eeqmaxo<model.dmax);
                case 'exponential'
                    d=(1-exp(-eeqmaxo*(1./model.b)))*model.a;
                    %                dd=diag(sparse(eeqmaxo))*(exp(-eeqmaxo*(1./model.b))*diag(1./model.b))*model.a;
                    %                 figure
                    %                 plot(eeqmaxo,1-d,'bo')
                    %                 hold on
                    %                 plot(eeqmaxo,1-d-dd,'rx')
                case 'ulb'
                    ko=model.ko;
                    a=model.a;
                    b=model.b;
                    d=(1-ko*(1-a+a*exp(-b*(eeqmaxo-ko)))./eeqmaxo).*double(eeqmaxo>ko);
                otherwise
                    error ('INVALID DAMAGE LAW');
                    
            end
            D=max(0,1-d);
            %        DD=1-d-dd.*(eeq>eeqmaxo);
            D=diag(sparse(D));
            %        DD=diag(sparse(DD));
            DD=1;
            eeqmax=eeqmaxo;
            save(fullfile('TMP',sprintf('%d_matmod',nmod)),'D','DD','eeqmaxo','eeqmax','-append');
            if ~isempty(step)
                load(fullfile('TMP',sprintf('%d_matmod',nmod)),'S');
                save(sprintf('%s_%04d',param0.result_file,step),'model','d','eeqmax','S','-v7.3');
            end
        case 'elastic_plastic_homogeneous_isotropic'
            if ~isempty(step)
                load(fullfile('TMP',sprintf('%d_matmod',nmod)),'model','E','S','Ep','Eeqp');
                save(sprintf('%s_%04d',param0.result_file,step),'model','E','S','Ep','Eeqp','-v7.3');
            end
    end
    if isfield(param,'interface_model')
        matmod=param.interface_model;
        switch matmod
            case 'cohesive'
                load(fullfile('TMP',sprintf('%d_intmod',nmod)),'model','ueqmax');
                if restart
                    load(fullfile('TMP',sprintf('%d_dphin_%d',nmod,iscale-1)),'dphincz');
                    load(fullfile('TMP',sprintf('%d_dphit_%d',nmod,iscale-1)),'dphitcz');
                end
                un=dphincz*U;
                ut=dphitcz*U;
                
                ueq=sign(un).*abs(pix2m*un+i*model.alpha*model.beta*pix2m*ut);
                ueqmax=max(ueq,ueqmax);
                % figure
                % plot(ueq,'ro');
                % hold on
                % plot(ueqmax,'bx')
                % title('update internal stress')
                save(fullfile('TMP',sprintf('%d_intmod',nmod)),'ueq','ueqmax','-append');
            otherwise
                error ('INVALID INTERFACE MODEL');
                
        end
    end
end

end
