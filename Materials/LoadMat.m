function LoadMat(nmod)

load(fullfile('TMP','params'));
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if isfield(param,'integration')
    integration=param.integration;
else
    integration='Gauss_points';
end
iscale=1;

if isfield(param,'material_model')
    matmod=param.material_model;
else
    matmod='elastic_homogeneous_isotropic';
disp(sprintf('    Fictitious material...'));
                    nuo=0.25;
           param.material_model=matmod;         
    param.material_parameters.mu=0.5/(1+nuo);
    param.material_parameters.lambda=nuo/((1.+nuo)*(1.-2.*nuo));%DP ou 3D
save(fullfile('TMP',sprintf('%d_params',nmod)),'param','-append');

end

tic;
iscale=1;
roi=param0.roi;

switch matmod
    case 'aline'
    case 'hyperelastic_homogeneous_isotropic'
    case 'elastic_homogeneous_isotropic'
        if isfield(param,'material_file')
            filemat=param.material_file;
            load(filemat,'mu','lambda');
        else
            model=param.material_parameters;
        if isfield(model,'young')
            model.mu=0.5*model.young/(1+model.nu);
            model.lambda=model.nu*model.young/((1.+model.nu)*(1.-2.*model.nu));
            if isfield(param,'plane_stress')
                if param.plane_stress
                     model.lambda=model.nu*model.young/((1+model.nu)*(1-model.nu));%CP
                end
            end
%      model.lambda=model.nu*model.young/((1+model.nu)*(1-model.nu));%CP
%      'CP'
        end
            mu=model.mu;
            lambda=model.lambda;
        end
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambda');
        if isfield(model,'lc')
            lc=model.lc;
            save(fullfile('TMP',sprintf('%d_matmod',nmod)),'lc','-append');
        end
    case 'elastic_homogeneous_isotropic_damage'
        model=param.material_parameters;
        if isfield(model,'young')
            model.mu=0.5*model.young/(1+model.nu);
            model.lambda=model.nu*model.young/((1.+model.nu)*(1.-2.*model.nu));
        end
        D=1;DD=1;eeqmaxo=0;
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'model','D','DD','eeqmaxo','-v7.3');
    case 'elastic_plastic_homogeneous_isotropic'
        model=param.material_parameters;
        Eo=model.young;
        nuo=model.nu;

        mu=0.5*Eo/(1+nuo);
        lambda = nuo*Eo/((1.+nuo)*(1.-2.*nuo));
        model.mu=mu;
        model.lambda=lambda;
        const=Eo*(1-nuo)/((1+nuo)*(1-2*nuo));
        Eelast(1)=const;
        Eelast(2)=const*nuo/(1-nuo);
        Eelast(3)=const*(1-2*nuo)/(2*(1-nuo));
        model.Eelast=Eelast;
        Dmatx=sparse([2*mu+lambda,lambda,lambda,0,0,0,...
            lambda,2*mu+lambda,lambda,0,0,0,...
            lambda,lambda,2*mu+lambda,0,0,0,...
            0,0,0,mu,0,0,...
            0,0,0,0,mu,0,...
            0,0,0,0,0,mu...
            ]);
        E=0;S=0;Sp=0;Ui=0;Ep=0;Eeqp=0;R=0;Rp=0;
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'model','Dmatx','E','S','Sp','Rp','Ep','Eeqp','Ui','R','-v7.3');
    case 'elastic_heterogeneous_isotropic'
        load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
        filemat=param.material_file;
        load(filemat,'mu','lambda');
        mu=mu(roi(1):roi(2),roi(3):roi(4));
        lambda=lambda(roi(1):roi(2),roi(3):roi(4));
        if strcmp(param.basis,'fem')
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
            phi=CreateFiniteElementBasis(mesh_file,sizeim,1,unmasked_nodes,integration);
            load(mesh_file,'xo','yo');
            [Yo,Xo]=meshgrid(yo-0.5,xo-0.5);
            if isempty(unmasked_nodes)
                Xo=Xo(:);Yo=Yo(:);
            else
                Xo=Xo(unmasked_nodes);Yo=Yo(unmasked_nodes);
            end
            Xg=phi*Xo;
            Yg=phi*Yo;
            lambda=interp2(lambda,Yg,Xg,'*linear');
            mu=interp2(mu,Yg,Xg,'*linear');
            %             lambda=interp2(lambda,Yo,Xo,'*linear');
            %             mu=interp2(mu,Yo,Xo,'*linear');
            %                          if ~isempty(unmasked_nodes)
            %                              lambda=lambda(unmasked_nodes);
            %                             mu=mu(unmasked_nodes);
            %                          end
            %             lambda=phi*lambda(:);
            %             mu=phi*mu(:);
        end
        mu=diag(sparse(mu(:)));
        lambda=diag(sparse(lambda(:)));
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambda');

    case 'elastic_heterogeneous_orthotropic'
        filemat=param.material_file;
        load(filemat,'mu','lambdaxx','lambdaxy','lambdayx','lambdayy');
        mu=mu(roi(1):roi(2),roi(3):roi(4));
        lambdaxx=lambdaxx(roi(1):roi(2),roi(3):roi(4));
        lambdaxy=lambdaxy(roi(1):roi(2),roi(3):roi(4));
        lambdayx=lambdayx(roi(1):roi(2),roi(3):roi(4));
        lambdayy=lambdayy(roi(1):roi(2),roi(3):roi(4));
        if strcmp(param.basis,'fem')
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
            phi=CreateFiniteElementBasis(mesh_file,sizeim,1,unmasked_nodes,integration);
            load(sprintf(mesh_file,nmod,iscale-1),'xo','yo');
            [Yo,Xo]=meshgrid(max(yo-0.5,1),max(xo-0.5,1));
            if isempty(unmasked_nodes)
                Xo=Xo(:);Yo=Yo(:);
            else
                Xo=Xo(unmasked_nodes);Yo=Yo(unmasked_nodes);
            end
            Xg=phi*Xo;
            Yg=phi*Yo;
            lambdaxx=interp2(lambdaxx,Yg,Xg,'*linear');
            lambdaxy=interp2(lambdaxy,Yg,Xg,'*linear');
            lambdayx=interp2(lambdayx,Yg,Xg,'*linear');
            lambdayy=interp2(lambdayy,Yg,Xg,'*linear');
            mu=interp2(mu,Yg,Xg,'*linear');
        end
        mu=diag(sparse(mu(:)));
        lambdaxx=diag(sparse(lambdaxx(:)));
        lambdaxy=diag(sparse(lambdaxy(:)));
        lambdayx=diag(sparse(lambdayx(:)));
        lambdayy=diag(sparse(lambdayy(:)));
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambdaxx','lambdaxy','lambdayx','lambdayy','-v7.3');
    case 'elastic_heterogeneous_anisotropic'
        filemat=param.material_file;
        load(filemat,'C');%xx,yy,zz,xy,yz,xz
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'C');
        if isfield(param,'enrichment')
        load(filemat,'Cs');%xx,yy,zz,xy,yz,xz
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'Cs','-append');
        end
    case 'elastic_homogeneous_orthotropic'
            model=param.material_parameters;
            mu=model.mu;
            lambdaxx=model.lambdaxx;
            lambdayy=model.lambdayy;
            lambdaxy=model.lambdaxy;
            lambdayx=model.lambdayx;

        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambdaxx','lambdaxy','lambdayx','lambdayy');

    case 'hyperelastic_heterogeneous_isotropic'
        model=param.material_parameters;
        filemat=model.section_file;
        sec=readim(filemat);
        sec=double(sec(roi(1):roi(2),roi(3):roi(4)));

        load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix','inde','wdetJ','sizeim');
        phix=phix(:,1:(size(phix,2)/2));
        mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
        load(sprintf(mesh_file,nmod,iscale-1),'conn','xo','yo');

        Xi=phix*xo;
        Yi=phix*yo;
        sec=interp2(sec,Yi(:),Xi(:),'*linear');
        M=phix'*wdetJ*phix;
        F=phix'*wdetJ*double(sec(:));
        Sn=M\F;
        Sn-255-Sn;
        Sn=Sn/mean(Sn);
        Se=1+0*sum(Sn(conn),2)./sum(conn>0,2);
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'Se','Sn','-v7.3');

    otherwise
        error ('INVALID MATERIAL MODEL');

end


if isfield(param,'interface_model')
    matmod=param.interface_model;
    switch matmod
        case {'cohesive','plastic'}
            model=param.interface_parameters;
            ueqmax=-1;
            ueq=0;
            save(fullfile('TMP',sprintf('%d_intmod',nmod)),'model','ueq','ueqmax','-v7.3');

        otherwise
            error ('INVALID INTERFACE MODEL');

    end

end



%disp(sprintf('Loading mat for model %d',nmod));


end
