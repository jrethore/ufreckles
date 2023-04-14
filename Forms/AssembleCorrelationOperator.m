function AssembleCorrelationOperator(iscale,nmod,dflag)
if nargin<3, dflag=0;end
load(fullfile('TMP','params'));
onflight=0;
if isfield(param,'onflight')
    onflight=param.onflight;
end
if onflight
    global phidf wdetJ
end
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
tic();
if iscell(param0.reference_image)
    ncams=length(param0.reference_image);
else
    ncams=1;
end
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
lc=(min(sizeim)/10);
lm=lc;
if isfield(param0,'regularization_parameter')
    lm=param0.regularization_parameter/2^(iscale-1);
    if isfield(param0,'detect')&&iscale>1
        if param0.detect
            lm=0.1*lm;
        end
    end
end
reg_type='none';
if isfield(param0,'regularization_type')
    if strcmp(param0.regularization_type,'none')
        reg_type=param0.regularization_type;
    else
        if (iscale==1)
            reg_type=param0.regularization_type;
        else
            if strcmp(param0.regularization_type,'median')
                reg_type='median';
            else
            reg_type='tiko';
            end
        end
    end
end
if (strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri'))&&(iscale>1)
    param.basis='fem';
end
disp(sprintf('    Correlation operator...'));
for ncam=1:ncams
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
    if ~onflight
        load(fullfile('TMP',sprintf('%d_phidf_%d',nmod*10^(ncam-1),iscale-1)),'phidf');
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(ncam-1),iscale-1)),'wdetJ');
    end
    M=phidf'*mask*wdetJ*phidf;
    nM=norm(diag(M));
    R=sparse(size(M,1),size(M,2));
    P=sparse(size(M,1),size(M,2));
    if strcmp(param.basis,'fem')&&~isfield(param,'gp_file')
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(ncam-1),iscale-1)),'xo','yo','selected','Nnodes');
        load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
        ki=1/lc;
        V=[repmat(cos(2*pi*(xo)*ki).*cos(2*pi*(yo)*ki),2,1);zeros(size(M,1)-2*prod(Nnodes),1)];
        if mean(Nnodes(1:2))<2
            reg_type='none';
        end
        switch reg_type
            case 'tiko'
                %                 if iscale==1
                %                     keyboard
                %                 end
                %if iscale>1
                if ~((iscale==1)&&isfield(param,'coupling_parameter'))
                    Rall=AssembleRegularizationOperator([],iscale,nmod,2);
                    a=(V'*M*V)/(V'*Rall*V);
                    a=a*(2*lm/lc)^2;
                    R=R+a*Rall;
                    disp(sprintf('    Tikhonov regularization...'));
                else
                    if dflag
                        Nddl_tot=size(M,1);
                        save(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,iscale-1)),'Nddl_tot');
                        
                        Ri=AssembleRegularizationOperator3D(nmod,1);
                        
                    else
                        selected_nodes=find(~selected(:));
                        Ri=AssembleRegularizationOperator(selected_nodes,iscale,nmod,2);
                        
                        disp(sprintf('    Tikhonov boundary regularization...'));
                    end
                    a=(V'*M*V)/(V'*Ri*V);
                    a=0.1*a*(2*lm/lc)^2;
                    R=R+a*Ri;
                end
                load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'cut_nodes');
                %     if numel(cut_nodes)
                %         Ri=AssembleRegularizationForMaskedElements(cut_nodes,iscale,nmod);
                %         R=R+10*nM*Ri/norm(diag(Ri));
                %     end
                
                if (iscale==1)&&isfield(param,'enrichment')&&(~isfield(param,'coupling_parameter'))
                    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(ncam-1),iscale-1)),'Nnodes');
                    CreateGradBasisFunction(iscale,nmod*10^(ncam-1));
                    load(fullfile('TMP',sprintf('%d_xepsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_xepsyy_%d',nmod,iscale-1)),'epsyy');
                    load(fullfile('TMP',sprintf('%d_xepsxy_%d',nmod,iscale-1)),'Uxy','Uyx');
                    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'masks');
                    Nddl_tot=size(epsxx,2);
                    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
                    if isempty(unmasked_nodes)
                        Nddls=2*prod(Nnodes);
                    else
                        Nddls=2*numel(unmasked_nodes);
                    end
                    Ri=epsxx'*masks*wdetJ*epsxx+epsyy'*masks*wdetJ*epsyy+Uxy'*masks*wdetJ*Uxy+Uyx'*masks*wdetJ*Uyx;
                    
                    [indi,indj,val]=find(Ri);
                    keep=~((indi<=Nddls)&(indj<=Nddls));
                    Ri=sparse(indi,indj,val.*keep,Nddl_tot,Nddl_tot);
                    %                    R=R+Po*Ri*Po;
                    R=R+a*Ri;
                    %       if ~iscell(param0.levelset_file)
                    %             nbfis=1;
                    %         else
                    %             nbfis=numel(param0.levelset_file);
                    %         end
                    %         for ic=1:nbfis
                    %             load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'face_nodes');
                    %             Ri=AssembleRegularizationForMaskedElements(face_nodes,iscale,nmod);
                    %             size(Ri)
                    % display(sprintf('Assembling penlaty terms for enriched nodes...'));
                    % %ind=[face_nodes,face_nodes+prod(Nnodes)];
                    % %Ri=sparse(ind,ind,1,size(M,1),size(M,2));
                    %
                    %
                    %
                    %
                    %
                    % %            R=R+1.e1*nM*Ri/norm(diag(Ri));
                    %         end
                end
                
                
                
            case {'equilibrium_gap','constitutive_gap'}
                if (param0.detect==0)
                    
                    selected_nodes=find(~selected(:));
                    Ri=AssembleRegularizationOperator(selected_nodes,iscale,nmod,2);
                    disp(sprintf('    Tikhonov boundary regularization...'));
                    %                                  if (norm(diag(Ri)))>0
                    %                                      R=R+1*nM*Ri/norm(diag(Ri));
                    %                                  end
                    
                    S=diag(sparse([repmat(selected(:),2,1);ones(size(M,1)-2*prod(Nnodes),1)]));
                    a=((S*V)'*M*(S*V))/(V'*Ri*V);
                    a=a*(0.1*2*lm/lc)^2;
                    R=R+a*Ri;
                end
                if ~isfield(param,'coupling_parameter')
                    wdetJo=wdetJ;
                    CreateGradBasisFunction(iscale,nmod*10^(ncam-1));
                    
                    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
                    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(ncam-1),iscale-1)),'Nnodes','selected','conn');
                    load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
                    load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'epsxy');
                    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'maskg');
                    
                    %                      mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
                    %  [phig]=CreateFiniteElementBasis(mesh_file,1,1,[],'Gauss_points');
                    %                             outer_nodes=~selected(:);
                    %  fouter_nodes=AddOneNodeLayer(conn,find(outer_nodes));
                    %  outer_nodes(fouter_nodes)=1;
                    %               maskg=diag(sparse((phig*(1-outer_nodes))>0.1));
                    
                    if isfield(param,'material_model')
                        [sxx,syy,sxy]=ComputeTangentStress(epsxx,epsyy,epsxy,nmod);
                    else
                        nuo=0.3;
                        mu=0.5/(1+nuo);
                        lambda = nuo/((1.-nuo^2));%CP
                        ltr=lambda*(epsxx+epsyy);
                        sxx=2*mu*epsxx+ltr;
                        syy=2*mu*epsyy+ltr;
                        sxy=mu*epsxy;
                    end
                    Ri= epsxx'*maskg*(wdetJ*sxx)+epsyy'*maskg*(wdetJ*syy)+epsxy'*maskg*(wdetJ*sxy);
                    if isfield(param,'enrichment')&&(iscale==1)
                        load(fullfile('TMP',sprintf('%d_xepsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
                        load(fullfile('TMP',sprintf('%d_xepsyy_%d',nmod,iscale-1)),'epsyy');
                        load(fullfile('TMP',sprintf('%d_xepsxy_%d',nmod,iscale-1)),'epsxy');
                        load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'masks');
                        Nddl_tot=size(epsxx,2);
                        %             phis=CreateFiniteElementBasis(mesh_file,1,1,enriched_nodes,'sub_cells',true);
                        %               masks=diag(sparse((phis*(1-outer_nodes))>0.0001));
                        if isfield(param,'material_model')
                            [sxx,syy,sxy]=ComputeTangentStress(epsxx,epsyy,epsxy,nmod);
                        else
                            ltr=lambda*(epsxx+epsyy);
                            sxx=2*mu*epsxx+ltr;
                            syy=2*mu*epsyy+ltr;
                            sxy=mu*epsxy;
                        end
                        Rx= epsxx'*masks*wdetJ*sxx+epsyy'*masks*wdetJ*syy+epsxy'*masks*wdetJ*sxy;
                        [indi,indj,val]=find(Rx);
                        Nddls=size(Ri,1);
                        keep=~((indi<=Nddls)&(indj<=Nddls));
                        Rx=sparse(indi,indj,val.*keep,Nddl_tot,Nddl_tot);
                        if Nddl_tot>size(Ri,2)
                            Ro=sparse(Nddl_tot-size(Ri,1),Nddl_tot-size(Ri,2));
                            Ri=blkdiag(Ri,Ro);
                        end
                        Ri=Ri+Rx;
                        
                    end
                    select=ones(size(Ri,1),1);
                    outer_nodes=~selected;
                    Nnods=prod(Nnodes);
                    if ~isempty(unmasked_nodes)
                        outer_nodes=outer_nodes(unmasked_nodes);
                        Nnods=length(unmasked_nodes);
                    end
                    reject=find(outer_nodes(:));
                    select(reject)=0;
                    select(prod(Nnodes)+reject)=0;
                    if isfield(param,'enrichment')&&(iscale==1)
                        if ~iscell(param0.levelset_file)
                            nbfis=1;
                        else
                            nbfis=numel(param0.levelset_file);
                        end
                        for ic=1:nbfis
                            
                            load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'face_nodes');
                            select(face_nodes)=1;
                            select(Nnods+face_nodes)=1;
                        end
                    end
                    select=diag(sparse(select));
                    switch reg_type
                        case 'equilibrium_gap'
                            disp(sprintf('    Equilibrium Gap regularization...'));
                            Ri=select*Ri;
                            Ri=Ri'*Ri;
                            Pi=sparse(size(M,1),size(M,2));
                            a=(V'*M*V)/(V'*Ri*V);
                            a=a*(2*lm/lc)^4;
                        case 'constitutive_gap'
                            disp(sprintf('    Constitutive Equation Gap regularization...'));
                            Pi=1e9*diag(diag(Ri).*(~diag(select)));
                            a=(V'*M*V)/(V'*Ri*V);
                            a=a*(2*lm/lc)^2;
                            
                    end
                    R=R+a*Ri;
                    P=P+a*Pi;
                    %                    R=R+Po*Ri*Po;
                    
                    wdetJ=wdetJo;
                end
            case 'none'
                disp(sprintf('    No regularization...'));
            case 'median'
                disp(sprintf('    Median regularization...'));
                
            otherwise
                error('REGULARIZATION TYPE NOT CODED YET');
        end
    end
    if strcmp(param.basis,'vic-nurbs')
        if param.Nelems>1
            load(fullfile('TMP',sprintf('sample0_%d',1-1)),'ls1');
            load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim','nband','on');
            ls1=ls1(nband(on));
            load(fullfile('TMP',sprintf('%d_phi_%d',nmod,iscale-1)),'dphii','phii');
            load(fullfile('TMP',sprintf('%d_vicmesh_%d',nmod,iscale-1)),'mesh_size');
            if lm==lc
                lm=mean(mesh_size);
            end
            Ri=dphii'*dphii;
            ki=1/lc;
            vi=cos(2*pi*ls1*ki);
            L=phii'*wdetJ*phii;
            b=phii'*wdetJ*vi;
            V=L\b;
            a=(V'*M*V)/(V'*Ri*V);
            a=a*(2*lm/lc)^2;
            R=R+a*Ri;
            
            %             val=(1*max(diag(M))-mean(diag(M)))/mean(diag(Ri));
            %             %        Po=diag(sparse(ones(size(M,1),1)*sqrt(val)));
            %             Po=sqrt(val);
            %             R=R+penalty*Po*Ri*Po;
            disp(sprintf('    Tikhonov regularization...'));
        end
        
    end
    if strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri')
        switch reg_type
            case 'tiko'
                disp(sprintf('    Tikhonov regularization...'));
                load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(ncam-1),iscale-1)),'xo','yo','Nnodes');
                load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'sizeim');
                load(fullfile('TMP',sprintf('%d_phio_%d',nmod,iscale-1)),'phio');
                ki=1/lc;
                Vn=cos(2*pi*(xo)*ki).*cos(2*pi*(yo)*ki);
                L=phio'*phio;
                b=phio'*Vn;
                V=L\b;
                V=repmat(V,2,1);
                
                CreateGradBasisFunction(iscale,nmod);
                
                load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
                load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uxy','Uyx');
                Ri= epsxx'*wdetJ*epsxx+Uxy'*wdetJ*Uxy...
                    +epsyy'*wdetJ*epsyy+Uyx'*wdetJ*Uyx;
                a=(V'*M*V)/(V'*Ri*V);
                a=a*(2*lm/lc)^2;
                R=R+a*Ri;
                
            case 'none'
                disp(sprintf('    No regularization...'));
                
            otherwise
                error('REGULARIZATION TYPE NOT CODED YET');
        end
        
        
    end
    save(fullfile('TMP',sprintf('%d_operator_%d',nmod*10^(ncam-1),iscale-1)),'M','R','P');
    
end
disp(sprintf('Computing correlation operator for model %d...%6.2f s',nmod,toc()));



end
