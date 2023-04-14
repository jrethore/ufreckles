function [sxx,syy,sxy,sxz,syz,szz]=ComputeStress(exx,eyy,exy,nmod,exz,eyz,ezz)
if nargin<5, exz=sparse(size(exx,1),size(exx,2));eyz=exz;ezz=exz;end

load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if ~isfield(param,'material_model')
    matmod='elastic_homogeneous_isotropic';
else
    matmod=param.material_model;
end

switch matmod
    case 'elastic_homogeneous_isotropic'
        if isfield(param,'mu')
        mu=param.mu;
        lambda=param.lambda;
        else
        load(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambda');
        end
        ltr=lambda*(exx+eyy+ezz);
        sxx=2*mu*exx+ltr;
        syy=2*mu*eyy+ltr;
        szz=2*mu*ezz+ltr;
        sxy=mu*exy;
        sxz=mu*exz;
        syz=mu*eyz;
    case 'elastic_homogeneous_isotropic_damage'
        load(fullfile('TMP',sprintf('%d_matmod',nmod)),'model','D');
        mu=D*(model.mu);
        lambda=D*(model.lambda);
        ltr=lambda*(exx+eyy+ezz);
        sxx=2*mu*exx+ltr;
        syy=2*mu*eyy+ltr;
        szz=2*mu*ezz+ltr;
        sxy=mu*exy;
        sxz=mu*exz;
        syz=mu*eyz;
        S=[sxx,syy,szz,sxy,syz,sxz];
        save(fullfile('TMP',sprintf('%d_matmod',nmod)),'S','-append');
        case 'elastic_plastic_homogeneous_isotropic'
        load(fullfile('TMP',sprintf('%d_matmod',nmod)),'S');
        if numel(S)==1
            S=sparse(size(exx,1),6);
        end
        sxx=S(:,1);syy=S(:,2);szz=S(:,3);
        sxy=S(:,4);syz=S(:,5);sxz=S(:,6);
    case 'elastic_heterogeneous_isotropic'
        load(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambda');
        l2m=lambda+2*mu;
        sxx=l2m*exx+lambda*eyy;
        syy=lambda*exx+l2m*eyy;
        sxy=mu*exy;
    case {'elastic_heterogeneous_orthotropic','elastic_homogeneous_orthotropic'}
        load(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambdaxx','lambdaxy','lambdayx','lambdayy');
%         sxx=(lambdaxx+2*mu)*exx+lambdaxy*eyy;
%         syy=lambdayx*exx+(lambdayy+2*mu)*eyy;
%         sxy=mu*exy;
        sxx=(lambdaxx)*exx+lambdaxy*eyy;
        syy=lambdayx*exx+(lambdayy)*eyy;
        sxy=mu*exy;
        szz=sparse(size(ezz,1),size(ezz,2));
        sxz=sparse(size(ezz,1),size(ezz,2));
        syz=sparse(size(ezz,1),size(ezz,2));
    case 'elastic_heterogeneous_anisotropic'
                  load(fullfile('TMP',sprintf('%d_matmod',nmod)),'C');
if ~size(C,1)==size(exx,1)
              load(fullfile('TMP',sprintf('%d_matmod',nmod)),'Cs');
              C=Cs;
end

%xx,yy,zz,xy,yz,xz
        sxx=diag(C(:,1+6*(1-1)))*exx+diag(C(:,2+6*(1-1)))*eyy+diag(C(:,3+6*(1-1)))*ezz+...
            diag(C(:,4+6*(1-1)))*exy+diag(C(:,5+6*(1-1)))*eyz+diag(C(:,6+6*(1-1)))*exz;
        syy=diag(C(:,1+6*(2-1)))*exx+diag(C(:,2+6*(2-1)))*eyy+diag(C(:,3+6*(2-1)))*ezz+...
            diag(C(:,4+6*(2-1)))*exy+diag(C(:,5+6*(2-1)))*eyz+diag(C(:,6+6*(2-1)))*exz;
        szz=diag(C(:,1+6*(3-1)))*exx+diag(C(:,2+6*(3-1)))*eyy+diag(C(:,3+6*(3-1)))*ezz+...
            diag(C(:,4+6*(3-1)))*exy+diag(C(:,5+6*(3-1)))*eyz+diag(C(:,6+6*(3-1)))*exz;
        sxy=diag(C(:,1+6*(4-1)))*exx+diag(C(:,2+6*(4-1)))*eyy+diag(C(:,3+6*(4-1)))*ezz+...
            diag(C(:,4+6*(4-1)))*exy+diag(C(:,5+6*(4-1)))*eyz+diag(C(:,6+6*(4-1)))*exz;
        syz=diag(C(:,1+6*(5-1)))*exx+diag(C(:,2+6*(5-1)))*eyy+diag(C(:,3+6*(5-1)))*ezz+...
            diag(C(:,4+6*(5-1)))*exy+diag(C(:,5+6*(5-1)))*eyz+diag(C(:,6+6*(5-1)))*exz;
        sxz=diag(C(:,1+6*(6-1)))*exx+diag(C(:,2+6*(6-1)))*eyy+diag(C(:,3+6*(6-1)))*ezz+...
            diag(C(:,4+6*(6-1)))*exy+diag(C(:,5+6*(6-1)))*eyz+diag(C(:,6+6*(6-1)))*exz;

    otherwise
        error ('INVALID MATERIAL MODEL');

end




end