function run_fea_job(param,model)
clear functions
nmod=0;
param.onflight=1;
sizeim=param.roi([2,4])-param.roi([1,3])+1;
model.nscale=1;
LoadParameters(param);
LoadParameters(model,nmod);
save(fullfile('TMP','sample0_0.mat'),'sizeim')
save(fullfile('TMP','sample0.mat'),'sizeim')
LoadMask(nmod);
LoadMeshes(nmod);

filres=param.result_file;
%%
load(fullfile('TMP','0_mesh_0.mat'),'xo','yo','Nnodes')
extract=0;
selected=[];
load(strrep(param.result_file,'.res','.dat'),'-mat','selected')
if isempty(selected)

Fo=zeros(2*prod(Nnodes),1);
fromdic=0;
    indi=[];Up=[];
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        switch zone{4}
            case 6
                loads=zone{7};
                nodes=zone{6};
                for ii=1:2
                    Fo(nodes+(ii-1)*prod(Nnodes))=loads(ii,2)/length(nodes);
                    up=loads(ii,1);
                    if ~isnan(up)
                        indi=[indi;nodes+(ii-1)*prod(Nnodes)];
                        Up(end+(1:length(nodes)))=up;
                    end
                end
            case 5
                if (zone{8}>0)
                    extract=1;
                end
        end
    end
    Up=Up(:);
else
    fromdic=1;
    load(strrep(param.result_file,'.res','.dat'),'-mat','U')
    indi=find(~selected);
    indi=[indi(:);indi(:)+prod(Nnodes)];
    Up=U(indi,:);
    
Fo=zeros(2*prod(Nnodes),size(Up,2));
    clear U
end
C=sparse(indi,1:length(indi),1,2*prod(Nnodes),length(indi));
L=1;
Nddl=2*prod(Nnodes);
tips=zeros(size(model.zone,2),2);
cracks=[];
if extract
    roi=param.roi;
    freenodes=ones(length(xo),1);
    xo=xo+roi(1)-1;
    yo=yo+roi(3)-1;
    L=zeros(prod(Nnodes),1);
    modes=1:2;
    harms=0:7;
    nu=model.material_parameters.nu;
    kappa=(3-4*nu);
    zo=xo+1i*yo;
    nwi=length(modes)*length(harms);
    nw=0;
    ntip=0;
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        if zone{4}==5
            if (zone{8}>0)
                tips(iz,1)=ntip+1;
                indc=zone{10};
                cnodes=indc{1};
                ztip=zo(cnodes(1,1));
                tt=-diff(zo(cnodes(1:2,1)));
                tt=tt/abs(tt);
                rtip=mean(abs(zo(indc{2})-ztip));
                in=abs(zo-ztip)<=rtip;
                in(indc{2})=1;
                freenodes(in)=0;
                zzone=(zo(in)-ztip)*tt';
                idin=1:prod(Nnodes);
                idin=idin(in);
                for ip=1:length(zzone)
                    if any(cnodes(:,1)==idin(ip))
                        zzone(ip)=abs(zzone(ip))*exp(1i*pi);
                    end
                    if any(cnodes(:,2)==idin(ip))
                        zzone(ip)=abs(zzone(ip))*exp(-1i*pi);
                    end
                end
                phi=Williams(zzone,modes,harms,kappa);
                phi=phi*tt;
                L(in,nw+(1:size(phi,2)))=phi;
                ntip=ntip+1;
                nw=nw+nwi;
                if  (zone{9}>0)
                    tips(iz,2)=ntip+1;
                    ztip=zo(cnodes(end,1));
                    tt=diff(zo(cnodes(end+(-1:0),1)));
                    tt=tt/abs(tt);
                    rtip=mean(abs(zo(indc{3})-ztip));
                    in=abs(zo-ztip)<=rtip;
                    in(indc{3})=1;
                    freenodes(in)=0;
                    zzone=(zo(in)-ztip)*tt';
                    idin=1:prod(Nnodes);
                    idin=idin(in);
                    for ip=1:length(zzone)
                        if any(cnodes(:,1)==idin(ip))
                            zzone(ip)=abs(zzone(ip))*exp(-1i*pi);
                        end
                        if any(cnodes(:,2)==idin(ip))
                            zzone(ip)=abs(zzone(ip))*exp(1i*pi);
                        end
                    end
                    phi=Williams(zzone,modes,harms,kappa);
                    phi=phi*tt;
                    L(in,nw+(1:size(phi,2)))=phi;
                    nw=nw+size(phi,2);
                    ntip=ntip+1;
                end
            end
            
        end
    end
    nfree=sum(freenodes);
    indK=2*nfree+(1:size(L,2));
    Ln=sparse(find(freenodes),1:nfree,1,prod(Nnodes),nfree);
    L=[real(L);imag(L)];
    L=[blkdiag(Ln,Ln),L];
    Nddl=size(L,2);
end
%%
LoadMat(nmod);
CreateGradBasisFunction(1,nmod);
K=AssembleStiffnessMatrix(1,nmod);
Kc=[L'*(K*L),L'*C;C'*L,sparse(size(C,2),size(C,2))];
Fc=[L'*Fo;Up];
X=Kc\Fc;
U=L*X(1:Nddl,:);
Fl=K*U;
if extract
    Ks=X(indK,:);
else
    Ks=[];
end
%%
for iim=1:size(U,2)
UpdateInternalState(nmod,U(:,iim),iim);
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
model=param;
%%
for iz=1:size(model.zone,2)
    uu=zeros(2,size(U,2));
    ff=zeros(2,size(U,2));
    zone=model.zone(:,iz);
    switch zone{4}
        case 6
            nodes=zone{6};
            for ii=1:2
                indo=nodes+(ii-1)*prod(Nnodes);
                ff(ii,:)=sum(Fl(indo,:),1);
                uu(ii,:)=mean(U(indo,:),1);
            end
            model.zone{7,iz}=[uu;ff];
    end
end
load(fullfile('TMP','params'),'param');
load(fullfile('TMP','0_mesh_0'),'Nnodes','Nelems','xo','yo','conn','elt','ng','rflag','rint');
try
    tips(~(cell2mat(model.zone(4,:))==5),:)=[];
    model.zone(:,~((cell2mat(model.zone(4,:))==5)|(cell2mat(model.zone(4,:))==6)))=[];
    cracks=find(cell2mat(model.zone(4,:))==5);
catch
end
save(filres,'U','Ks','tips','cracks','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
if fromdic
    ReferenceImage(nmod);
    param.iter_max=2;
    param.analysis='correlation';
    LoadParameters(param)
    [U1]=Solve2D(U,1,nmod);
    copyfile(fullfile('TMP','0_error_0.mat'),strrep(filres,'.res','-error.res'));
    save(filres,'param','-append')
    postproVTK([filres,''],0);
else
    postproVTK([filres,''],0,0);
end

end