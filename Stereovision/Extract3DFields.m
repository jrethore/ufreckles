function [U,X3,Y3,Z3]=Extract3DFields(Uo,nmod,iscale)
if nargin<3,iscale=1;end
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
nurbs=strcmp(param.basis,'nurbs');
if iscell(param0.reference_image)
    ncamr=length(param0.reference_image);
    docrosscorrelation=0;
    if isfield(param0,'cross_correlation')
        docrosscorrelation=param0.cross_correlation;
    end
else
    ncamr=1;
    docrosscorrelation=0;
end
if iscell(param0.deformed_image)
    nim=size(param0.deformed_image,2);
    ncamd=size(param0.deformed_image,1);
else
    nim=1;
    ncamd=1;
end
if ncamr==1
    indcam=ncamd:-1:1;
else
    if docrosscorrelation
        indo=(ncamd:-1:1)';
        indcam=indo';
        for ic=1:ncamd-1
            indo=circshift(indo,1);
            indcam=[indcam;indo'];
        end
    else
        indcam=(1:ncamd)';
    end
end
pscale=2^(iscale-1);
Nddl=size(Uo,1)/2;
dotopo=~(isfield(param,'topography'));%||isfield(param,'mesh_file'));
if ~dotopo
    if nurbs
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(1-1),iscale-1)),'PX','PY','PZ');
        X3=PX(:);Y3=PY(:);Z3=PZ(:);dec=0;
    else
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(1-1),iscale-1)),'Xo','Yo','Zo');
        X3=Xo;Y3=Yo;Z3=Zo;dec=0;
    end
else
    dec=1;
end
U=zeros(3*Nddl,nim-dec);
for ijm=1:nim
    
    phi=[];
    X=[];
    for icamr=1:ncamr
        if ijm==1||ncamr>1
            %  load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icamr-1),10*(1-1))),'phix','Xi','Yi');
            %  load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icamr-1),10*(1-1))),'phiy');
            if nurbs
                load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'Px','Py');
                Xi=Px(:);Yi=Py(:);
            else
                load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'xo','yo');
                Xi=(xo-0.5)*pscale+0.5;
                Yi=(yo-0.5)*pscale+0.5;
            end
            phix=[diag(sparse(ones(length(Xi),1))),sparse(length(Xi),length(Xi))];
            phiy=[sparse(length(Xi),length(Xi)),diag(sparse(ones(length(Xi),1)))];
            phio=phix(:,1:Nddl);
            load(fullfile('TMP',sprintf('sample%d',(icamr-1))),'sizeim','roi');
        end
        for icamd=1:size(indcam,2)
            iz=indcam(icamr,icamd)+(icamr-1)*size(indcam,2)*docrosscorrelation;
            Ux=phix*Uo(:,ijm,iz);
            Uy=phiy*Uo(:,ijm,iz);
            xsi=roi(1)-1+Ux+Xi-0.5*(sizeim(1)+1);
            ysi=roi(3)-1+Uy+Yi-0.5*(sizeim(2)+1);
            C=param0.calibration_data{indcam(icamr,icamd)};
            fx=C(1);fy=C(2);
            cx=C(3);cy=C(4);
            tx=C(5);ty=C(6);tz=C(7);
            %correct distortion
            if length(C)>16
                k1=C(17);k2=C(18);k3=C(19);
                xsio=xsi;
                ysio=ysi;
                res=1;
                while res>1.e-6
                    zs=(xsi-cx+1i*(ysi-cy));
                    rs=abs(zs);
                    d=rs.*(k1+rs.*(k2+k3*rs)).*exp(1i*angle(zs));
                    dxsi=xsio-real(d)-xsi;
                    dysi=ysio-imag(d)-ysi;
                    xsi=xsi+dxsi;
                    ysi=ysi+dysi;
                    res=max(max(abs(dxsi./xsi)),max(abs(dysi./ysi)));
                end
            end
            R=reshape(C(8:16),3,3);
            Xxi=-(xsi-cx)*tz-fx*tx;
            Xyi=-(ysi-cy)*tz-fy*ty;
            %Xxi=(xsi-cx)*tz-fx*tx;
            %Xyi=(ysi-cy)*tz-fy*ty;
            X=[X;Xxi;Xyi];
            phixi=[diag(sparse(fx*R(1,1)+R(3,1)*(xsi-cx)))*phio,diag(sparse(fx*R(1,2)+R(3,2)*(xsi-cx)))*phio,diag(sparse(fx*R(1,3)+R(3,3)*(xsi-cx)))*phio];
            phiyi=[diag(sparse(fy*R(2,1)+R(3,1)*(ysi-cy)))*phio,diag(sparse(fy*R(2,2)+R(3,2)*(ysi-cy)))*phio,diag(sparse(fy*R(2,3)+R(3,3)*(ysi-cy)))*phio];
            %  phixi=[diag(sparse(fx*R(1,1)-R(3,1)*(xsi-cx)))*phio,diag(sparse(fx*R(1,2)-R(3,2)*(xsi-cx)))*phio,diag(sparse(fx*R(1,3)-R(3,3)*(xsi-cx)))*phio];
            %  phiyi=[diag(sparse(fy*R(2,1)-R(3,1)*(ysi-cy)))*phio,diag(sparse(fy*R(2,2)-R(3,2)*(ysi-cy)))*phio,diag(sparse(fy*R(2,3)-R(3,3)*(ysi-cy)))*phio];
            phi=[phi;phixi;phiyi];
        end
    end
    M=phi'*phi;
    F=phi'*X;
    XYZo=M\F;
    Xo=XYZo((1:Nddl));
    Yo=XYZo(Nddl+(1:Nddl));
    Zo=XYZo(2*Nddl+(1:Nddl));
    if ijm==1&&dotopo
        X3=Xo;Y3=Yo;Z3=Zo;
        if ijm==1
            if nurbs
                PX=Xo;PY=Yo;PZ=Zo;
                save(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'PX','PY','PZ','-append');
            else
                save(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(icamr-1),iscale-1)),'Xo','Yo','Zo','-append');
            end
        end
        
    else
        U((1:Nddl),ijm-dec)=Xo-X3;
        U(Nddl+(1:Nddl),ijm-dec)=Yo-Y3;
        U(2*Nddl+(1:Nddl),ijm-dec)=Zo-Z3;
    end
end
end