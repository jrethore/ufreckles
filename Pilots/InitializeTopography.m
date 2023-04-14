function [Uini]=InitializeTopography(Uo,iscale,nmod)
load(fullfile('TMP','params'),'param');
param0=param;
roi=param0.roi;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');

 if strcmp(param.basis,'nurbs')&&iscale==1
 load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nbsnodes');
    Nnodes=Nbsnodes;
 else
load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nnodes');
 end
Nddl=prod(Nnodes);

if size(Uo,1)>prod(Nnodes)
    load(fullfile('TMP','sample0'),'sizeim');
    nurbs=strcmp(param.basis,'nurbs');
    pscale=2^(iscale-1);
    if nurbs&&iscale==1
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Px','Py');
        Xi=Px(:);Yi=Py(:);
    else
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','Nnodes');
        Xi=(xo-0.5)*pscale+0.5;
        Yi=(yo-0.5)*pscale+0.5;
    end
    phi=[];
    X=[];
    nim=size(Uo,3);
    for iim=1:nim
        Ux=Uo((1:Nddl),1,iim);
        Uy=Uo(Nddl+(1:Nddl),1,iim);
        xsi=roi(1)-1+Ux+Xi-0.5*(sizeim(1)+1);
        ysi=roi(3)-1+Uy+Yi-0.5*(sizeim(2)+1);
        C=param0.calibration_data{iim};
        fx=C(1);fy=C(2);
        cx=C(3);cy=C(4);
        tx=C(5);ty=C(6);tz=C(7);
        if length(C)>16
            k1=C(17);k2=C(18);k3=C(19);
            xsio=xsi;
            ysio=ysi;
            res=1;
            while res>1.e-6
                zs=(xsi-cx+i*(ysi-cy));
                rs=abs(zs);
                d=rs.*(k1+rs.*(k2+k3*rs)).*exp(i*angle(zs));
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
        % Xxi=(xsi-cx)*tz-fx*tx;
        % Xyi=(ysi-cy)*tz-fy*ty;
        X=[X;Xxi;Xyi];
        xsi=sparse(xsi);
        ysi=sparse(ysi);
        phixi=[diag(fx*R(1,1)+R(3,1)*(xsi-cx)),diag(fx*R(1,2)+R(3,2)*(xsi-cx)),diag(fx*R(1,3)+R(3,3)*(xsi-cx))];
        phiyi=[diag(fy*R(2,1)+R(3,1)*(ysi-cy)),diag(fy*R(2,2)+R(3,2)*(ysi-cy)),diag(fy*R(2,3)+R(3,3)*(ysi-cy))];
        %   phixi=[diag(fx*R(1,1)-R(3,1)*(xsi-cx)),diag(fx*R(1,2)-R(3,2)*(xsi-cx)),diag(fx*R(1,3)-R(3,3)*(xsi-cx))];
        %   phiyi=[diag(fy*R(2,1)-R(3,1)*(ysi-cy)),diag(fy*R(2,2)-R(3,2)*(ysi-cy)),diag(fy*R(2,3)-R(3,3)*(ysi-cy))];
        phi=[phi;phixi;phiyi];
    end
    M=phi'*phi;
    F=phi'*X;
    XYZo=M\F;
    Xo=XYZo((1:Nddl));
    Yo=XYZo(Nddl+(1:Nddl));
    Zo=XYZo(2*Nddl+(1:Nddl));
    C=param0.calibration_data{1};
    tz=C(7);
    R=reshape(C(8:16),3,3);
    XYZo=[Xo,Yo,Zo];
    Z1=(XYZo*R(3,:)')+tz;
    Uini=Z1;
    if nurbs&&iscale==1
        PX=Xo;PY=Yo;PZ=Zo;
        Xo=CPToNodes(PX,nmod,iscale);
        Yo=CPToNodes(PY,nmod,iscale);
        Zo=CPToNodes(PZ,nmod,iscale);
        save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'PX','PY','PZ','-append');
    end
    save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Xo','Yo','Zo','-append');
    %  figure
    %  plot3(Xo(1),Yo(1),Zo(1),'rx')
    %  hold on
    %  plot3(Xo,Yo,Zo,'o')
    %  xlabel('x')
    %  ylabel('y')
    %  axis equal
else
    
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale)),'unmasked_nodes');
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale)),'xo','yo','Xo','Yo','Zo','Nnodes');
    xg=2*(xo-1)+1;
    yg=2*(yo-1)+1;
    
    xg=reshape(xg,Nnodes);
    yg=reshape(yg,Nnodes);
    Xg=reshape(Xo,Nnodes);
    Yg=reshape(Yo,Nnodes);
    Zg=reshape(Zo,Nnodes);
    if isempty(unmasked_nodes)
        Ug=reshape(Uo,Nnodes);
    else
        Ug=repmat(0,Nnodes);
        Ug(unmasked_nodes)= Uo;
    end
    
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','Nnodes');
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
    
    
    
    
    xo=min(xo,max(xg(:)));xo=max(xo,min(xg(:)));
    yo=min(yo,max(yg(:)));yo=max(yo,min(yg(:)));
    
    
    Uini=interp2(yg,xg,Ug,yo(:),xo(:),'linear');
    Xn=interp2(yg,xg,Xg,yo(:),xo(:),'linear');
    Yn=interp2(yg,xg,Yg,yo(:),xo(:),'linear');
    Zn=interp2(yg,xg,Zg,yo(:),xo(:),'linear');
    
    if ~isempty(unmasked_nodes)
        Uini=Uini(unmasked_nodes);
    end
    if strcmp(param.basis,'nurbs')&&iscale==1
        load(fullfile('TMP',sprintf('%d_phio_%d',nmod,iscale-1)),'phio');
        phi=phio;
        M=phi'*phi;
        Fu=phi'*Uini;
        Uini=M\Fu;
    end
    %    CompleteMeshes(Zn,iscale,nmod,Xn,Yn);
    CompleteMeshes(Uini,iscale,nmod);
end
end