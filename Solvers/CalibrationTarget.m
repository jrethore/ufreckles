function [T,residus]=CalibrationTarget(nmod)
clear GetTargetPosition
load(fullfile('TMP','params'),'param');
param0=param;
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
ngrid=param.grid_size;
lgrid=param.grid_step;
if iscell(param0.deformed_image)
    ncam=size(param0.deformed_image,1);
    nim=size(param0.deformed_image,2);
    if nim<2
        error('YOU MAY GIVE AT LEAST 2 DEFORMED IMAGES FOR CALIBRATION');
    end
else
    error('YOU MAY GIVE AT LEAST 2 DEFORMED IMAGES FOR CALIBRATION');
end
disp(sprintf('Starting calibration...'));




[Yw,Xw]=meshgrid(1:ngrid(2),1:ngrid(1));
Xw=(Xw(:)-mean(Xw(:)))*lgrid;
Yw=(Yw(:)-mean(Yw(:)))*lgrid;
phi0=sparse(prod(ngrid),1);
phi1=sparse(ones(prod(ngrid),1));
phixo=[phi1,phi0,Xw,Yw,phi0,phi0,phi0,phi0];
phiyo=[phi0,phi1,phi0,phi0,Xw,Yw,phi0,phi0];
phis=[phi0,phi0,phi0,phi0,phi0,phi0,Xw,Yw];

oks=ones(nim,ncam);
Ut=zeros(8,nim,ncam);
Xt=zeros(numel(Xw),nim,ncam);
Yt=zeros(numel(Yw),nim,ncam);
scrsz = get(0,'ScreenSize');
fig=figure(100);
set(fig,'NumberTitle','off','Name','Calibration spot detection','Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2],'Renderer','painters');
ax=axes('Parent',fig,'drawmode','fast','FontSize',16);
gim=[];U=[];V=[];pscale=1;
for icam=1:ncam
    for ijm=1:nim
        if icam==1||~any(~oks(ijm,:))
            fildef=param0.deformed_image{icam,ijm};
            im1=double(readim(fildef));
            if length(size(im1))==3
                im1=mean(im1,3);
            end
            if reverse
                im1=im1';
            end
            %            display(sprintf('Camera %d Image %d : %s',icam,ijm,fildef));
            if isempty(gim)
                gim=imagesc(im1');
                colormap('gray')
                hold on
                axis equal
                axis xy;
                axis off
            end
            set(fig,'Name',sprintf('Calibration spot detection for Camera %d Image %d: %s',icam,ijm,fildef));
            % if ~(icam==1&&ijm==1)
            %     [U,V]=rbt(im0,MCoarseImage(im1(roi(1):roi(2),roi(3):roi(4)),nscale))
            % end
            [Xs,Ys,ok]=  GetTargetPosition(nmod,im1,U*pscale,V*pscale,gim);
            %            if icam==1&&ijm==1
            %                roi=round([min(size(im1,1),max(1,min(Xs)+0.5*(size(im1,1)-1)-5)),...
            %                    min(size(im1,1),max(1,max(Xs)+0.5*(size(im1,1)-1)+5)),...
            %                min(size(im1,2),max(1,min(Ys)+0.5*(size(im1,2)-1)-5)),...
            %                    min(size(im1,2),max(1,max(Ys)+0.5*(size(im1,1)-2)+5))])
            %                size(im1)
            % nscale=max(3,ceil(log(max([roi(2)-roi(1),roi(4)-roi(3)])/512)/log(2)));
            % nscale=1;
            % pscale=2^(nscale-1);
            %                im0=MCoarseImage(im1(roi(1):roi(2),roi(3):roi(4)),nscale);
            %            end
            
            phix=phixo-diag(sparse(Xs))*phis;
            phiy=phiyo-diag(sparse(Ys))*phis;
            M=[phix;phiy]'*[phix;phiy];
            F=[phix;phiy]'*[Xs;Ys];
            Ui=M\F;
            %     res=norm(F-M*Ui)/norm(F);
            %     display(sprintf('Normalized positioning error for camera %d image %d : %f %%',icam,ijm,res*100));
            
            if any(isnan(Ui))||any(isnan(Xs))||any(isnan(Ys))
                ok=0;
            end
            %mean(abs([phix;phiy]*Ui-[Xs;Ys]))
            % figure
            % plot(Xs,Ys,'bx')
            % hold on
            % plot(phix*Ui,phiy*Ui,'r+')
            % pause
            oks(ijm,icam)=ok;
            Ut(:,ijm,icam)=Ui;
            Xt(:,ijm,icam)=Xs;
            Yt(:,ijm,icam)=Ys;
        else
            oks(ijm,icam)=0;
        end
    end
end
delete(fig);
%%
% for ijm=1:size(Ut,2)
%     if any(~oks(ijm,:))
% Ut(:,ijm,:)=[];
% Xt(:,ijm,:)=[];
% Yt(:,ijm,:)=[];
%     end
% end
T=zeros(16,ncam);
residus=zeros(ncam,1);
for icam=1:ncam
    U=Ut(:,:,icam);
    U=U(:,oks(:,icam)>0);
    toremove=sum(isnan(U),1)>0;
    U(:,toremove)=[];
    A=[U(3,:).*U(8,:)+U(4,:).*U(7,:),2*(U(3,:).*U(7,:)-U(4,:).*U(8,:));...
        U(5,:).*U(6,:)               ,   U(5,:).*U(5,:)-U(6,:).*U(6,:) ;...
        U(5,:).*U(8,:)+U(6,:).*U(7,:),2*(U(5,:).*U(7,:)-U(6,:).*U(8,:));...
        U(7,:).*U(8,:)               ,   U(7,:).*U(7,:)-U(8,:).*U(8,:) ]';
    b=[-U(3,:).*U(4,:),U(4,:).*U(4,:)-U(3,:).*U(3,:)]';
    xsi=A\b;
    res=norm(b-A*xsi)/norm(b);
    residus(icam)=res;
    display(sprintf('Normalized residual for camera %d : %f %%',icam,res*100));
    cx=-xsi(1);
    cy=-xsi(3)/xsi(2);
    fx=(sqrt(xsi(4)-(xsi(3))^2/xsi(2)-(xsi(1))^2));
    fy=fx/(sqrt(xsi(2)));
    tz=-2*(sqrt(1/(((U(3,1)-cx*U(7,1))^2+(U(4,1)-cx*U(8,1))^2)/fx^2+((U(5,1)-cy*U(7,1))^2+(U(6,1)-cy*U(8,1))^2)/fy^2+U(7,1)^2+U(8,1)^2)));
    tx=-(U(1,1)-cx)*tz/fx;
    ty=-(U(2,1)-cy)*tz/fy;
    R=zeros(3,3);
    R(1,1)=-tz*(U(3,1)-cx*U(7,1))/fx;
    R(1,2)=-tz*(U(4,1)-cx*U(8,1))/fx;
    R(2,1)=-tz*(U(5,1)-cy*U(7,1))/fy;
    R(2,2)=-tz*(U(6,1)-cy*U(8,1))/fy;
    R(3,1)=tz*U(7,1);
    R(3,2)=tz*U(8,1);
    for ir=1:3
        for jr=1:3
            for kr=1:3
                eijk=0.5*(ir-jr)*(jr-kr)*(kr-ir);
                R(ir,3)=R(ir,3)+eijk*R(jr,1)*R(kr,2);
            end
        end
    end
    [UR,SR,VR]=svd(R);
    Rn=UR*VR';
    
    T(:,icam)=[fx;fy;cx;cy;tx;ty;tz;Rn(:)];
end
%%

C10=T(5:7,1);
R10=reshape(T(8:16,1),3,3);
R01=inv(R10);
OC1=-R01*C10;

C20=T(5:7,2);
R20=reshape(T(8:16,2),3,3);
R02=inv(R20);
OC2=-R02*C20;


x1=R01(:,1);
y1=R01(:,2);
z1=R01(:,3);
xn=x1;
yn=y1;
zn=z1;
ON=OC1;

% C1C2=-OC1+OC2;
% xn=C1C2/norm(C1C2);
% ON=OC1+0.5*C1C2;
% x1=R01(:,1);
% x2=R02(:,1);
% a=(-x1'*xn)/(x2'*xn);
% tmp=-(x1+a*x2)/norm(x1+a*x2);
% if find(tmp==max(abs(tmp)))==2
% yn=-tmp/max(tmp);
% zn=cross(xn,yn);
% else
%     zn=tmp;
% yn=cross(zn,xn);
% end
RON=[xn,yn,zn];
R1N=R10*RON;
R2N=R20*RON;
C1N=C10+R10*ON;
C2N=C20+R20*ON;
%%

T(8:16,1)=R1N(:);
T(5:7,1)=C1N;
T(8:16,2)=R2N(:);
T(5:7,2)=C2N;

for icam=3:ncam
    Ci0=T(5:7,icam);
    Ri0=reshape(T(8:16,icam),3,3);
    RiN=Ri0*RON;
    CiN=Ci0+Ri0*ON;
    T(8:16,icam)=RiN(:);
    T(5:7,icam)=CiN;
end
%% correction
icamr=1;
phio=diag(sparse(ones(prod(ngrid),1)));
d=0;
for ijm=1:nim
    phi=[];
    X=[];
    for icam=1:ncam
            xsi=Xt(:,ijm,icam);
            ysi=Yt(:,ijm,icam);
            C=T(:,icam);
            fx=C(1);fy=C(2);
            cx=C(3);cy=C(4);
            tx=C(5);ty=C(6);tz=C(7);
            R=reshape(C(8:16),3,3);
            Xxi=-(xsi-cx)*tz-fx*tx;
            Xyi=-(ysi-cy)*tz-fy*ty;
            X=[X;Xxi;Xyi];
            phixi=[diag(sparse(fx*R(1,1)+R(3,1)*(xsi-cx)))*phio,diag(sparse(fx*R(1,2)+R(3,2)*(xsi-cx)))*phio,diag(sparse(fx*R(1,3)+R(3,3)*(xsi-cx)))*phio];
            phiyi=[diag(sparse(fy*R(2,1)+R(3,1)*(ysi-cy)))*phio,diag(sparse(fy*R(2,2)+R(3,2)*(ysi-cy)))*phio,diag(sparse(fy*R(2,3)+R(3,3)*(ysi-cy)))*phio];
            phi=[phi;phixi;phiyi];        
    end
    M=phi'*phi;
    F=phi'*X;
    XYZo=M\F;
    [dx2,dx1]=gradient(reshape(XYZo((1:prod(ngrid))),ngrid));
    [dy2,dy1]=gradient(reshape(XYZo((1:prod(ngrid))),ngrid));
    [dz2,dz1]=gradient(reshape(XYZo((1:prod(ngrid))),ngrid));
    d=d+sum(sum(sqrt(dx1.^2+dy1.^2+dz1.^2))+sum(sqrt(dx2.^2+dy2.^2+dz2.^2)));
end
d=d/(2*prod(ngrid)*nim);
%%
f_corr=lgrid/d;
f_corr=1;

C10=T(5:7,1)*f_corr;
R10=reshape(T(8:16,1),3,3)*f_corr;
R01=inv(R10);
OC1=-R01*C10;

C20=T(5:7,2)*f_corr;
R20=reshape(T(8:16,2),3,3)*f_corr;
R02=inv(R20);
OC2=-R02*C20;


x1=R01(:,1);
y1=R01(:,2);
z1=R01(:,3);
xn=x1;
yn=y1;
zn=z1;
ON=OC1;

RON=[xn,yn,zn];
R1N=R10*RON;
R2N=R20*RON;
C1N=C10+R10*ON;
C2N=C20+R20*ON;
%%

T(8:16,1)=R1N(:);
T(5:7,1)=C1N;
T(8:16,2)=R2N(:);
T(5:7,2)=C2N;

for icam=3:ncam
    Ci0=T(5:7,icam)*f_corr;
    Ri0=reshape(T(8:16,icam),3,3)*f_corr;
    RiN=Ri0*RON;
    CiN=Ci0+Ri0*ON;
    T(8:16,icam)=RiN(:);
    T(5:7,icam)=CiN;
end

end
