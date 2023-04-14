function [tn,tt]=ComputeTangentCohesiveStress(un,ut,nmod)

load(fullfile('TMP','params'),'param');
pix2m=param.pixel_size;

load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
matmod=param.interface_model;

switch matmod
    case 'cohesive'
        load(fullfile('TMP',sprintf('%d_intmod',nmod)),'model','ueq','ueqmax');
        uo=model.uo;
        to=model.to;
        duo=[diff(uo(:)'),max(diff(uo(:)'))];
        unn=[uo(:)';uo(:)'+0.999999*duo];
        dun=repmat(diff(uo(:)'),2,1);
        dtn=repmat(diff(to(:)'),2,1);
        dtndun=[dtn(:)./dun(:);0;0];
        dtdui=interp1(unn(:),dtndun,max(0,ueqmax),'linear','extrap');
        tin=interp1(model.uo,model.to,sign(ueq).*max(0,ueqmax),'linear','extrap');
        tit=interp1(model.uo,model.to,max(0,ueqmax),'linear','extrap');
% kn=dtdui.*(ueq>=ueqmax)+tin./abs(ueqmax).*(ueq<ueqmax);
% kt=dtdui.*(ueq>=ueqmax)+tit./abs(ueqmax).*(ueq<ueqmax);
kn=dtdui.*(ueq>=ueqmax)+tin./abs(ueqmax).*(ueq<ueqmax);
kt=dtdui.*(ueq>=ueqmax)+tit./abs(ueqmax).*(ueq<ueqmax);
        tn=diag(sparse(kn*pix2m))*un;
        tt=diag(sparse(kt*pix2m*(model.beta)))*ut;
%          figure
%          plot(ueq,'bo')
%          hold on
%          plot(max(0,ueqmax),'rx')
% title('tangent Cohesive stress')
% %         figure
% %         plot(dtdui)
%         figure
%         plot(unn(:),dtndun,'k-o')
%         hold on
%         plot(ueq,tin,'r')
%         plot(ueq,(kn),'ks')
%         plot(ueq,(kt),'bs')
%  title('tangent Cohesive stress')
%        pause
    case 'plastic'
         load(fullfile('TMP',sprintf('%d_intmod',nmod)),'model','ueqmax');
      ko=model.stiffness;
 kt=ko*ones(size(ueqmax));
 kn=00000*kt;
         tn=diag(sparse(kn*pix2m))*un;
        tt=diag(sparse(kt*pix2m*(model.beta)))*ut;
    otherwise
        error ('INVALID INTERFACE MODEL');

end




end