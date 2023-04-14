function [tn,tt]=ComputeCohesiveStress(un,ut,nmod)

load(fullfile('TMP','params'),'param');
pix2m=param.pixel_size;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
matmod=param.interface_model;

switch matmod
    case 'cohesive'
        load(fullfile('TMP',sprintf('%d_intmod',nmod)),'model','ueqmax');
un=pix2m*un;
ut=pix2m*ut;
ueq=abs(max(0,un)+i*model.alpha*model.beta*ut);
tmax=interp1(model.uo,model.to,max(ueq,ueqmax),'linear','extrap');
tcon=interp1(model.uo,model.to,un,'linear','extrap');
tn=tcon.*(un<=0)+tmax.*((ueq>=ueqmax)&(un>0)).*un./ueq+tmax.*un./abs(ueqmax).*((ueq<ueqmax)&(un>0));
tt=(model.beta)*(tmax.*(ueq>=ueqmax).*ut./ueq+tmax.*ut./abs(ueqmax).*(ueq<ueqmax));
if any(isnan(tn))
    found=find(isnan(tn));
    tn(found)=0;
end
if any(isnan(tt))
    found=find(isnan(tt));
    tt(found)=0;
end
% figure
% plot(ueq,'ro')
% hold on
% plot(max(0,ueqmax),'bx')
% title('Cohesive stress')
% figure
% plot(model.uo,model.to,'k')
% hold on
% plot(un,tn,'rx')
% plot(ut,tt,'bx')
% plot(max(ueq,ueqmax),tmax,'sg')
% legend({'model','n','t','max'})
% title('Cohesive stress')
% pause
    case 'plastic'
        load(fullfile('TMP',sprintf('%d_intmod',nmod)),'model','ueqmax');
un=pix2m*un;
ut=pix2m*ut;
ueq=abs(un+i*model.alpha*model.beta*ut);
ur=max(0,ueqmax-model.up);
kn=model.stiffness;
tn=0*(un-ur)*kn;
tt=-(ut-ur)*kn;

    otherwise
        error ('INVALID INTERFACE MODEL');

end




end