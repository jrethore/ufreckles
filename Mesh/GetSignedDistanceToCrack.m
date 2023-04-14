function [lvl,lvlf,nx,ny]=GetSignedDistanceToCrack(zone,xy,tips)
if nargin<3,tips=1;end
npx=length(xy);
if npx==1
    xy=repmat(xy,2,1);
end
lvlf=0*xy;
lvl=0*xy;
if nargout>2
    nx=lvl;
    ny=lvl;
end
xypo=zone*[1;1i];
seg=diff(xypo(1+(0:1)));
xyp=[xypo(1)-100*max(abs(diff(xypo)))*seg;xypo];
seg=diff(xypo(end-1:end));
xyp=[xyp;xypo(end)+100*max(abs(diff(xypo)))*seg;];

seg=diff(xyp);
n=(-1i*seg)./abs(seg);
np=0.5*(n(1:end-1)+n(2:end));
np=1i*[n(1);np;n(end)];

frontp=cumsum(abs(diff(xyp)));
frontp=[frontp(1);-frontp+frontp(1)];
frontpp=flipud(cumsum(abs(diff(flipud(xyp)))));
frontpp=[-frontpp+frontpp(end);frontpp(end)];
for ip=1:length(xyp)-1
    seg=diff(xyp(ip+(0:1)));
    n=(-1i*seg)/abs(seg);
    t=(seg)/abs(seg);
    cracki=(real((xy-xyp(ip))'*n))';
    fronti=real((xy-0.5*(xyp(ip+1)+xyp(ip)))'*t)';
    frontii=fronti+0.5*(frontpp(ip)+frontpp(ip+1));
    fronti=-fronti+0.5*(frontp(ip)+frontp(ip+1));
    along=(real((xy-xyp(ip))'*np(ip))>=0)&(real((xy-xyp(ip+1))'*np(ip+1))<=0);
    lvl(along)=cracki(along);
    if tips
    lvlf(along)=max(fronti(along),frontii(along));
    else
    lvlf(along)=fronti(along);
    end
    if nargout>2
        nx(along)=real(n);
        ny(along)=imag(n);
    end
%       figure
%       plot(xypo,'k')
%       hold on
%       scatter(real(xy(along)),imag(xy(along)),[],fronti(along))
%       axis equal
%       figure
%       plot(xypo,'k')
%       hold on
%       scatter(real(xy(along)),imag(xy(along)),[],frontii(along))
%       axis equal
    
end
if nargout>2

    for ip=3:length(xyp)-2
       
        segm=diff(xyp(ip+(-1:0)));
        segp=diff(xyp(ip+(0:1)));
        da=-angle(segp*exp(-1i*angle(segm)));
        if abs(da)>0
            lc=10*min(abs(segp),abs(segm))/(5);
            zc=xyp(ip)-sign(da)*np(ip)*1i*lc;
%             figure
%             plot(xypo,'k')
%        hold on
%        axis equal
%              plot(xyp(ip),'ro')
%               plot(zc,'rx')

            xyloc=(xy-zc)*(sign(da)*np(ip)*1i)';
            
%             figure
%             plot(xyloc,'xk')
%        hold on
%        axis equal
%              plot((xypo-zc)*(np(ip)*1i)','k')
%            
%           keyboard
            
            ac=atan2(lc/3,lc);
            nx(abs(angle(xyloc))<ac)=-sign(da)*real((xy(abs(angle(xyloc))<ac)-zc))./abs(xy(abs(angle(xyloc))<ac)-zc);
            ny(abs(angle(xyloc))<ac)=-sign(da)*imag((xy(abs(angle(xyloc))<ac)-zc))./abs(xy(abs(angle(xyloc))<ac)-zc);
%             figure
%        hold on
%        quiver(real(xyloc(abs(angle(xyloc))<ac)),imag(xyloc(abs(angle(xyloc))<ac)),nx(abs(angle(xyloc))<ac),ny(abs(angle(xyloc))<ac))
%        axis equal
%             figure
%             plot(xypo,'k')
%        hold on
%        quiver(real(xy(abs(angle(xyloc))<ac)),imag(xy(abs(angle(xyloc))<ac)),nx(abs(angle(xyloc))<ac),ny(abs(angle(xyloc))<ac))
%        axis equal
        end
    end
end

if npx==1
    lvl=lvl(1,:);
    lvlf=lvlf(1,:);
    if nargout>2
        nx=nx(1,:);
        ny=ny(1,:);
    end
end
end