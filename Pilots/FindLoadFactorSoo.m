function [lf,Smax,Oc]=FindLoadFactorSoo(Soos,Soot,d,Sc,ic,lfo)
if nargin<5, ic=1;end

load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'dist');
xc=find(dist(:)==min(dist(:)));
[xfo,yfo]=ind2sub(size(Soos),xc);
xfo=mean(xfo);
yfo=mean(yfo);
thet=-pi/2:pi/180:pi/2;

xon=d*cos(thet)+xfo;
yon=d*sin(thet)+yfo;


yt=interp2(Soot,yon,xon);
ys=interp2(Soos,yon,xon);


if nargin<6

Smax=max(yt);
 lf=Sc/Smax;
 dlf=1;
 iter=0;
%   figure
%   hold on
while abs(dlf)>1.e-6&&iter<1000
    Smax=max(ys+lf*yt);
%plot(lf,Smax,'o');
%plot(thet,ys+lf*yt)
dlf=0.025*(1-Smax/Sc)/lf;
lf=lf+dlf*lf;
iter=iter+1;
end
     Smax=max(ys+lf*yt);
Oc=mean(thet(find(ys+lf*yt==Smax)));
% figure
% plot(thet,ys+lf*yt)
if iter==1000
    display('WARNING: FAILURE TO CONVERGE IN 1000 ITERATIONS IN FindLoadFactor.m');
end
else
lf=lfo;
     Smax=max(ys+lf*yt);
Oc=mean(thet(find(ys+lf*yt==Smax)));
end
end