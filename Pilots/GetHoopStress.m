function [soo]=GetHoopStress(X,nmod,ic)

load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
mu=param.mu;
lambda=param.lambda;

load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'angl');
nx=-sin(angl);
ny=cos(angl);
load(fullfile('TMP',sprintf('%d_epsxx_0',nmod)),'epsxx','sizeim');
exx=reshape((epsxx*X),sizeim);
load(fullfile('TMP',sprintf('%d_epsyy_0',nmod)),'epsyy');
eyy=reshape((epsyy*X),sizeim);
load(fullfile('TMP',sprintf('%d_epsxy_0',nmod)),'epsxy');
exy=reshape((epsxy*X),sizeim);

sxx=(lambda+2*mu)*exx+lambda*eyy;
syy=lambda*exx+(lambda+2*mu)*eyy;
sxy=mu*exy;
%szz=lambda*exx+lambda*eyy;
soo=nx.*(sxx.*nx+sxy.*ny)+ny.*(sxy.*nx+syy.*nx);
%soo=ny.*(sxx.*nx+sxy.*ny)-nx.*(sxy.*nx+syy.*nx);
%soo=sqrt(0.5*((sxx-syy).^2+(syy-szz).^2+(szz-sxx).^2+6*sxy.^2));
%soo=0.5*(sxx+syy)+sqrt((0.5*(sxx-syy)).^2+sxy.^2);

% figure
% imagesc(syy)
% colorbar
% figure
% imagesc(sxy)
% colorbar
% figure
% imagesc(sxx)
% colorbar
end