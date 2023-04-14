function dg=GetEffectiveGain(anglc)
load('cell-mesh','n','zp','zpp','indc','homo')
zpp=zpp(:);
angl=angle(zpp(homo)-repmat(zpp(1:size(homo,1)),1,size(homo,2)));

% dg=zeros(size(anglc));
% for ip=1:numel(angl(:))
%    dg=max(dg,normpdf(anglc,angl(ip),2*pi/(4*numel(angl(:))))); 
% end
% dg=1-0.8*dg*sqrt(2*pi)*2*pi/(4*numel(angl(:)));
dg=zeros(size(anglc));
for ip=1:numel(angl(:))
   dg=max(dg,abs(anglc-angl(ip))<10*pi/180);
end
dg=1-0.99999999*dg;
end