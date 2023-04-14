function [phi]=CreateKMBasis(modes,ind,kappa,ipix,ic,dec,z)

if nargin<4, ipix=[];end
if nargin<5, ic=1;end
if nargin<6, dec=0;end
if nargin<7, z=[];end
load(fullfile('TMP','sample0_0'),'sizeim');
sizeim=prod(sizeim);
if ~isempty(z)
    dist=abs(z);
    angl=angle(z);
    clear z
    load(sprintf('TMP/%d_levelsets_cylco',ic),'theta');
    thetap=theta;
else
if dec==0
    load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'dist','angl','theta','thetap');
    if ~isempty(ipix)
        dist=dist(ipix);
        angl=angl(ipix);
            thetap=thetap(ipix);
    end
else
    load(fullfile('TMP',sprintf('%d_levelsets',ic)),'crack','front');
 load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'thetap');
    if ~isempty(ipix)
        front=front(ipix);
        crack=crack(ipix);
            thetap=thetap(ipix);
    end
    z=front-dec+i*crack;
    clear front crack
    dist=abs(z);
    angl=angle(z);
    clear z
    load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'theta');

end
end
% display('WARNING : normalized distance');
% dist=dist/max(dist(:));

if isempty(ipix)

    phi=zeros(length(dist(:)),length(modes)*length(ind));

else
    val=zeros(length(ipix)*length(modes)*length(ind),1);
    indp=zeros(length(ipix)*length(modes)*length(ind),1);
    indm=zeros(length(ipix)*length(modes)*length(ind),1);
end
inds=0;
icon=0;
for m=1:length(modes)
    for ii=1:length(ind)
        icon=icon+1;
        meq=modes(m);
        heq=ind(ii);
        [harm1,amp1]=KMPotFourrier(meq,heq,kappa);
        feq=0*dist(:);
        for kk=1:3
            feq=feq+amp1(kk)*exp(harm1(kk)*angl(:)*i/2);
        end
        feq=feq.*(dist(:).^(heq/2));
%        feq=feq*exp(i*theta);
        feq=feq.*exp(i*thetap(:));
        if isempty(ipix)
            phi(:,icon)=feq;
        else
            val(inds+(1:length(ipix)))=feq;
            indp(inds+(1:length(ipix)))=ipix(:);
            indm(inds+(1:length(ipix)))=icon;
            inds=inds+length(ipix);

        end

    end
end

if ~isempty(ipix)

    phi=sparse(indp,indm,val,sizeim,length(modes)*length(ind));

end


end
