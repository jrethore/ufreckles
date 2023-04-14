function [dphidx,dphidy]=CreateGradKMBasis(modes,ind,kappa,ipix,ic,dec,z)

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
        dist=max(abs(z),1);
        angl=angle(z);
        clear z
        load(fullfile('TMP',sprintf('%d_levelsets_cylco',ic)),'theta');

    end
end
% display('WARNING : normalized distance');
% rmax=1/max(dist(:));
% dist=dist*rmax;
if isempty(ipix)
    dphidx=zeros(length(dist(:)),length(modes)*length(ind));
    dphidy=zeros(length(dist(:)),length(modes)*length(ind));

else

    valx=zeros(length(ipix)*length(modes)*length(ind),1);
    valy=zeros(length(ipix)*length(modes)*length(ind),1);
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

        ft=0*dist(:);
        fr=0*dist(:);

        for kk=1:3
            fr=fr+amp1(kk)*exp(harm1(kk)*angl(:)*i/2);
            ft=ft+harm1(kk)*i/2*amp1(kk)*exp(harm1(kk)*angl(:)*i/2);
        end
        ft=ft.*(dist(:).^(heq/2-1));
        fr=heq/2*fr.*(dist(:).^(heq/2-1));

        dfndth=imag(ft);
        dfndr=imag(fr);
        dftdth=real(ft);
        dftdr=real(fr);

        drdt=cos(angl(:));
        drdn=sin(angl(:));
        dthdt=-sin(angl(:));
        dthdn=cos(angl(:));

        dfndt=dfndr.*drdt+dfndth.*dthdt;
        dfndn=dfndr.*drdn+dfndth.*dthdn;
        dftdt=dftdr.*drdt+dftdth.*dthdt;
        dftdn=dftdr.*drdn+dftdth.*dthdn;

        dftdx=dftdt.*cos(thetap(:))-dftdn.*sin(thetap(:));
        dftdy=dftdt.*sin(thetap(:))+dftdn.*cos(thetap(:));
        dfndx=dfndt.*cos(thetap(:))-dfndn.*sin(thetap(:));
        dfndy=dfndt.*sin(thetap(:))+dfndn.*cos(thetap(:));

        dfxdx=(dftdx.*cos(thetap(:))-dfndx.*sin(thetap(:)));
        dfxdy=(dftdy.*cos(thetap(:))-dfndy.*sin(thetap(:)));
        dfydx=(dftdx.*sin(thetap(:))+dfndx.*cos(thetap(:)));
        dfydy=(dftdy.*sin(thetap(:))+dfndy.*cos(thetap(:)));

        %         dftdx=dftdt*cos(theta)-dftdn*sin(theta);
        %         dftdy=dftdt*sin(theta)+dftdn*cos(theta);
        %         dfndx=dfndt*cos(theta)-dfndn*sin(theta);
        %         dfndy=dfndt*sin(theta)+dfndn*cos(theta);
        %
        %         dfxdx=dftdx*cos(theta)-dfndx*sin(theta);
        %         dfxdy=dftdy*cos(theta)-dfndy*sin(theta);
        %         dfydx=dftdx*sin(theta)+dfndx*cos(theta);
        %         dfydy=dftdy*sin(theta)+dfndy*cos(theta);


        %         if m==1
        %                 feq=0*dist(:);
        %
        %         for kk=1:3
        %                        feq=feq+amp1(kk)*exp(harm1(kk)*angl(:)*i/2);
        %
        %         end
        %         feq=feq.*(dist(:).^(heq/2));
        %         [gdfndy,gdfndx]=gradient(reshape(imag(feq),200,200));
        %         [gdftdy,gdftdx]=gradient(reshape(real(feq),200,200));
        %
        %         figure
        %         imagesc(gdftdx)
        %         colorbar;title('dftdx grad')
        %         figure
        %         imagesc(reshape(dftdx,200,200))
        %         colorbar;title('dftdx')
        %         end


        %                 if m==1
        %        feq=feq*exp(i*theta);
        %         [gdfydy,gdfydx]=gradient(reshape(imag(feq),200,200));
        %         [gdfxdy,gdfxdx]=gradient(reshape(real(feq),200,200));
        %
        %         figure
        %         imagesc(gdfydx)
        %         colorbar;title('dfydx grad')
        %         figure
        %         imagesc(reshape(dfydx,200,200))
        %         colorbar;title('dfydx')
        %         end




        if isempty(ipix)
            dphidx(:,icon)=dfxdx+i*dfydx;
            dphidy(:,icon)=dfxdy+i*dfydy;
        else
            valx(inds+(1:length(ipix)))=dfxdx+i*dfydx;
            valy(inds+(1:length(ipix)))=dfxdy+i*dfydy;
            indp(inds+(1:length(ipix)))=ipix(:);
            indm(inds+(1:length(ipix)))=icon;
            inds=inds+length(ipix);

        end

    end
end

if ~isempty(ipix)

    dphidx=sparse(indp,indm,valx,sizeim,length(modes)*length(ind));
    dphidy=sparse(indp,indm,valy,sizeim,length(modes)*length(ind));

end


end
