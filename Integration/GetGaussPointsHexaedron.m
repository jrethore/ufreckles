function [xg,yg,zg,wg]=GetGaussPointsHexaedron(nb_gauss_point,nb_sub_cell,xn,yn,zn,xpix,ypix,zpix)
if nargin==1, nb_sub_cell=1;end
check=0;
switch nb_gauss_point
    case 1
        xg=0;
        yg=0;
        zg=0;
        wg=8;

    case 8

        xg=[-1   1  -1   1  -1   1  -1   1]'/sqrt(3);
        yg=[-1  -1   1   1  -1  -1   1   1]'/sqrt(3);
        zg=[-1  -1  -1  -1   1   1   1   1]'/sqrt(3);
        wg = [1.,1.,1.,1.,1.,1.,1.,1.]';
    case 888

        xi = [-.9602898564975362;-.7966664774136267;-.5255324099163290;-.1834346424956498;.1834346424956498;...
            .5255324099163290;.7966664774136267;.9602898564975362];
        wi= [.1012285362903763;.2223810344533745;.3137066458778873;.3626837833783620;...
            .3626837833783620;.3137066458778873;.2223810344533745;.1012285362903763];
        [yg,xg,zg]=meshgrid(xi,xi,xi);
        wg=wi*wi';
        wg=wg(:)*wi';
        wg=wg(:);
        xg=xg(:);
        yg=yg(:);
        zg=zg(:);
    case 151515

        xi = [-.9879925180204854;-.9372733924007059;-.8482065834104272;...
            -.7244177313601700;-.5709721726085388;-.3941513470775634;...
            -.2011940939974345;0.; .2011940939974345;...
            .3941513470775634;.5709721726085388;.7244177313601700;...
            .8482065834104272;.9372733924007059;.9879925180204854];
        wi = [.03075324199611807;.07036604748811134;.1071592204671351;...
            .1395706779261761;.1662692058169852;.1861610000155741;...
            .1984314853271374;.2025782419255562;.1984314853271374;...
            .1861610000155741;.1662692058169852;.1395706779261761;...
            .1071592204671351;.07036604748811134;.03075324199611807];
        [yg,xg,zg]=meshgrid(xi,xi,xi);
        wg=wi*wi';
        wg=wg(:)*wi';
        wg=wg(:);
        xg=xg(:);
        yg=yg(:);
        zg=zg(:);

    case 0
        % xg=0*xpix;
        % yg=0*ypix;
        % zg=0*zpix;
        %         for ig=1:length(xg)
        %             xgi=0;ygi=0;zgi=0;
        %             xi=xpix(ig);yi=ypix(ig);zi=zpix(ig);
        %             res=1;
        %             while res>1.e-6
        %                 N=[0.125*(1-xgi).*(1-ygi).*(1-zgi),0.125*(1+xgi).*(1-ygi).*(1-zgi),0.125*(1+xgi).*(1+ygi).*(1-zgi),0.125*(1-xgi).*(1+ygi).*(1-zgi),...
        %                     0.125*(1-xgi).*(1-ygi).*(1+zgi),0.125*(1+xgi).*(1-ygi).*(1+zgi),0.125*(1+xgi).*(1+ygi).*(1+zgi),0.125*(1-xgi).*(1+ygi).*(1+zgi)];
        %                 N_r=[-0.125*(1-ygi).*(1-zgi),0.125*(1-ygi).*(1-zgi),0.125*(1+ygi).*(1-zgi),-0.125*(1+ygi).*(1-zgi),...
        %                     -0.125*(1-ygi).*(1+zgi),0.125*(1-ygi).*(1+zgi),0.125*(1+ygi).*(1+zgi),-0.125*(1+ygi).*(1+zgi)];
        %                 N_s=[-0.125*(1-xgi).*(1-zgi),-0.125*(1+xgi).*(1-zgi),0.125*(1+xgi).*(1-zgi),0.125*(1-xgi).*(1-zgi),...
        %                     -0.125*(1-xgi).*(1+zgi),-0.125*(1+xgi).*(1+zgi),0.125*(1+xgi).*(1+zgi),0.125*(1-xgi).*(1+zgi)];
        %                 N_t=[-0.125*(1-xgi).*(1-ygi),-0.125*(1+xgi).*(1-ygi),-0.125*(1+xgi).*(1+ygi),-0.125*(1-xgi).*(1+ygi),...
        %                     0.125*(1-xgi).*(1-ygi),0.125*(1+xgi).*(1-ygi),0.125*(1+xgi).*(1+ygi),0.125*(1-xgi).*(1+ygi)];
        %
        %                 dxdr=N_r*xn;
        %                 dydr=N_r*yn;
        %                 dzdr=N_r*zn;
        %                 dxds=N_s*xn;
        %                 dyds=N_s*yn;
        %                 dzds=N_s*zn;
        %                 dxdt=N_t*xn;
        %                 dydt=N_t*yn;
        %                 dzdt=N_t*zn;
        %                 J=[dxdr,dxds,dxdt;...
        %                     dydr,dyds,dydt;...
        %                     dzdr,dzds,dzdt];
        %                 xp=N*xn;
        %                 yp=N*yn;
        %                 zp=N*zn;
        %                 dx=J\[xi-xp;yi-yp;zi-zp];
        %                 res=dx'*dx;
        %                 xgi=xgi+dx(1);
        %                 ygi=ygi+dx(2);
        %                 zgi=zgi+dx(3);
        %             end
        %             xg(ig)=xgi;yg(ig)=ygi;zg(ig)=zgi;
        %         end

npix=numel(xpix);
deleted=0*xpix;
        xg=0*xpix;
        yg=0*ypix;
        zg=0*zpix;
        res=1;
        while res>1.e-6
            N=[0.125*(1-xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1+yg).*(1-zg),0.125*(1-xg).*(1+yg).*(1-zg),...
                0.125*(1-xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1+yg).*(1+zg),0.125*(1-xg).*(1+yg).*(1+zg)];
            N_r=[-0.125*(1-yg).*(1-zg),0.125*(1-yg).*(1-zg),0.125*(1+yg).*(1-zg),-0.125*(1+yg).*(1-zg),...
                -0.125*(1-yg).*(1+zg),0.125*(1-yg).*(1+zg),0.125*(1+yg).*(1+zg),-0.125*(1+yg).*(1+zg)];
            N_s=[-0.125*(1-xg).*(1-zg),-0.125*(1+xg).*(1-zg),0.125*(1+xg).*(1-zg),0.125*(1-xg).*(1-zg),...
                -0.125*(1-xg).*(1+zg),-0.125*(1+xg).*(1+zg),0.125*(1+xg).*(1+zg),0.125*(1-xg).*(1+zg)];
            N_t=[-0.125*(1-xg).*(1-yg),-0.125*(1+xg).*(1-yg),-0.125*(1+xg).*(1+yg),-0.125*(1-xg).*(1+yg),...
                0.125*(1-xg).*(1-yg),0.125*(1+xg).*(1-yg),0.125*(1+xg).*(1+yg),0.125*(1-xg).*(1+yg)];

            dxdr=N_r*xn;
            dydr=N_r*yn;
            dzdr=N_r*zn;
            dxds=N_s*xn;
            dyds=N_s*yn;
            dzds=N_s*zn;
            dxdt=N_t*xn;
            dydt=N_t*yn;
            dzdt=N_t*zn;
            DetJ =dxdr.*dyds.*dzdt + dxdt .*dydr.*dzds +...
                dxds.*dydt.*dzdr - dxdt .*dyds.*dzdr -...
                dxdr.*dydt.*dzds - dxds .*dydr.*dzdt;
             if any(DetJ<0)
%                 warning(sprintf('DetJ<0 at %d points !!! THEY ARE ELIMINATED !!!',sum(DetJ<0)))
             if ~any(deleted)
                 deleted(DetJ<0)=1;
             else
                 fnd=find(~deleted);
                 deleted(fnd(DetJ<0))=1;
             end
             xpix(DetJ<0)=[];
             ypix(DetJ<0)=[];
             zpix(DetJ<0)=[];
             xg(DetJ<0)=[];
             yg(DetJ<0)=[];
             zg(DetJ<0)=[];
             invJ(DetJ<0,:)=[];
   N=[0.125*(1-xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1+yg).*(1-zg),0.125*(1-xg).*(1+yg).*(1-zg),...
                0.125*(1-xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1+yg).*(1+zg),0.125*(1-xg).*(1+yg).*(1+zg)];
            N_r=[-0.125*(1-yg).*(1-zg),0.125*(1-yg).*(1-zg),0.125*(1+yg).*(1-zg),-0.125*(1+yg).*(1-zg),...
                -0.125*(1-yg).*(1+zg),0.125*(1-yg).*(1+zg),0.125*(1+yg).*(1+zg),-0.125*(1+yg).*(1+zg)];
            N_s=[-0.125*(1-xg).*(1-zg),-0.125*(1+xg).*(1-zg),0.125*(1+xg).*(1-zg),0.125*(1-xg).*(1-zg),...
                -0.125*(1-xg).*(1+zg),-0.125*(1+xg).*(1+zg),0.125*(1+xg).*(1+zg),0.125*(1-xg).*(1+zg)];
            N_t=[-0.125*(1-xg).*(1-yg),-0.125*(1+xg).*(1-yg),-0.125*(1+xg).*(1+yg),-0.125*(1-xg).*(1+yg),...
                0.125*(1-xg).*(1-yg),0.125*(1+xg).*(1-yg),0.125*(1+xg).*(1+yg),0.125*(1-xg).*(1+yg)];

            dxdr=N_r*xn;
            dydr=N_r*yn;
            dzdr=N_r*zn;
            dxds=N_s*xn;
            dyds=N_s*yn;
            dzds=N_s*zn;
            dxdt=N_t*xn;
            dydt=N_t*yn;
            dzdt=N_t*zn;
            DetJ =dxdr.*dyds.*dzdt + dxdt .*dydr.*dzds +...
                dxds.*dydt.*dzdr - dxdt .*dyds.*dzdr -...
                dxdr.*dydt.*dzds - dxds .*dydr.*dzdt;
             end
            invJ(:,1) = (dyds .* dzdt - dydt .* dzds) ./ DetJ;
            invJ(:,4) = -(dydr .* dzdt - dydt .* dzdr) ./ DetJ;
            invJ(:,7)  = (dydr .* dzds - dyds .* dzdr) ./ DetJ;

            invJ(:,2)  = -(dxds .* dzdt - dxdt .* dzds) ./ DetJ;
            invJ(:,5)  = (dxdr .* dzdt - dxdt .* dzdr) ./ DetJ;
            invJ(:,8)  = -(dxdr .* dzds - dxds .* dzdr) ./ DetJ;

            invJ(:,3)  = (dxds .* dydt - dxdt .* dyds) ./ DetJ;
            invJ(:,6)  = -(dxdr .* dydt - dxdt .* dydr) ./ DetJ;
            invJ(:,9)  = (dxdr .* dyds - dxds .* dydr) ./ DetJ;

            xp=N*xn;
            yp=N*yn;
            zp=N*zn;

            dxg=invJ(:,1).*(xpix-xp)+invJ(:,2).*(ypix-yp)+invJ(:,3).*(zpix-zp);
            dyg=invJ(:,4).*(xpix-xp)+invJ(:,5).*(ypix-yp)+invJ(:,6).*(zpix-zp);
            dzg=invJ(:,7).*(xpix-xp)+invJ(:,8).*(ypix-yp)+invJ(:,9).*(zpix-zp);


            res=dxg'*dxg+dyg'*dyg+dzg'*dzg;

            xg=xg+dxg;
            yg=yg+dyg;
            zg=zg+dzg;


            if check
                figure
                subplot(3,2,1)
                scatter(xg,yg,50+0*yg)
                hold on
                rect1=line([-1,1,1,-1,-1],[-1,-1,1,1,-1]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',2)
                subplot(3,2,3)
                scatter(yg,zg,50+0*yg)
                hold on
                rect1=line([-1,1,1,-1,-1],[-1,-1,1,1,-1]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',5)
                subplot(3,2,1)
                scatter(zg,xg,50+0*yg)
                hold on
                rect1=line([-1,1,1,-1,-1],[-1,-1,1,1,-1]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',2)
                subplot(3,2,2)
                scatter(xp,yp,50+0*yg,'o')
                hold on
                scatter(xpix,ypix,50+0*yg,'x')
                rect1=line([xn;xn(1)],[yn;yn(1)]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',2)
                subplot(3,2,4)
                scatter(yp,zp,50+0*yg,'o')
                hold on
                scatter(ypix,zpix,50+0*yg,'x')
                rect1=line([yn;yn(1)],[zn;zn(1)]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',2)
                subplot(3,2,6)
                scatter(zp,xp,50+0*yg,'o')
                hold on
                scatter(zpix,xpix,50+0*yg,'x')
                rect1=line([zn;zn(1)],[xn;xn(1)]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',2)

                pause
            end
        end
        
        wg=~(abs(xg)>1|abs(yg)>1|abs(zg)>1);
if any(deleted)
    tmp=zeros(npix,1);
    tmp(~deleted)=xg;
    mp=mean(xg);
    xg=tmp;xg(deleted==1)=mp;
    tmp(~deleted)=yg;
    mp=mean(yg);
    yg=tmp;yg(deleted==1)=mp;
    tmp(~deleted)=zg;
    mp=mean(zg);
    zg=tmp;zg(deleted==1)=mp;
    tmp(~deleted)=wg;
    wg=tmp;
end
end
if prod(nb_sub_cell)>1&&nb_gauss_point>0
    dx=2/nb_sub_cell(1);
    dy=2/nb_sub_cell(2);
    dz=2/nb_sub_cell(2);
        xc=(-1+dx/2):dx:(1-dx/2);
        yc=(-1+dy/2):dy:(1-dy/2);
        zc=(-1+dz/2):dz:(1-dz/2);
    [yc,xc,zc]=meshgrid(yc,xc,zc);
        % quad sub cells

        xg=repmat(xc(:)',numel(xg),1)+repmat(dx*xg/2,1,numel(xc));
        yg=repmat(yc(:)',numel(yg),1)+repmat(dy*yg/2,1,numel(yc));
        zg=repmat(zc(:)',numel(zg),1)+repmat(dz*zg/2,1,numel(zc));
        xg=xg(:);
        yg=yg(:);
        zg=zg(:);
        wg=repmat(wg/numel(xc),numel(xc),1);
end



end
