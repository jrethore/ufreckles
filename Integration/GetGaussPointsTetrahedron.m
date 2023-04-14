function [xg,yg,zg,wg]=GetGaussPointsTetrahedron(nb_gauss_point,nb_sub_cell,xn,yn,zn,xpix,ypix,zpix)
if nargin==1, nb_sub_cell=1;end
check=0;
switch nb_gauss_point
    case 1
        xg=[1/4];
        yg=[1/4];
        zg=[1/4];
        wg=[1/6];
    case 0
        
        npix=numel(xpix);
        deleted=0*xpix;
        xg=0*xpix;
        yg=0*ypix;
        zg=0*zpix;
        res=1;
        while res>1.e-6
            
            N=[1-xg-yg-zg,xg,yg,zg];
            N_r=[-1+0*xg,1+0*xg,0*yg,0*zg];
            N_s=[-1+0*xg,0*xg,1+0*yg,0*zg];
            N_t=[-1+0*xg,0*xg,0*yg,1+0*zg];
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
                N=[1-xg-yg-zg,xg,yg,zg];
                N_r=[-1+0*xg,1+0*xg,0*yg,0*zg];
                N_s=[-1+0*xg,0*xg,1+0*yg,0*zg];
                N_t=[-1+0*xg,0*xg,0*yg,1+0*zg];
                
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
        
        wg=~(xg<0|yg<0|zg<0|1-xg-yg-zg<0);
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
if prod(nb_sub_cell)>1
    dx=1/nb_sub_cell(1);
    dy=1/nb_sub_cell(2);
    dz=1/nb_sub_cell(3);
    xg=dx/2:dx:1;
    yg=dy/2:dy:1;
    zg=dz/2:dz:1;
    dw=1/prod(nb_sub_cell);
    [yg,xg,zg]=meshgrid(yg,xg,zg);
    keep=(1-xg(:)-yg(:)-zg(:))>0;
    xg=xg(keep);
    yg=yg(keep);
    zg=zg(keep);
    wg=ones(numel(xg),1)*dw;
    if check
        figure
        plot3(xg,yg,zg,'rx')
    end
    
end
end
