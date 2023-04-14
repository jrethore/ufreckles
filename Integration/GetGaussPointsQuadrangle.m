function [xg,yg,wg]=GetGaussPointsQuadrangle(nb_gauss_point,nb_sub_cell,xn,yn,xpix,ypix)
if nargin==1, nb_sub_cell=1;end
check=0;
checkc=0;
switch nb_gauss_point
    case 0
        xg=0*xpix;
        yg=0*ypix;
        res=1;
        while res>1.e-6
            
            N=[0.25*(1-xg).*(1-yg),0.25*(1+xg).*(1-yg),0.25*(1+xg).*(1+yg),0.25*(1-xg).*(1+yg)];
            N_r=[-0.25*(1-yg),0.25*(1-yg),0.25*(1+yg),-0.25*(1+yg)];
            N_s=[-0.25*(1-xg),-0.25*(1+xg),0.25*(1+xg),0.25*(1-xg)];
            dxdr=N_r*xn;
            dydr=N_r*yn;
            dxds=N_s*xn;
            dyds=N_s*yn;
            detJ=(dxdr.*dyds-dydr.*dxds);
            invJ=[dyds./detJ,-dxds./detJ,-dydr./detJ,dxdr./detJ];
            xp=N*xn;
            yp=N*yn;
            dxg=invJ(:,1).*(xpix-xp)+invJ(:,2).*(ypix-yp);
            dyg=invJ(:,3).*(xpix-xp)+invJ(:,4).*(ypix-yp);
            
            res=dxg'*dxg+dyg'*dyg;
            if check
                figure
                subplot(1,2,1)
                scatter(xg,yg,50+0*yg)
                hold on
                rect1=line([-1,1,1,-1,-1],[-1,-1,1,1,-1]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',2)
                subplot(1,2,2)
                scatter(xp,yp,50+0*yg,'o')
                hold on
                scatter(xpix,ypix,50+0*yg,'x')
                rect1=line([xn;xn(1)],[yn;yn(1)]);
                set(rect1,...
                    'Color','red',...
                    'LineStyle','-',...
                    'LineWidth',2)
                display(sprintf('Residual %f\n',res));
                pause
            end
            xg=xg+dxg;
            yg=yg+dyg;
        end
        wg=~(abs(xg)>1|abs(yg)>1);
        %         ind=find(wg);
        %         xg=xg(ind);yg=yg(ind);wg=wg(ind);
    case 1
        xg=0;
        yg=0;
        wg = 4;
    case 4
        
        xg=[-1   1 -1    1]'/sqrt(3);
        yg=[-1 -1    1   1]'/sqrt(3);
        wg = [1.,1.,1.,1.]';
    case 16
        
        a16 = 0.861136311594053;
        b16 = 0.339981043584856;
        xg = [-a16, -b16, b16, a16, -a16, -b16, b16, a16, -a16, -b16, b16, a16, -a16, -b16, b16, a16]';
        yg = [-a16, -a16, -a16, -a16, -b16, -b16, -b16, -b16, b16, b16, b16, b16, a16, a16, a16, a16]';
        pe2 = 0.347854845137454 * 0.347854845137454;
        pf2 = 0.652145154862546 * 0.652145154862546;
        pef = 0.347854845137454 * 0.652145154862546;
        wg = [ pe2, pef,  pef, pe2, pef, pf2, pf2, pef, pef, pf2, pf2, pef, pe2, pef, pef, pe2]';
    case 64
        x1=       [ -.9602898564975362,
            -.7966664774136267,
            -.5255324099163290,
            -.1834346424956498,
            .1834346424956498,
            .5255324099163290,
            .7966664774136267,
            .9602898564975362];
        [yg,xg]=meshgrid(x1,x1);
        
        w1 = [.1012285362903763,
            .2223810344533745,
            .3137066458778873,
            .3626837833783620,
            .3626837833783620,
            .3137066458778873,
            .2223810344533745,
            .1012285362903763];
        wg=w1'*w1;
        xg=xg(:);
        yg=yg(:);
        wg=wg(:);
    case 100
        
        x1=       [ -0.97390653,
            -0.86506337,
            -0.67940957,
            -0.43339539,
            -0.14887434,
            0.14887434,
            0.43339539,
            0.67940957,
            0.86506337,
            0.97390653];
        [yg,xg]=meshgrid(x1,x1);
        
        
        w1=[0.06667134,
            0.14945135,
            0.21908636,
            0.26926672,
            0.29552422,
            0.29552422,
            0.26926672,
            0.21908636,
            0.14945135,
            0.06667134];
        
        wg=w1'*w1;
        xg=xg(:);
        yg=yg(:);
        wg=wg(:);
        
end
if prod(nb_sub_cell)>1&&nb_gauss_point>0
    dx=2/nb_sub_cell(1);
    dy=2/nb_sub_cell(2);
    xc=(-1+dx/2):dx:(1-dx/2);
    yc=(-1+dy/2):dy:(1-dy/2);
    [yc,xc]=meshgrid(yc,xc);
    if nb_gauss_point==4
        % quad sub cells
        xgc=[-1   1 -1    1]'/sqrt(3);
        ygc=[-1 -1    1   1]'/sqrt(3);
        wgc = [1.,1.,1.,1.]';
        xg=repmat(xc(:)',4,1)+dx*repmat(xgc,1,prod(nb_sub_cell))/2;
        yg=repmat(yc(:)',4,1)+dy*repmat(ygc,1,prod(nb_sub_cell))/2;
        wg=repmat(wgc,1,prod(nb_sub_cell))/prod(nb_sub_cell);
        xg=xg(:);
        yg=yg(:);
        wg=wg(:);
        %triangle subcells
        % [xgt,ygt,wgt]=GetGaussPointsTriangle(1);
        % xgo=repmat(0,4*length(xgt),1);
        % ygo=repmat(0,4*length(xgt),1);
        % wgo=repmat(2*wgt,4,1);
        % for ia=0:3
        %     alph=ia*pi/2+pi/4;
        %     xgo(ia*length(xgt)+(1:length(xgt)))=xgt*cos(alph)-ygt*sin(alph);
        %     ygo(ia*length(xgt)+(1:length(xgt)))=xgt*sin(alph)+ygt*cos(alph);
        % end
        %          xg=repmat(xc(:)',numel(xgo),1)+repmat(dx*xgo/2,1,numel(xc));
        %          yg=repmat(yc(:)',numel(ygo),1)+repmat(dy*ygo/2,1,numel(yc));
        %          xg=xg(:);
        %          yg=yg(:);
        %          wg=repmat(wgo/numel(xc),numel(xc),1);
        
    else
        % quad sub cells
        
        xg=repmat(xc(:)',numel(xg),1)+repmat(dx*xg/2,1,numel(xc));
        yg=repmat(yc(:)',numel(yg),1)+repmat(dy*yg/2,1,numel(yc));
        xg=xg(:);
        yg=yg(:);
        wg=repmat(wg/numel(xc),numel(xc),1);
    end
    
    if checkc
        figure
        plot(xg,yg,'rx')
        hold on
        plot(xc,yc,'bo')
        rect1=line([-1,1,1,-1,-1],[-1,-1,1,1,-1]);
        set(rect1,...
            'Color','red',...
            'LineStyle','-',...
            'LineWidth',2)
        title(sprintf('Gauss point and sub cell centers, S=%g',sum(wg)))
    end
    
    
    
    
end
end
