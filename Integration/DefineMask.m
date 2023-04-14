function [mask]=DefineMask(image);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Defines masks as polys and circles
%
% Called by correli_q4.m
%
% SR don't remember when
% HF 16/08/06
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[nx ny]=size(image); % image size
%
nptmax=30;    % maximum number of points in a poly
npolymax=20;  % maximum number of exclusion polys
xpoly=zeros(nptmax,npolymax); % initialize x coordinates of poly points
ypoly=zeros(nptmax,npolymax); % initialize x coordinates of poly points
%
figure
set(gcf,'Name','MASK MANAGER');
set(gca,'FontName','Arial','FontSize',12);
colormap(gray)
imagesc(image) % draw picture...
axis equal     % ... with square pixels
hold on;
%
rect1=rectangle; % draw ROI
% set(rect1,'position',[yroi(1),xroi(1),yroi(2)-yroi(1),xroi(2)-xroi(1)],...
%    'EdgeColor','yellow',...
%    'LineStyle','--',...
%    'LineWidth',2)
%
npoly=0; ncirce=0; ncirci=0; finish=0; cire = []; ciri = [];
npt=0;% initialize
%
while (finish==0) 
   %
   item=menu('Mask Definition','Define exclusion polygon','Remove polygon',...
      'Define exclusion circle','Remove exclusion circle',...
      'Define inclusion circle','Remove inclusion circle','Redraw','Exit');
   %
   switch item
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 1
      % Define exclusion polygon
      if (npoly<=npolymax)
         npoly=npoly+1;
         title('Select an arbitrary number of points and press RETURN when finished');
         [yp xp]=ginput;
         xp=round(xp);xp=min(xp,nx);xp=max(xp,1);
         yp=round(yp);yp=min(yp,ny);yp=max(yp,1);
         npt(npoly)=length(xp)+1;
         xp(npt(npoly))=xp(1);yp(npt(npoly))=yp(1);
         xpoly(1:npt(npoly),npoly)=xp;ypoly(1:npt(npoly),npoly)=yp;
         col=rand(1,3);
         fill(yp,xp,col);
         clear xp yp;
      else
         menu(sprintf('Maximum number of exclusion polygons reached: \n Enough is enough'),'OK');
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 2
      % remove polygon
      if (npoly>0)
         found=0;
         while (found==0)
            title('Select a point within the polygon you want to remove');
            [y0 x0]=ginput(1);
            ipoly=0;
            while (found==0)&(ipoly<npoly)
               ipoly=ipoly+1;
               found=inpolygon(y0,x0,ypoly(1:npt(ipoly),ipoly),xpoly(1:npt(ipoly),ipoly));
            end
            if (found==0)
               menu(sprintf('Sorry, polygon not found! \n Please try again'),'OK');
            end
         end
         for jpoly=ipoly+1:npoly
            npt(jpoly-1)=npt(jpoly);
            xpoly(:,jpoly-1)=xpoly(:,jpoly);
            ypoly(:,jpoly-1)=ypoly(:,jpoly);
         end
         npoly=npoly-1;
         DrawPolyCirc(image,xroi,yroi,npoly,npt,xpoly,ypoly,ncirce,cire,ncirci,ciri);
      else
         menu(sprintf('No exclusion polygon'),'OK');
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 3
      % Define exclusion circle
      ncirce=ncirce+1;
      title('Select two diametricaly opposed points');
      [yp xp]=ginput(2);
      xcenter=(xp(1)+xp(2))*.5;ycenter=(yp(1)+yp(2))*.5;
      d2=(yp(2)-yp(1))^2+(xp(2)-xp(1))^2;r1=sqrt(d2)*.5;
      cire(1:3,ncirce)=[xcenter ycenter r1]';
      col=rand(1,3);
      rectangle('Position',...
         [(cire(2,ncirce)-cire(3,ncirce)),(cire(1,ncirce)-cire(3,ncirce)),...
            (2*cire(3,ncirce)),(2*cire(3,ncirce))],...
         'Curvature',[1,1],...
         'FaceColor',col);
      clear xp yp;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 4
      % remove exclusion circle
      if (ncirce>0)
         found=0;
         while (found==0)
            title('Select a point within the circle you want to remove');
            [y0 x0]=ginput(1);
            icirce=0;
            while (found==0)&(icirce<ncirce)
               icirce=icirce+1;
               d2=(x0-cire(1,icirce))^2+(y0-cire(2,icirce))^2-(cire(3,icirce))^2;
               found=(d2<0);
            end
            if (found==1)
               for jcirce=icirce+1:ncirce
                  cire(:,jcirce-1)=cire(:,jcirce);
               end
               ncirce=ncirce-1;
            else
               respon = menu(sprintf('Sorry, exclusion circle not found!'),'OK','EXIT');
               if (respon==2) 
                  break 
               end
            end
         end
         DrawPolyCirc(image,xroi,yroi,npoly,npt,xpoly,ypoly,ncirce,cire,ncirci,ciri);
      else
         menu(sprintf('No exclusion circle'),'OK');
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 5
      % Define inclusion circle
      ncirci=ncirci+1;
      title('Select two diametricaly opposed points. BEWARE: inclusion boundary depicted');
      [yp xp]=ginput(2);
      xcenter=(xp(1)+xp(2))*.5;ycenter=(yp(1)+yp(2))*.5;
      d2=(yp(2)-yp(1))^2+(xp(2)-xp(1))^2;r1=sqrt(d2)*.5;
      ciri(1:3,ncirci)=[xcenter ycenter r1]';
      col=rand(1,3);
      rectangle('Position',...
         [(ciri(2,ncirci)-ciri(3,ncirci)),(ciri(1,ncirci)-ciri(3,ncirci)),...
            (2*ciri(3,ncirci)),(2*ciri(3,ncirci))],...
         'Curvature',[1,1],...
         'EdgeColor',col,...  % exclusion zone not depicted
         'LineStyle','--',... % only its contour
         'LineWidth',2)
      clear xp yp;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 6
      % remove inclusion circle
      if (ncirci>0)
         found=0;
         while (found==0)
            icirci=0;
            title('Select a point within the circle you want to remove');
            [y0 x0]=ginput(1);
            while (found==0)&(icirci<ncirci)
               icirci=icirci+1;
               d2=(x0-ciri(1,icirci))^2+(y0-ciri(2,icirci))^2-(ciri(3,icirci))^2;
               found=(d2<0);
            end
            if (found==1)
               for jcirci=icirci+1:ncirci
                  ciri(:,jcirci-1)=ciri(:,jcirci);
               end
               ncirci=ncirci-1;
            else
               respon = menu(sprintf('Sorry, inclusion circle not found!'),'OK','EXIT');
               if (respon==2) 
                  break 
               end
            end
         end
         DrawPolyCirc(image,xroi,yroi,npoly,npt,xpoly,ypoly,ncirce,cire,ncirci,ciri);
      else
         menu(sprintf('No inclusion circle'),'OK');
      end
      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 7
      % Redraw
      DrawPolyCirc(image,xroi,yroi,npoly,npt,xpoly,ypoly,ncirce,cire,ncirci,ciri);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 8
      % Exit
      mask=ones(nx,ny);
      [Y X]=meshgrid((1:ny),(1:nx));
%      mask(xroi(1):xroi(2),yroi(1):yroi(2))=repmat(1,xroi(2)-xroi(1)+1,yroi(2)-yroi(1)+1); % ROI considered
      for ipoly=1:npoly
         miny = max(1,floor(min(ypoly(:,ipoly))));
         minx = max(1,floor(min(xpoly(:,ipoly))));
         maxy = min(ny,ceil(max(ypoly(:,ipoly))));
         maxx = min(nx,ceil(max(xpoly(:,ipoly))));
         IN=~inpolygon(Y(minx:maxx,miny:maxy),X(minx:maxx,miny:maxy),ypoly(1:npt(ipoly),ipoly),xpoly(1:npt(ipoly),ipoly));
         mask(minx:maxx,miny:maxy)=mask(minx:maxx,miny:maxy)&IN;
      end
      for icirce=1:ncirce
         miny = max(1,floor(min(cire(2,icirce)-cire(3,icirce))));
         minx = max(1,floor(min(cire(1,icirce)-cire(3,icirce))));
         maxy = min(ny,ceil(max(cire(2,icirce)+cire(3,icirce))));
         maxx = min(nx,ceil(max(cire(1,icirce)+cire(3,icirce))));
         d2=(X(minx:maxx,miny:maxy)-cire(1,icirce)).^2+(Y(minx:maxx,miny:maxy)-cire(2,icirce)).^2-(cire(3,icirce))^2;
         mask(minx:maxx,miny:maxy)=mask(minx:maxx,miny:maxy)&(d2>0);
      end
      for icirci=1:ncirci
         miny = max(1,floor(min(ciri(2,icirci)-ciri(3,icirci))));
         minx = max(1,floor(min(ciri(1,icirci)-ciri(3,icirci))));
         maxy = min(ny,ceil(max(ciri(2,icirci)+ciri(3,icirci))));
         maxx = min(nx,ceil(max(ciri(1,icirci)+ciri(3,icirci))));
         d2=(X(xroi(1):xroi(2),yroi(1):yroi(2))-ciri(1,icirci)).^2+(Y(xroi(1):xroi(2),yroi(1):yroi(2))-ciri(2,icirci)).^2-(ciri(3,icirci))^2;
         mask(xroi(1):xroi(2),yroi(1):yroi(2))=mask(xroi(1):xroi(2),yroi(1):yroi(2))&(d2<0);
      end
      close
      figure
      set(gcf,'Name','MASK MANAGER');
      set(gca,'FontName','Arial','FontSize',12);
      colormap(gray)
      imagesc(uint16(double(image).*mask))
      axis equal
      hold on;
      %index = find(mask(:)==1);
      %maskd = ones(nx,ny);
      %masdd(index) = NaN;
      %imagesc(maskd)
      %clear maskd
      %pause
      %imagesc(mask)
%       rect1=rectangle;
%       set(rect1,'position',[yroi(1),xroi(1),yroi(2)-yroi(1),xroi(2)-xroi(1)],...
%          'EdgeColor','yellow',...
%          'LineStyle','--',...
%          'LineWidth',2)
      nver = version;
      %
         [c,h] = contour(mask,[0.5 0.5]);
         set(h,'EdgeColor','red','LineWidth',2);
 %       title(sprintf('Strike any key to continue'));
%       pause
%      close; pause(0.1);
      title(sprintf('Mask'));
      menu(sprintf('Number of exclusion polygons %g\n Number of exclusion circles  %g\n Number of exclusion circles  %g',...
         npoly,ncirce,ncirci),'OK');
      finish=1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end
end
end