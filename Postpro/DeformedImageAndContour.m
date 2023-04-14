function DeformedImageAndContour(xn,yn,im,filexp,epsflag,invcb)
fac=1;
if nargin<5, epsflag=false;end
if nargin<6, invcb=false;end
figure
cb=colormap(hot);
cb=cb(:,[3,2,1]);
colormap(cb);
if invcb
   cb=flipud(cb);
   colormap(cb);
end
             colormap(gray)
             imagesc(im')

hold on;
axis equal
axis xy;
axis off;
% [Yi,Xi]=meshgrid(yn,xn);
% h=contour(Xi',Yi',lso',[0 0],'Color','blue','LineWidth',2);
   h=plot(xn,yn,'b.','LineWidth',1);
    title(['Detected contour'],'FontSize',24);
         print ('-djpeg', fullfile('FIG',[filexp,'.jpg']));
         if epsflag
        print ('-depsc2', fullfile('FIG',[filexp,'.eps']));
         end
end
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    