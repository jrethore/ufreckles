function DeformedMesh25D(xn,yn,zn,dmap,lim,fac,filexp,epsflag,invcb)
if nargin<8, epsflag=false;end
if nargin<9, invcb=false;end
figure
cb=colormap(hot);
cb=cb(:,[3,2,1]);
colormap(cb);
if invcb
   cb=flipud(cb);
   colormap(cb);
end
if isempty(lim)
axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
    'XTick',zeros(1,0),...
    'FontSize',20);
else
axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
    'XTick',zeros(1,0),...
    'FontSize',20,'CLim',lim);
end
hold on;
axis equal
%axis xy;
axis off;


    
mail2=mesh(xn',yn',zn',dmap','Parent',axes1);
    
      set(mail2,'EdgeColor','k',...
        'FaceColor','interp',...
        'Marker','none','LineWidth',1.5)
  
  
  
  colorbar('peer',axes1,'FontSize',20);
    title(['Deformed mesh (x',num2str(fac),')'],'FontSize',24);
         print ('-djpeg', fullfile('FIG',[filexp,'.jpg']));
         if epsflag
        print ('-depsc', fullfile('FIG',[filexp,'.eps']));
         end
end
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    