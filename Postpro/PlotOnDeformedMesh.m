function PlotOnDeformedMesh(xn,yn,unmap,vnmap,dmap,titre,lim,fac,filexp,epsflag,invcb)
if nargin<10, epsflag=false;end
if nargin<11, invcb=false;end
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
axis xy;
axis off;


[Yo,Xo]=meshgrid(yn(1):yn(length(yn)),xn(1):xn(length(xn)));
[Yn,Xn]=meshgrid(yn-0.5,xn-0.5);
mail2=mesh((Xo+fac*unmap)',(Yo+fac*vnmap)',0*Yo',dmap','Parent',axes1);
    
      set(mail2,'EdgeColor','none',...
        'FaceColor','interp',...
        'Marker','none')
    if max(length(xn),length(yn))>20
    linew=1;
    else
        linew=1.5;
    end
 hline1=  line( (Xo(xn,:)+fac*unmap(xn,:))',(Yo(xn,:)+fac*vnmap(xn,:))',1+0*Yo(xn,:)','Color','k','LineWidth',linew,'Parent',axes1);
  hline2= line( (Xo(:,yn)+fac*unmap(:,yn)),(Yo(:,yn)+fac*vnmap(:,yn)),1+0*Yo(:,yn),'Color','k','LineWidth',linew,'Parent',axes1);
  
  
  
  
  
  colorbar('peer',axes1,'FontSize',20);
    title([titre,' on deformed mesh (x',num2str(fac),')'],'FontSize',24);
         print ('-djpeg', fullfile('FIG',[filexp,'.jpg']));
         if epsflag

         print ('-depsc', fullfile('FIG',[filexp,'.eps']));
% unix(['epstopdf FIG/',filexp,'.eps']);
% unix(['pdftops FIG/',filexp,'.pdf']);
% unix(['ps2eps -f FIG/',filexp,'.ps']);
% unix(['rm  FIG/',filexp,'.pdf']);
% unix(['rm  FIG/',filexp,'.ps']);
%  unix(['jpeg2eps FIG/',filexp,'.jpg >',filexp,'.eps']);
         end
end
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    