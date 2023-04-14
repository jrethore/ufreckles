function DeformedImageAndMesh(xn,yn,unmap,vnmap,lim,im,filexp,epsflag,invcb)
fac=1;
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
             colormap(gray)
             imagesc(im')

hold on;
axis equal
axis xy;
axis off;

if length(xn)==size(unmap,1)
    
    [Yo,Xo]=meshgrid(yn,xn);

    mail2=mesh((Xo+fac*unmap)',(Yo+fac*vnmap)',0*Yo',abs(unmap+i*vnmap)','Parent',axes1);
    
      set(mail2,'EdgeColor','interp',...
        'FaceColor','none',...
        'Marker','none');
    
else
[Yo,Xo]=meshgrid(yn(1):yn(length(yn)),xn(1):xn(length(xn)));
%Xo=Xo-mean(Xo(:))+max(xn);
xn=xn-xn(1)+1;
yn=yn-yn(1)+1;
    if max(length(xn),length(yn))>20
    linew=1;
    else
        linew=1.5;
    end
linew=2;
%                    xn(1)=round(mean(xn));
%                    xn(2)=xn(1)+1;
if length(xn)==2
                   xn(1)=xn(2)-1;
else
                   yn(1)=yn(2)-1;
end

hline1=  line( (Xo(xn,:)+fac*unmap(xn,:))',(Yo(xn,:)+fac*vnmap(xn,:))',1+0*Yo(xn,:)','Color','b','LineWidth',linew);
  hline2= line( (Xo(:,yn)+fac*unmap(:,yn)),(Yo(:,yn)+fac*vnmap(:,yn)),1+0*Yo(:,yn),'Color','b','LineWidth',linew);
  
  
end
  
  
    title(['Deformed mesh'],'FontSize',24);
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
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    