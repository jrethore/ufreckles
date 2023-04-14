function DeformedMesh(xn,yn,unmap,vnmap,lim,fac,filexp,epsflag,invcb)
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
axis xy;
axis off;

if length(xn)==size(unmap,1)
    
    [Yo,Xo]=meshgrid(yn,xn);

    mail2=mesh((Xo+fac*unmap)',(Yo+fac*vnmap)',0*Yo',abs(unmap+i*vnmap)','Parent',axes1);
    
      set(mail2,'EdgeColor','k',...
        'FaceColor','interp',...
        'Marker','none');
    
else
    if (size(unmap,1)==(max(xn)-min(xn)+1))&&(size(unmap,2)==(max(yn)-min(yn)+1))
        [Yo,Xo]=meshgrid(yn(1):yn(length(yn)),xn(1):xn(length(xn)));
        xn=xn-xn(1)+1;
        yn=yn-yn(1)+1;
        mail2=mesh((Xo+fac*unmap)',(Yo+fac*vnmap)',0*Yo',abs(unmap+i*vnmap)','Parent',axes1);

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
    else
        xo=[ceil(min(xn));floor(max(xn))-1];
        yo=[ceil(min(yn));floor(max(yn))-1];
        [Yo,Xo]=meshgrid(yo(1):yo(length(yo)),xo(1):xo(length(xo)));
        mail2=mesh((Xo+fac*unmap)',(Yo+fac*vnmap)',0*Yo',abs(unmap+i*vnmap)','Parent',axes1);

        set(mail2,'EdgeColor','none',...
            'FaceColor','interp',...
            'Marker','none')
        if max(length(xn),length(yn))>20
            linew=1;
        else
            linew=1.5;
        end
        xe=round(xn)-xo(1)+1;
        ye=round(yn)-yo(1)+1;
        xe=max(1,min(xe,max(xo-xo(1)+1)));
        ye=max(1,min(ye,max(yo-yo(1)+1)));
        ind=sub2ind(size(unmap),xe,ye);
        hline1=  plot( xn+fac*unmap(ind),yn+fac*vnmap(ind),'x','Color','k','LineWidth',linew,'MarkerSize',10*linew,'Parent',axes1);

    end

end



colorbar('peer',axes1,'FontSize',20);
title(['Deformed mesh (x',num2str(fac),')'],'FontSize',24);
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
























