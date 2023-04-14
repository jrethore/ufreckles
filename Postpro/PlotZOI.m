function PlotZOI(roi,x,y,dmap,titre,filexp,epsflag)

if nargin<7, epsflag=false;end


figure % CREER UNE FONCTION
colormap(gray)
imagesc(dmap')
axis xy
axis image
axis off;
hold on;
rect1=rectangle;
set(rect1,'position',[roi(1),roi(3),roi(2)-roi(1),roi(4)-roi(3)],...
    'EdgeColor','yellow',...
    'LineStyle','--',...
    'LineWidth',2)
%            [c,h] = contour(mask,[0.5 0.5]);
%            set(h,'EdgeColor','blue','LineWidth',2);
%    plot(x,y,'w+','LineWidth',1.5,'MarkerSize',10)
if ~(length(x)==length(y))
    [mrx mry]=ndgrid(x,y);
    mail1=mesh(mrx',mry',0*mry');
    set(mail1,'EdgeColor','white',...
        'FaceColor','none');%,...
end
title(titre,'FontSize',24);

if epsflag
    print ('-depsc', fullfile('FIG',[filexp]));
end
print ('-djpeg', fullfile('FIG',[filexp]));

end


