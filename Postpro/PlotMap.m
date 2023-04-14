function PlotMap(dmap,titre,lim,filexp,epsflag,invcb)

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
if isempty(lim)
axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
    'XTick',zeros(1,0),...
    'FontSize',20);
else
axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
    'XTick',zeros(1,0),...
    'FontSize',20,'CLim',lim);
end
box('on');
hold 'all';
image(dmap','Parent',axes1,'CDataMapping','scaled');
axis off;
axis xy;
axis image;
colorbar('peer',axes1,'FontSize',20);
title(titre,'FontSize',24);

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







