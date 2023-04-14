function CreateLevelSets(xf,yf,sizeim,filnam)
check=1;
revers=std(yf)>std(xf);
revers=0;
dx=norm(sizeim,2);
tx=xf(length(xf))-xf(length(xf)-1);
ty=yf(length(yf))-yf(length(yf)-1);
tnorm=abs(tx+i*ty);
tx=tx/tnorm;
ty=ty/tnorm;
xf=[xf,xf(length(xf))+tx*dx];
yf=[yf,yf(length(yf))+ty*dx];
if revers
    yc=xf;
    xc=yf;
     yi=(1:sizeim(1));
xi=(1:sizeim(2));
   
else
    xi=(1:sizeim(1));
yi=(1:sizeim(2));
    xc=xf;
    yc=yf;

    
    
end
[Yi,Xi]=meshgrid(yi,xi);

yon=interp1(xf,yf,xi,'linear','extrap');
%yon=pchip(xf,yf,xi);

% figure
% plot(xf,yf,'o')
% hold on
% plot(xi,yon,'-')

nx=-gradient(yon);
Nnorm=repmat(abs(nx+i)',1,size(Yi,2));
Yon=repmat(yon',1,size(Yi,2));

crack=(Yi-Yon)./Nnorm;
if tx<0
    crack=-crack;
end
% front=cumsum(Nx,2)-cumsum(Ny,1);
% front0=front;%-((Xi-xf(length(xf)-1)).*Ny(round(xf(length(xf)-1)),round(yf(length(yf)-1)))-(Yi-yf(length(yf)-1)).*Nx(round(xf(length(xf)-1)),round(yf(length(yf)))));
% front=front0-front0(round(xf(length(xf)-1)),round(yf(length(yf)-1)));
front=(Xi-xf(length(xf)-1))*tx+(Yi-yf(length(yf)-1))*ty;

dx=4;
xo=1:dx:size(front,1);
        yo=1:dx:size(front,2);
        frontn=front(xo,yo);
        crackn=crack(xo,yo);
    frontn=LSOrtho(frontn,crackn,20*dx,dx);
    frontn=LSReinit(frontn,10*dx,dx);
   front=interp2(yo,xo,frontn,Yi,Xi,'*linear'); 
if check
x2=xi(1):(xi(length(xi))-xi(1))/20:xi(length(xi));
y2=yi(1):(yi(length(yi))-yi(1))/10:yi(length(yi));
nnorm=abs(nx+i);
Nx=repmat((nx./nnorm)',1,size(Yi,2));
Ny=repmat((1./nnorm)',1,size(Yi,2));

        figure
        cm=colormap(hot);
        cm=cm(:,[3,2,1]);
        colormap(cm);
        axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
            'XTick',zeros(1,0),...
            'FontSize',20);
        box('on');
        hold 'all';
        image(crack','Parent',axes1,'CDataMapping','scaled');
        [c,h1] = contour(crack',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
        set(h1,'EdgeColor','black','LineWidth',2);
        [c,h1] = contour(front',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
        set(h1,'EdgeColor','black','LineWidth',2);
quiver(x2,y2,Nx(round(x2),round(y2))',Ny(round(x2),round(y2))','Color','w')
quiver(x2,y2,Ny(round(x2),round(y2))',-Nx(round(x2),round(y2))','Color','w')
        axis off;
        axis xy;
        axis image;
        colorbar('peer',axes1,'FontSize',20);


        figure
        cm=colormap(hot);
        cm=cm(:,[3,2,1]);
        colormap(cm);
        axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
            'XTick',zeros(1,0),...
            'FontSize',20);
        box('on');
        hold 'all';
        image(front','Parent',axes1,'CDataMapping','scaled');
        [c,h1] = contour(crack',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
        set(h1,'EdgeColor','black','LineWidth',2);
        [c,h1] = contour(front',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
        set(h1,'EdgeColor','black','LineWidth',2);
quiver(x2,y2,Nx(round(x2),round(y2))',Ny(round(x2),round(y2))','Color','w')
quiver(x2,y2,Ny(round(x2),round(y2))',-Nx(round(x2),round(y2))','Color','w')

        axis off;
        axis xy;
        axis image;
        colorbar('peer',axes1,'FontSize',20);
end

save(filnam,'crack','front');



end