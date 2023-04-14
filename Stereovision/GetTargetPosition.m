function [xp,yp,calibok]=GetTargetPosition(nmod,im0,U,V,gim)
persistent dist xyi dmin

sizeim0=size(im0);
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
ngrid=param.grid_size;

nim0=255-im0;
nx=ngrid(1);
ny=ngrid(2);
calibok=0;
ntry=0;
set(gim,'Cdata',im0')

%while ~calibok&&ntry<5
    calibok=1;
if 1||isempty(xyi)
title('Pick up the lower left dot')
    [xo,yo] = ginput(1);
    if isempty(dist)
title('Pick up a neighboor')
    [xn,yn] = ginput(1);
dist=abs(xo-xn+1i*(yo-yn));
    end
else
    xo=xyi(1)+U;
    yo=xyi(2)+V;
end
%gini= plot(xo,yo,'sm','MarkerSize',20)

%xo=316;yo=964;
%%
%nx=9;
%ny=9;
[indy,indx]=meshgrid(1:ny,1:nx);

%if isempty(dmin),
    dmin=round(0.2*dist);
%end
res=Inf;
conv=1.e-3;
xo=round(xo);
yo=round(yo);
in=-5:5;
xin=in+xo;
yin=in+yo;
    try 
sub=nim0(xin,yin);
    catch
        calibok=0;
ntry=ntry+1;
%        continue
    end
thres=0.9*mean(sub(:));
nim0s=nim0>thres;
subold=0;
while res>conv
    in=-dmin:dmin;
    xo=round(xo);
    yo=round(yo);
    xin=in+xo;
    yin=in+yo;
    try 
    sub=nim0s(xin,yin);
    catch
        calibok=0;
ntry=ntry+1;
        break
        
    end
    [yin,xin]=meshgrid(in,in);
    dx=xin.*sub;
    dy=yin.*sub;

    
    sub=sum(sub(:));
    
    dx=sum(dx(:))/sub;
    dy=sum(dy(:))/sub;
    res=abs(subold-sub)/dmin^2;
    subold=sub;
    dmin=dmin+1;
    xo=xo+dx;
    yo=yo+dy;
    if dmin/round(0.2*dist)>10||dmin/round(0.2*dist)<.1
                calibok=0;
ntry=ntry+1;
        break
    end
end
in=-2*dmin:2*dmin;
    xo=round(xo);
    yo=round(yo);
    xin=in+xo;
    yin=in+yo;
    sub=nim0(xin,yin);
in=(-dmin:dmin)+2*dmin+1;
select=0*sub;
select(in,in)=1;
thres=0.5*(mean(sub(select==1))+mean(sub(select==0)));

nim0s=nim0>thres;

%dmin=dmin-1;
lx=dist;
ly=dist;
xx=1;xy=0;
yx=0;yy=1;
dmin=round(1.25*dmin);
conv=0.5;
Xp=xo+(indx-1)*lx;
Yp=yo+(indy-1)*ly;
dXYp=max(indx,indy);
[dXYps,ordre]=sort(dXYp(:));
if calibok
% figure
% imagesc(im0)
% colormap('gray')
% hold on
%gini= plot(Xp,Yp,'bo','MarkerSize',10);
donex=zeros(nx,ny-1);
doney=zeros(nx-1,ny);
%          figure
%          imagesc(nim0)
%          colormap('gray')
%          hold on
for ii=1:length(ordre)
    [ix,iy]=ind2sub([nx,ny],ordre(ii));
    if (ix>1)||(iy>1)
        reso=0;
        res=Inf;
        xo=Xp(1,1)+(ix-1)*lx*xx+(iy-1)*ly*yx;
        yo=Yp(1,1)+(ix-1)*lx*xy+(iy-1)*ly*yy;
        
        
        in=-2*dmin:2*dmin;
    xo=round(xo);
    yo=round(yo);
    xin=in+xo;
    yin=in+yo;
    sub=nim0(xin,yin);
in=(-dmin:dmin)+2*dmin+1;
select=0*sub;
select(in,in)=1;
thres=0.5*(mean(sub(select==1))+mean(sub(select==0)));



        
%         plot(yoo,xoo,'w+','MarkerSize',10,'LineWidth',2)
%         plot(yo,xo,'rx','MarkerSize',10,'LineWidth',2)
iter=0;
        while (res>conv)&&(abs(res-reso)>0)&&(iter<100)
            reso=res;
            in=-dmin:dmin;
            xo=round(xo);
            yo=round(yo);
            xin=min(max(1,in+xo),size(nim0s,1));
            yin=min(max(1,in+yo),size(nim0s,2));
            sub=nim0(xin,yin)>thres;
            [yin,xin]=meshgrid(in,in);
            dx=xin.*sub;
            dy=yin.*sub;
            
            
            sub=sum(sub(:));
            
            dx=sum(dx(:))/sub;
            dy=sum(dy(:))/sub;
            res=abs(dx+1i*dy);
            xo=xo+dx;
            yo=yo+dy;
            iter=iter+1;
%             if iter>5
%                 keyboard
%             end
%             plot(yo,xo,'bx','MarkerSize',10)
        end
%        pause
%            plot(yo,xo,'bo','MarkerSize',10)
%              pause(.1)
        Xp(ix,iy)=xo;
        Yp(ix,iy)=yo;
    end
    if (ix>1)||(iy>1)
        donex(ix,max(1,iy-1))=1;
        doney(max(1,ix-1),iy)=1;
    end
    if any(donex(:))
        xx=diff(Xp,1,1);
        xy=diff(Yp,1,1);
        yx=diff(Xp,1,2);
        yy=diff(Yp,1,2);
        
        xx=mean(xx(doney>0));
        xy=mean(xy(doney>0));
        yx=mean(yx(donex>0));
        yy=mean(yy(donex>0));
        lx=abs(xx+1i*xy);
        xx=xx/lx;
        xy=xy/lx;
        ly=abs(yx+1i*yy);
        yx=yx/ly;
        yy=yy/ly;
    end
end
ntry=ntry+1;
end
xp=Xp(:)-0.5*(sizeim0(1)-1);
yp=Yp(:)-0.5*(sizeim0(2)-1);
if isempty(xyi)
    xyi=[Xp(1),Yp(1)];
end
if any(isnan(xp))||any(isnan(yp))
    calibok=0;
end
if ~calibok
title('Detected dots: failure ! Image will be ignored....')
    display('Calibration failed')
gp=plot(Xp,Yp,'rx','MarkerSize',10,'LineWidth',2);
else
    title('Detected dots')

gp=plot(Xp,Yp,'gx','MarkerSize',10,'LineWidth',2);

end


    pause(1)
delete(gp)
try delete(gini),catch ;end
end