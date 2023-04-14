function PropagateLevelSets(dec,thp,imod,step,dx,dec2)
tic();
check=0;
fmax=20;
load(fullfile('TMP',sprintf('%d_params',imod)),'param');
if isfield(param,'cz_length')
    fmax  =fmax+param.cz_length;
end
if nargin<3
    ic=1;
else
    if isfield(param,'crack_id')
        ic=param.crack_id;
    else
        ic=1;
    end
end
if nargin<4, step=0;end
if nargin<5, dx=1;end
if nargin<6, dec2=0;end
load(fullfile('TMP','params'),'param');



fid=fullfile('TMP',sprintf('%d_levelsets',ic));
fid1=fullfile('TMP',sprintf('%d_levelsets_grads',ic));
fid2=fullfile('TMP',sprintf('%d_levelsets_cylco',ic));
load(fid,'crack')
load(fid,'front')
sizeim=size(crack);
if thp==0
    front=front-dec;
else
    
    if dx>1
        xo=1:dx:size(front,1);
        yo=1:dx:size(front,2);
        front=front(xo,yo);
        crack=crack(xo,yo);
    end
    
    if dec2>0
        
       front=front-dec2; 
    end
    
    dist=sqrt(sum(size(crack).^2));
    dist=1*dec;
    vcrack=dec*sin(thp)+0*crack;
    vfront=dec*cos(thp)+0*front;

    vcrack=LSExtend(vcrack,front,dist,dx);
    vcrack=LSExtend(vcrack,crack,dist,dx);

    vfront=LSExtend(vfront,front,dist,dx);
    vfront=LSExtend(vfront,crack,dist,dx);
    vcrack=LSAdjust(vcrack,vfront,front);

    crack=LSPropagate(vcrack,crack,dx);
    front=LSPropagate(vfront,front,dx);


    crack=LSReinit(crack,dist,dx);
    front=LSOrtho(front,crack,dist,dx);
    front=LSReinit(front,dist,dx);
if dx>1
    [Yi,Xi]=meshgrid(yo(1):yo(length(yo)),xo(1):xo(length(xo)));
   front=interp2(yo,xo,front,Yi,Xi,'*linear'); 
   crack=interp2(yo,xo,crack,Yi,Xi,'*linear'); 
    clear Xi Yi
    
end
    if dec2>0
        
       front=front+dec2; 
    end


end
if check
    figure
cb=colormap(hot);
cb=cb(:,[3,2,1]);
   colormap(cb);
    axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
        'XTick',zeros(1,0),...
        'FontSize',20);
    box('on');
    hold 'all';
    image(front','Parent',axes1,'CDataMapping','scaled');
    [c,h1] = contour(crack',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
    set(h1,'EdgeColor','black','LineWidth',2);
    [c,h2] = contour(front',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
    set(h2,'EdgeColor','black','LineWidth',2);
    axis off;
    axis xy;
    axis image;
    colorbar('peer',axes1,'FontSize',20);
    title('Front levelset','FontSize',24);
    if step
        print ('-djpeg', fullfile('FIG',[param.result_file,'-levelsets-',num2str(ic),'-front-',num2str(step)]));


    end
    figure
cb=colormap(hot);
cb=cb(:,[3,2,1]);
   colormap(cb);
    axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
        'XTick',zeros(1,0),...
        'FontSize',20);
    box('on');
    hold 'all';
    image(crack','Parent',axes1,'CDataMapping','scaled');
    [c,h1] = contour(crack',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
    set(h1,'EdgeColor','black','LineWidth',2);
    [c,h2] = contour(front',[0.,0.],'Parent',axes1,'CDataMapping','scaled');
    set(h2,'EdgeColor','black','LineWidth',2);
    axis off;
    axis xy;
    axis image;
    colorbar('peer',axes1,'FontSize',20);
    title('Crack levelset','FontSize',24);
    if step
        print ('-djpeg', fullfile('FIG',[param.result_file,'-levelsets-',num2str(ic),'-crack-',num2str(step)]));


    end


end



save(fid,'crack','front','sizeim','-append');

nx=FDgradient(crack,1);
save(fid1,'nx','-append');

ny=FDgradient(crack,2);
save(fid1,'ny','-append');
z=i*crack;
clear crack

nn=nx+i*ny;
clear nx ny

nnorm=abs(nn);
nn=nn./nnorm;
clear nnorm
nn=nn*exp(-i*pi/2);
thetap=angle(nn);
theta=mean(thetap(:));
clear nn
save(fid2,'theta','thetap','-append');



tx=FDgradient(front,1);
save(fid1,'tx','-append');
clear tx

ty=FDgradient(front,2);
save(fid1,'ty','-append');
clear ty
z=z+front;
clear front

dist=max(abs(z),1);
save(fid2,'dist','-append');
clear dist

angl=angle(z);
save(fid2,'angl','-append');
clear angl
onfaces=(imag(z)<=0)&(imag(z)>-1)&(real(z)<=fmax);
clear z
save(fid,'onfaces','-append');

disp(sprintf('Propagating levelsets...%6.2f s',toc()));

end

