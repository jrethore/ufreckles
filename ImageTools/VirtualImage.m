function VirtualImage(nmod)
if nargin<1, nmod=1;end
load(fullfile('TMP','params'),'param');
tic;
given=[];
check=1;
filref=param.reference_image;
roi=param.roi;
reverse=0;
if isfield(param,'reverse_image')
    reverse=param.reverse_image;
end
tau=param.transition_length;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
thickness=0;
if isfield(param,'line_thickness')
    thickness=round(0.5*(param.line_thickness));
end
nscale=param.nscale;
mim=0;
if isfield(param,'min_grey_level')
    mim=param.min_grey_level;
end
dim=255;
if isfield(param,'max_grey_level')
    dim=(param.max_grey_level)-mim;
end
if length(roi)==4
    [patho, filo, exto] = fileparts(filref);
    if strcmp(exto,'.mat')
        load(filref);
        if reverse
            lso=lso';
            ls1=ls1';
        end
        if exist('nx')
            if reverse
                nx=nx';
                ny=ny';
            end
        end
        if exist('im0')
            given=im0;
        end
        sizeim=size(lso);
    else
        ls1=[];
        lso=double(readim(filref));
        if length(sizeim)==3
            lso=mean(lso,3);
            sizeim=size(lso);
        end
        if reverse
            lso=lso';
        end
        sizeim=size(lso);
        lso=lso-min(lso(:));
        lso=lso/max(lso(:));
        lso=2*lso-1;
        lso=LSReinit(lso,tau*2^(nscale-1),1);
    end
    if isempty(given)
    switch param.contour_type
        case 'edge'
            im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
        case 'line'
            im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
    end
    else
        im0=given;
    end
    if check
        figure
        axes1 = axes('YTick',zeros(1,0),'YDir','reverse',...
            'XTick',zeros(1,0),...
            'FontSize',20,'CLim',[0,255]);
        box('on');
        image(im0','Parent',axes1,'CDataMapping','scaled');
        hold 'all';
        rect1=rectangle;
        set(rect1,'position',[roi(1),roi(3),roi(2)-roi(1),roi(4)-roi(3)],...
            'EdgeColor','yellow',...
            'LineStyle','--',...
            'LineWidth',2)
        axis xy;
        axis image;
        colorbar('peer',axes1,'FontSize',20);
        colormap('gray')
        
        title('Virtual Image','FontSize',24);
    end
    save(fullfile('TMP','sample0'),'im0','lso','ls1','sizeim','-v7.3');
    lso=(lso(roi(1):roi(2),roi(3):roi(4)));
    if ~isempty(ls1)
        ls1=(ls1(roi(1):roi(2),roi(3):roi(4)));
    end
    if isempty(given)
    switch param.contour_type
        case 'edge'
            im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
            nband=find((lso(:)>=0)&(lso(:)<=tau));
        case 'line'
            im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
            nband=find(abs(lso(:))<=(tau+thickness));
    end
    else
        
       im0=(given(roi(1):roi(2),roi(3):roi(4)));
    switch param.contour_type
        case 'edge'
            nband=find((lso(:)>=0)&(lso(:)<=tau));
        case 'line'
            nband=find(abs(lso(:))<=(tau+thickness));
    end
       
       
       
    end
    if exist('nx')
        nx=(nx(roi(1):roi(2),roi(3):roi(4)));
        ny=(ny(roi(1):roi(2),roi(3):roi(4)));
    else
        nx=FDgradient(lso,1);
        ny=FDgradient(lso,2);
    end
    nnorm=abs(nx+i*ny);
    nnorm=(nnorm<=eps)+nnorm;
    nx=nx./nnorm;
    ny=ny./nnorm;
    sizeim=size(im0);
    on=find((lso(nband)<1)&(lso(nband)>=0));
%      ls1on=ls1(nband(on));
%      keep=(ls1(nband)>=min(ls1on))&(ls1(nband)<=max(ls1on));
%      nband=nband(on);
%      on=find((lso(nband)<1)&(lso(nband)>=0));
    
    mean0=mean(im0(nband));
    std0=std(im0(nband));
    save(fullfile('TMP','sample0_0'),'im0','nx','ny','lso','ls1','sizeim','nband','on','mean0','std0','-v7.3');
else
    [patho, filo, exto] = fileparts(filref);
    assert(strcmp(exto,'.mat'))
    load(filref,'lso');
    sizeim=size(lso);
    
    %     switch param.contour_type
    %         case 'edge'
    %             im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
    %         case 'line'
    %             im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
    %     end
    %     save(fullfile('TMP','sample0'),'im0','lso','ls1','ls2','sizeim');
    lso=(lso(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)));
    switch param.contour_type
        case 'edge'
            im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
            nband=find((lso(:)>=0)&(lso(:)<=tau));
        case 'line'
            im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
            nband=find(abs(lso(:))<=(tau+thickness));
    end
    sizeim=size(im0);
    on=find((lso(nband)<1)&(lso(nband)>=0));
    %    save(fullfile('TMP','sample0_0'),'im0','nx','ny','nz','lso','ls1','ls2','sizeim','on','nband');
    mean0=mean(im0(nband));
    std0=std(im0(nband));
    save(fullfile('TMP','sample0_0'),'im0','lso','sizeim','on','nband','mean0','std0','-v7.3');
    nx=FDgradient(lso,1);
    nx=nx(nband);
    save(fullfile('TMP','sample0_0'),'nx','-append');
    clear nx
    ny=FDgradient(lso,2);
    ny=ny(nband);
    save(fullfile('TMP','sample0_0'),'ny','-append');
    clear ny
    nz=FDgradient(lso,3);
    nz=nz(nband);
    save(fullfile('TMP','sample0_0'),'nz','-append');
    clear nz lso
    
    %     nnorm=sqrt(nx.^2+ny.^2+nz.^2);
    %     nnorm=(nnorm<=eps)+nnorm;
    %     nx=nx./nnorm;
    %     ny=ny./nnorm;
    %     nz=nz./nnorm;
    load(filref,'ls1');
    ls1=(ls1(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)));
    ls1=ls1(nband);
    save(fullfile('TMP','sample0_0'),'ls1','-append');
    clear ls1
    load(filref,'ls2');
    ls2=(ls2(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)));
    ls2=ls2(nband);
    save(fullfile('TMP','sample0_0'),'ls2','-append');
    clear ls2
end

disp(sprintf('Loading reference image...%6.2f s',toc));
for iscale=2:nscale
    tic
    tau=2*tau;
    switch param.contour_type
        case 'edge'
            im0=mim+dim*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
            nband=find((lso(:)>=0)&(lso(:)<=tau));
        case 'line'
            im0=mim+dim*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
            nband=find(abs(lso(:))<=(tau+thickness));
    end
    %      figure
    %      plot(im0(:,1))
    on=find((lso(nband)<1)&(lso(nband)>=0));
    sizeim=size(im0);
    mean0=mean(im0(nband));
    std0=std(im0(nband));
    save(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim','nband','on','mean0','std0','-v7.3');
    disp(sprintf('   Coarsening level %d...%6.2f s',iscale,toc()));
end


% for iscale=2:nscale
%     tic();
%     NestedCoarseImage();
% im0=mim+dim*0.5*(1-cos(pi*min(ls,tau)/tau)).*double(ls>=0);
%     sizeim=size(im0);
% nx=FDgradient(ls,1);
% ny=FDgradient(ls,2);
% nnorm=abs(nx+i*ny);
% nnorm=(nnorm<=eps)+nnorm;
% nx=nx./nnorm;
% ny=ny./nnorm;
%     save(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','nx','ny','sizeim');
%     disp(sprintf('   Coarsening level %d...%6.2f s',iscale,toc()));
% end

    function NestedCoarseImage()
        
        scale=2;
        imsiz0=size(lso);
        imsiz1=floor(imsiz0/2);
        nn=2*imsiz1;
        ndim=length(nn);
        if ndim==2
            lso=lso(1:nn(1),1:nn(2));
        elseif ndim==3
            lso=lso(1:nn(1),1:nn(2),1:nn(3));
        end
        
        lso=reshape(lso,scale,prod(nn)/scale);
        lso=mean(lso,1);
        nn(1)=nn(1)/scale;
        lso=reshape(lso,nn);
        
        if ndim==2
            lso=lso';
            lso=reshape(lso,scale,prod(nn)/scale);
            lso=mean(lso,1);
            nn(2)=nn(2)/scale;
            lso=reshape(lso,nn([2,1]));
            lso=lso';
            %toc();
        elseif ndim==3
            lso=permute(lso,[2,3,1]);
            lso=reshape(lso,scale,prod(nn)/scale);
            lso=mean(lso,1);
            nn(2)=nn(2)/scale;
            lso=reshape(lso,nn([2,3,1]));
            lso=permute(lso,[3,1,2]);
            
            lso=permute(lso,[3,1,2]);
            lso=reshape(lso,scale,prod(nn)/scale);
            lso=mean(lso,1);
            nn(3)=nn(3)/scale;
            lso=reshape(lso,nn([3,1,2]));
            lso=permute(lso,[2,3,1]);
        end
        
    end




end

