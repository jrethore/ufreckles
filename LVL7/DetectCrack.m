function DetectCrack(fileres)
%%
printok=1;
filid='0';
nsca=10;
itermx=50;
filres='32Q4';
filres='16P1-5';
        load('TMP/sample0','sizeim');
        sizeim0=sizeim;
        load('TMP/sample0_0','sizeim');
load([filres,'.mat'])
roi=param.roi;

roi2(1)=xo(1)+roi(1)-1;roi2(2)=xo(length(xo))+roi(1)-1;
roi2(3)=yo(1)+roi(3)-1;roi2(4)=yo(length(yo))+roi(3)-1;
normal=1; % mean normal to the surface is x
%normal=2; % mean normal to the surface is y
%normal=3; % mean normal to the surface is z
fileres='TMP/1_error_0';
load(fileres);
disc=reshape(disc,sizeim);
res=abs(disc(xo(1):xo(length(xo))-1,yo(1):yo(length(yo))-1));
% figure
% imagesc(res)
clear disc;
% shift dimensions so that the mean normal is along z
if normal~=ndims(res)
    res=shiftdim(res,normal);
end

nx=size(res,1);ny=size(res,2);
nitermx=50;
%
%% Coarsen the residue field
fil=['TMP/res_' num2str(1) '.mat'];
save(fil,'res');
if nsca>1
    for isca=2:nsca
        res(:,2:ny-1)=(1/3)*(res(:,1:ny-2)+res(:,2:ny-1)+res(:,3:ny));
        fil=['TMP/res_' num2str(isca) '.mat'];
        save(fil,'res');
    end
end
%%
%
% Compute node locations for all grids
% (nodes are at integer voxel positions)
% mesh size is 2^n
% storage is such that x=[Xgrid(1,m):Xgrid(2,m);Xgrid(3,m)]
% and similarly for y.  We store in Xgrid(4,m) the number of NODES along x
%
close all
visco=250;
nmin=4
ngrid=-nmin+floor(log2(nx-1));
Xgrid=zeros(4,ngrid);
del=2.^(nmin+(0:ngrid-1));
Xgrid(2,:)=del;
ngx=ceil((nx-1)./del);
Xgrid(4,:)=ngx+1;
Xgrid(1,:)=1+floor((nx-1-ngx.*del)/2);
Xgrid(3,:)=Xgrid(1,:)+ngx.*del;

% Initialization of the surface as a straight cut
fil=['TMP/res_' num2str(1) '.mat'];
load(fil,'res');
sumres=sum(res);
[mi,imin]=max(sumres(:));
%restot=sum(res(:))-2*max(res(:))*(nx-1);
%restot=min(sumres(:))*length(sumres);
figure
plot(sumres(:));
sumres=sum(sumres(:));
Z0=imin
%%
%Z0=550;
Zmini=Z0-180;
Zmaxi=Z0+180;
Z=Z0*ones(nx-1,1);
% figure
% surf(Z([1,nx-1],[1,ny-1]));
%%
% loop over coarsening of the residuals AND grid refinement
ncoarse=max(ngrid-1,nsca);
oldisca=nsca+1;
Dini=0;critold=1;
jter=0;
clear LTObj LTcrit

for icoarse=ncoarse:-1:1
    isca=max(nsca+icoarse-ncoarse,1);
    if isca~=oldisca
        fil=['TMP/res_' num2str(isca) '.mat'];
        load(fil,'res');
    end
    res0=res;
    normres=(1./mean(res,2));
    normres=diag(sparse(normres));
    
    res=normres*res;
restot=sum(max(res,[],2));
    igrid=max(ngrid+icoarse-ncoarse,1);
    disp(sprintf('Coarsening %d : Res scale %d Grid scale %d',icoarse,isca,igrid));

    X=[Xgrid(1,igrid):Xgrid(2,igrid):Xgrid(3,igrid)];
    nX=Xgrid(4,igrid);
    ZZ=Z(min(max(1,X),nx-1)); % Get surface elevation at grid nodes
    %ZZ=Z(X,Y); % Get surface elevation at grid nodes
    xlin=((0:Xgrid(2,igrid)-1)+0.5)/Xgrid(2,igrid);
    shape1=(1-xlin);
    shape2=(xlin);
    visco=100;
    clear LObj Lcrit
    scrsz = get(0,'ScreenSize');
    figure1=figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)*0.75]) ;
    axes1 = axes('Parent',figure1,...
        'LineWidth',2,...
        'FontWeight','demi',...
        'FontSize',14);
    hold('all');
    if icoarse<ncoarse
    close(figure2);
    end
    for iter=1:nitermx
        F=zeros(nX,1);Obj=0;
        for iel=1:nX-1
            i1=iel;
            i2=i1+1;
            jx=[X(iel):X(iel+1)-1];
            z=ZZ(i1)*shape1+ZZ(i2)*shape2;
            %                 figure
            %                 imagesc(z);colorbar;
            %                 title(sprintf('Grid %d Element %d %d',icoarse,iel,jel));
            z0=floor(z-0.1);
            z1=ceil(z+.1);
            z2=floor(z);dz=z-z2;
            gz=repmat(0,Xgrid(2,igrid),1);
            for ijx=1:Xgrid(2,igrid)
                jjx=jx(ijx);
                if (jx(ijx)>0)&&(jx(ijx)<nx)

                    gz(ijx)=(res(jjx,z1(ijx))...
                        -res(jjx,z0(ijx)))...
                        /(z1(ijx)-z0(ijx));
                    Obj=Obj+res(jjx,z2(ijx))*...
                        (1-dz(ijx))+...
                        res(jjx,z2(ijx)+1)*...
                        (dz(ijx));
                end
            end
            F(i1)=F(i1)+(shape1*gz);
            F(i2)=F(i2)+(shape2*gz);
        end
        %       F(1,:)=0;F(nX,:)=0;F(:,1)=0;F(:,nY)=0;
        D=F;
        if Dini==0&&max(abs(D(:)))>0
            Dini=2*max(abs(D(:)))/1;
        end
        %         if max(abs(D(:)))>1
        %             D=D/max(abs(D(:)));
        %         end
        if Dini>0
%            D=D/Dini*(4*double(isca<2)+5)*visco;
            D=D/Dini*(1-isca/nsca)*visco;
        end
        if iter==1
            if norm(D(:))==0
                break
            end
            D0=norm(D(:));
            Dold=D0;
        end
        normD=norm(D(:));
        crit=abs(normD-Dold)/D0;
        climit=(max(isca,4)*1.e-4);
%        Obj=abs(1-(restot+Obj)/sumres);
        Obj=abs(1-Obj/restot);
        disp(sprintf('Iter %d Objective = %f',iter,Obj));
        disp(sprintf('Norm of move = %f',normD));
        if iter>1
            disp(sprintf('Convergence = %f',crit));
        end
        jter=jter+1;
        LTObj(jter)=Obj;
        LObj(iter)=Obj;
        if iter==1
            Lcrit(iter)=NaN;
        else
            Lcrit(iter)=crit;
        end
        % if jter==1
        % LTcrit(jter)=NaN;
        % else
        if iter>1
            LTcrit(iter,-icoarse+ncoarse+1)=crit;
        end
        %dcrit=(crit/critold);
        %         if iter>2&dcrit>0.5&(min(nX,nY)>3)
        %             visco=3;
        %         else
        % if iter==2&dcrit>0.9&(min(nX,nY)>3)
        %             visco=2;
        % %         else
        % %             visco=1;
        %         end

        % if iter>1&(min(nX,nY)>3)
        %    visco= 1/(isca*0.01*(norm(D(:)))/D0);
        % end
        ZZ=ZZ+D;
        Dold=normD;critold=crit;




        sub1=subplot(3,2,[3:6],'Parent',figure1,...
            'LineWidth',2,...
            'FontWeight','demi',...
            'FontSize',14);
        %       view([-37.5 30]);
        grid('on');
        hold('all');
        if iter>1
            delete(h);
        end
        h=image(res0);
        colormap(hot);
        h2=plot(ZZ,X,'w-','LineWidth',5);
        %       axis([Y(1) Y(nY) X(1) X(nX) Zmini Zmaxi]);
        title(sprintf(' Residual scale %d Grid step %d Iter %d ',isca,Xgrid(2,igrid),iter),...
            'FontSize',16,'FontWeight','demi');
        subplot(3,2,1,'Parent',figure1,...
            'LineWidth',2,...
            'FontWeight','demi',...
            'FontSize',14);
        hold('all');

        plot(LObj,'LineWidth',2,'Color',[0 0 0])
        title('Relative distance to the objective function','FontWeight','demi','FontSize',16)
        xlabel('Iteration','FontWeight','demi','FontSize',14);
        subplot(3,2,2,'Parent',figure1,'YScale','log','YMinorTick','on',...
            'LineWidth',2,...
            'FontWeight','demi',...
            'FontSize',14);
        hold('all');
        semilogy([Lcrit'],'LineWidth',2,'Color',[0 0 0])
        semilogy([0*Lcrit'+climit],'LineStyle','--','Color',[0 0 0])
        ylim([1e-4,1e0])
        title('Convergence criterion','FontWeight','demi','FontSize',16)
        xlabel('Iteration','FontWeight','demi','FontSize',14);



        pause(0.1)
        %         sortiefilename =['FIG/detect-scale-',filid,'-',num2str(icoarse),'-iter-',num2str(iter)];
        %          if printok
        %              print ('-djpeg', sortiefilename);
        %          end
        if iter>1&crit<climit
            break
        end
    end % on iter
figure2=figure1; 
    figure22=figure ;
    axes1 = axes('Parent',figure22,...
        'LineWidth',2,...
        'FontWeight','demi',...
        'FontSize',14);
    hold('all');

    h=image(res0);
    colormap(hot);
    h2=plot(ZZ,X,'w-','LineWidth',5);
    title(sprintf(' Residual scale %d Grid step %d',isca,Xgrid(2,igrid)),...
        'FontSize',16,'FontWeight','demi');
    sortiefilename =['FIG/detect-',filid,'-scale-',num2str(icoarse),'-converged'];
    if printok
        print ('-djpeg', sortiefilename);
    end
    for iel=1:nX-1
        i1=iel;i2=i1+1;
        jx=[X(iel):X(iel+1)-1];
        z=ZZ(i1)*shape1+ZZ(i2)*shape2;
        %Z(jx,jy)=z;
        for ijx=1:Xgrid(2,igrid)
            jjx=jx(ijx);
            if (jx(ijx)>0)&&(jx(ijx)<nx)
                Z(jx(ijx))=z(ijx);
            end
        end
    end
end
figure
plot(LTObj)
title('Relative distance to the objective function')
xlabel('Iteration');
sortiefilename =['FIG/detect-',filid,'-objective'];
if printok
    print ('-djpeg', sortiefilename);
end
figure
semilogy([LTcrit])
title('Convergence criterion')
xlabel('Iteration');
sortiefilename =['FIG/detect-',filid,'-convergence'];
if printok
    print ('-djpeg', sortiefilename);
end
LTcrit(find(isnan(LTcrit)))=-1;
resm=double([[1:length(LTObj)]' LTObj' ]);
eval(['save ',['detect-obj.dat'],' resm',' -ASCII']);
resm=double([[1:size(LTcrit,1)]' LTcrit ]);
eval(['save ',['detect-conv.dat'],' resm',' -ASCII']);



%% reconstruction

idim=circshift([1:2]',-normal);


x1=X+roi2(2*idim(1)-1);

ztmp=ZZ+roi2(2*idim(2)-1)-1;

%%
x11=1:sizeim0(idim(1));

zztmp=interp1(x1,ztmp,x11,'linear','extrap');
figure
hold on
plot(ztmp,x1,'bs')
plot(zztmp,x11,'b-')

%% crude level set
[Yp Xp]=meshgrid(1:sizeim0(idim(2)),1:sizeim0(idim(1)));
front=-Xp+max(Xp(:));
crack=(Yp-repmat(zztmp',[1 sizeim0(normal)]));
    figure10=figure;
    axes1 = axes('Parent',figure10,...
        'LineWidth',2,...
        'FontWeight','demi',...
        'FontSize',14);
    hold('all');
        imagesc(crack);
%        axis image
        colormap(hot);
        colorbar;
        title(sprintf('Crack Level Set'),...
            'FontSize',16,'FontWeight','demi');
     figure20=figure;
    axes1 = axes('Parent',figure20,...
        'LineWidth',2,...
        'FontWeight','demi',...
        'FontSize',14);
    hold('all');
        imagesc(front);
 %       axis image
        colormap(hot);
        colorbar;
        title(sprintf('Front Level Set'),...
            'FontSize',16,'FontWeight','demi');

chsign=menu('Change sign ?','Yes','No');
if chsign==1
crack=-crack;
front=-front;
front=front-max(front(:));
end


if normal~=ndims(res)
    front=shiftdim(front,ndims(res)-normal);
    crack=shiftdim(crack,ndims(res)-normal);
end
figure
imagesc(crack)
%crack=LSReinit(crack,100,1);
sauve=menu('Save this crack ?','Yes','No');
if sauve==1
    [fil2,path1]=uiputfile('*.mat','Select geometry file');
    save(fil2,'crack','front');

end


