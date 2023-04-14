function GenUnitCellMesh(n,ll,ep,phio)
if nargin<4,phio=0;end
%ll=100;
%ep=20;
check=1;
r=ep/3;
h=ep/4;
hc=0.25*h;
w=2*pi/n;
if n>2
    anglep=0;
    wp=linspace(0,w,3)-w/2;
    ep=ep/cos(w/2);
    ll=ll/2/sin(w/2)-ep;
    zp=ll*exp(1i*wp');
    zp=[zp;conj(zp)*(ll+ep)/ll];
    zp=[zp;zp(1)];
    zp=zp.*cos(angle(zp));
    mesh.zone{1,1}=1;
    mesh.zone{2,1}=[real(zp),imag(zp)];
    mesh.zone{3,1}=0;
    mesh.zone{4,1}=2;
    mesh.zone{5,1}=0;
    mesh.zone{6,1}=0;
    
    mesh.zone{1,2}=-2;
    mesh.zone{2,2}=[ll,0];
    mesh.zone{3,2}=0;
    mesh.zone{4,2}=3;
    mesh.zone{5,2}=[ll,0,r];
    mesh.zone{6,2}=hc;
    mesh.mesh_size=h/2;
    roi=[-ep,ll+2*ep,-ll-2*ep,ll+2*ep];
    ld=@(xy) GetSignedDistanceToZone(mesh,roi,xy(:,1),xy(:,2));
    lh=@(xy) (GetMeshDensity(mesh,xy(:,1),xy(:,2)));
    [d,~,indc]=GetMeshDensity(mesh,real(zp),imag(zp));
    zsup=linspace(ll,ll+ep,ep/h+2)'*exp(1i*max(wp))*cos(max(wp));
    zinf=linspace(ll,ll+ep,ep/h+2)'*exp(1i*min(wp))*cos(min(wp));
    zcircle=exp(1i*linspace(-pi/2-w/2,pi/2+w/2,round(r*abs(2*pi-w)/hc)+2)')*r+ll;
    zcorner=ll;
    figure
    plot(zp,'k')
    axis equal
    hold on
    plot(zinf,'ro')
    plot(zsup,'ro')
    plot(zcircle,'ro')
    plot(zcorner,'ro')
    zfix=[zcircle;zsup;zinf;zcorner];
    [xo,yo,conn,xyfix]=GenMeshFromLVL7(roi,h*min(d),ld,lh,[real(zfix),imag(zfix)],0,10000);
    triplot(conn,xo,yo)
    plot(xyfix(:,1),xyfix(:,2),'sm')
    indc=(1:numel(zcircle))';
    pause
    %xc=mean(xo(conn),2);yc=mean(yo(conn),2);
    %plot(xc,yc,'xg')
    zo=xo+1i*yo;
    conno=conn;
    indco=indc;
    for ii=1:n-1
        xo=[xo;real(zo*exp(1i*ii*w))];
        yo=[yo;imag(zo*exp(1i*ii*w))];
        indc=[indc;indco+ii*numel(zo)];
        conn=[conn;conno+ii*numel(zo)];
        anglep=[anglep;anglep(end)+w];
    end
    zp=(ll+ep)*exp(1i*linspace(0,2*pi-w,n));
    nn=[(n/2-1):(n/2-1):(n-1)];
    homo=mod(repmat((1:n)',1,numel(nn))-1+repmat(nn,n,1),n)+1;
    %     figure
    %     plot(zp,'k')
    %     axis equal
    %     hold on
    % for in=1:n
    %     hn=plot(zp(homo(in,:)),'ro');
    %     pause
    %     delete(hn)
    % end
    if n==6
        %     homo=[homo,mod(repmat((1:n)',1,1)-1+1,n)+1,mod(repmat((1:n)',1,1)-1-1,n)+1];
        
    end
    zpp=zp;
else
    zp=[ll;ll+ep;ll+ep+1i*(ll+2*ep);ep+1i*(ll+2*ep);ep+1i*(2*ll+2*ep);1i*(2*ll+2*ep);1i*(ll);ll+1i*(ll)];
    zp=[zp;zp(1)];
    mesh.zone{1,1}=1;
    mesh.zone{2,1}=[real(zp),imag(zp)];
    mesh.zone{3,1}=0;
    mesh.zone{4,1}=2;
    mesh.zone{5,1}=0;
    mesh.zone{6,1}=0;
    
    mesh.zone{1,end+1}=-2;
    mesh.zone{2,end}=[ep,2*ep+ll];
    mesh.zone{3,end}=0;
    mesh.zone{4,end}=3;
    mesh.zone{5,end}=[ep,2*ep+ll,r];
    mesh.zone{6,end}=hc;
    
    mesh.zone{1,end+1}=-2;
    mesh.zone{2,end}=[ll,ll];
    mesh.zone{3,end}=0;
    mesh.zone{4,end}=3;
    mesh.zone{5,end}=[ll,ll,r];
    mesh.zone{6,end}=hc;
    
    % zp=[ll+ep;ll+ep+1i*2*(ll+ep);-(ll+ep)+1i*2*(ll+ep);-(ll+ep)];
    % mesh.zone{1,1}=1;
    % mesh.zone{2,1}=[real(zp),imag(zp)];
    % mesh.zone{3,1}=0;
    % mesh.zone{4,1}=2;
    % mesh.zone{5,1}=0;
    % mesh.zone{6,1}=0;
    %
    % zp=[ep+1i*(ll+2*ep);ll+ep+1i*(ll+2*ep);ll+ep+1i*2*(ll+ep);ep+1i*2*(ll+ep)];
    % mesh.zone{1,end+1}=0;
    % mesh.zone{2,end}=[real(zp),imag(zp)];
    % mesh.zone{3,end}=0;
    % mesh.zone{4,end}=2;
    % mesh.zone{5,end}=0;
    % mesh.zone{6,end}=0;
    %
    % zp=[-ep+1i*(ll+2*ep);-ll-ep+1i*(ll+2*ep);-ll-ep+1i*2*(ll+ep);-ep+1i*2*(ll+ep)];
    % mesh.zone{1,end+1}=0;
    % mesh.zone{2,end}=[real(zp),imag(zp)];
    % mesh.zone{3,end}=0;
    % mesh.zone{4,end}=2;
    % mesh.zone{5,end}=0;
    % mesh.zone{6,end}=0;
    %
    % zp=[ll;ll+1i*ll;-ll+1i*ll;-ll];
    % mesh.zone{1,end+1}=0;
    % mesh.zone{2,end}=[real(zp),imag(zp)];
    % mesh.zone{3,end}=0;
    % mesh.zone{4,end}=2;
    % mesh.zone{5,end}=0;
    % mesh.zone{6,end}=0;
    
    % mesh.zone{1,end+1}=-2;
    % mesh.zone{2,end}=[ll,ll];
    % mesh.zone{3,end}=0;
    % mesh.zone{4,end}=3;
    % mesh.zone{5,end}=[ll,ll,r];
    % mesh.zone{6,end}=hc;
    %
    % mesh.zone{1,end+1}=-2;
    % mesh.zone{2,end}=[-ll,ll];
    % mesh.zone{3,end}=0;
    % mesh.zone{4,end}=3;
    % mesh.zone{5,end}=[-ll,ll,r];
    % mesh.zone{6,end}=hc;
    
    % mesh.zone{1,end+1}=-2;
    % mesh.zone{2,end}=[ep,ll+2*ep];
    % mesh.zone{3,end}=0;
    % mesh.zone{4,end}=3;
    % mesh.zone{5,end}=[ep,ll+2*ep,r];
    % mesh.zone{6,end}=hc;
    %
    % mesh.zone{1,end+1}=-2;
    % mesh.zone{2,end}=[-ep,ll+2*ep];
    % mesh.zone{3,end}=0;
    % mesh.zone{4,end}=3;
    % mesh.zone{5,end}=[-ep,ll+2*ep,r];
    % mesh.zone{6,end}=hc;
    
    mesh.mesh_size=h/2;
    
    roi=[-ep,ll+2*ep,-ep,2*ll+3*ep];
    ld=@(xy) GetSignedDistanceToZone(mesh,roi,xy(:,1),xy(:,2));
    lh=@(xy) (GetMeshDensity(mesh,xy(:,1),xy(:,2)));
    [d,~,indc]=GetMeshDensity(mesh,real(zp),imag(zp));
    zcircle=exp(1i*linspace(pi/2,2*pi,round(r*abs(3*pi/2)/hc)+2)')*r+(2*ep+ll)*1i+ep;
    zcirclei=exp(1i*linspace(-pi/2,pi,round(r*abs(3*pi/2)/hc)+2)')*r+(ll)*1i+ll;
    zcorner=(2*ep+ll)*1i+ep;
    zcorneri=(ll)*1i+ll;
    
    figure
    plot(zp,'k')
    axis equal
    hold on
    plot(zcircle,'ro')
    plot(zcorner,'ro')
    plot(zcirclei,'ro')
    plot(zcorneri,'ro')
    
    zfix=[zcircle;zcirclei;zcorner;zcorneri];
    [xo,yo,conn,xyfix]=GenMeshFromLVL7(roi,h*min(d),ld,lh,[real(zfix),imag(zfix)],0,10000);
    triplot(conn,xo,yo)
    plot(xyfix(1:numel(zcircle),1),xyfix(1:numel(zcircle),2),'sm')
    plot(xyfix((numel(zcircle)+1):numel(zcirclei),1),xyfix((numel(zcircle)+1):numel(zcirclei),2),'sr')
    
    indc=(1:(numel(zcircle)+numel(zcirclei)))';
    conno=conn;
    indco=indc;
    nn=numel(xo);
    xo=[xo;-xo];
    yo=[yo;yo];
    indc=[indc;indco+nn];
    conn=[conn;conno+nn];
    nn=numel(xo);
    xo=[xo;xo];
    yo=[yo;-yo];
    indc=[indc;indc+nn];
    conn=[conn;conn+nn];
    n=4;
    zp=[-(ll+ep)-1i*2*(ll+ep);(ll+ep)-1i*2*(ll+ep);(ll+ep)+1i*2*(ll+ep);-(ll+ep)+1i*2*(ll+ep)];
    
    zpp=[ep+1i*(2*ep+ll);-ep+1i*(2*ep+ll);ep-1i*(2*ep+ll);-ep-1i*(2*ep+ll)];
    zpp=[zpp;ll*(-1-1i);ll*(1-1i);ll*(-1+1i);ll*(1+1i)];
    anglep=[-3*pi/4;-pi/4;3*pi/4;pi/4;-3*pi/4;-pi/4;3*pi/4;pi/4];
    homo=[(5:8)';(1:4)'];
    resort=[1,8,2,7,3,6,4,5];
    zpp=zpp(resort);
    anglep=anglep(resort);
    homo=[8;7;6;5;4;3;2;1];
end

%%
d=0.1;
[~,ind,jnd]=unique(round([xo,yo]/d)*d,'rows','stable');
indc=jnd(indc);
conn=reshape(jnd(conn),size(conn));
xo=xo(ind);
yo=yo(ind);
indc=reshape(indc,numel(indc)/size(homo,1),size(homo,1));
zo=xo+1i*yo;
zo=zo*exp(1i*phio);
zp=zp*exp(1i*phio);
zpp=zpp*exp(1i*phio);
anglep=anglep+phio;
xo=real(zo);
yo=imag(zo);
if check
    for ii=1:numel(zpp)
        hh=figure;
        triplot(conn,xo,yo)
        hold on
        axis equal
        hn=plot(xo(indc(:,ii)),yo(indc(:,ii)),'ro');
        % hp=plot(real(zp(ii)),imag(zp(ii)),'ms','MarkerSize',20,'LineWidth',2);
        hpp=plot(real(zpp(ii)),imag(zpp(ii)),'kv','MarkerSize',20,'LineWidth',2);
        for jj=1:size(homo,2)
            plot([zpp(ii),zpp(homo(ii,jj))],'-k','LineWidth',2);
        end
        pause(1)
        delete(hh)
    end
end
hh=figure;
triplot(conn,xo,yo)
hold on
axis equal
plot(xo(indc),yo(indc),'ro')
plot(zp,'ms')
plot(zpp,'kv')
elt=3*ones(size(conn,1),1);
save('cell-mesh','xo','yo','conn','elt','n','zp','zpp','indc','homo','anglep')
writeVTKmesh('cell-mesh.mat')
end
