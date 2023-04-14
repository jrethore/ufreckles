function [xo,yo,conn,xyfixn]=GenMeshFromLVL7(roi,h,ld,lh,xyfixo,check,maxiter,dutol)
if nargin<6,check=0;end
if nargin<5,xyfixo=[];end
if nargin<7,maxiter=1000;end
if nargin<8,dutol=.001;end
utol=.1;
aa=1.2;
alpha=.2;
deps=.001*h;
dxy=sqrt(eps)*h;
meshok=0;ntry=0;
zxy=[1;1i];
tiers=[1;1;1]/3;
xyfixn=xyfixo;
if ~isempty(xyfixo)
    out=feval(ld,xyfixo)>=deps;
%    out=feval(ld,xyfixo)>0.75*h*feval(lh,xyfixo);
    if any(out)
        xyfixout=xyfixo(out,:);
        xyfixin=xyfixo(~out,:);
        in=find(~out);
        out=find(out);
        for ip=1:size(xyfixin,1)
            [dmin,id]=min(abs(xyfixout*zxy-xyfixin(ip,:)*zxy));
            ldn=feval(ld,xyfixin(ip,:));
            ldo=feval(ld,xyfixout(id,:));
            [dmin]=min(abs(xyfixin*zxy-xyfixout(id,:)*zxy));
            if (abs(ldo)<dmin)&&(abs(ldn)<dmin)
%                keyboard
                if (abs(ldn)<0.25*dmin)
                 dldnx=(feval(ld,[xyfixin(ip,1)+dxy,xyfixin(ip,2)])-ldn)/dxy;
                 dldny=(feval(ld,[xyfixin(ip,1),xyfixin(ip,2)+dxy])-ldn)/dxy;
                ndld=dldnx.^2+dldny.^2;
                dxyo=ldn/ndld*[dldnx,dldny];
                    
                                        xyfixo(in(ip),:)=xyfixo(in(ip),:)-dxyo;
                    
                else
                dldnx=(feval(ld,[xyfixout(id,1)+dxy,xyfixout(id,2)])-ldo)/dxy;
                dldny=(feval(ld,[xyfixout(id,1),xyfixout(id,2)+dxy])-ldo)/dxy;
                ndld=dldnx.^2+dldny.^2;
                dxyo=ldo/ndld*[dldnx,dldny];
                                        xyfixo(out(id),:)=xyfixo(out(id),:)-dxyo;

                    
                    
                end
%                 dldnx=(feval(ld,[xyfixo(ip,1)+dxy,xyfixo(ip,2)])-ldn)/dxy;
%                 dldny=(feval(ld,[xyfixo(ip,1),xyfixo(ip,2)+dxy])-ldn)/dxy;
%                 dldnx=(feval(ld,[xyfixout(id,1)+dxy,xyfixout(id,2)])-ldo)/dxy;
%                 dldny=(feval(ld,[xyfixout(id,1),xyfixout(id,2)+dxy])-ldo)/dxy;
%                 ndld=dldnx.^2+dldny.^2;
%                 dxyo=ldn/ndld*[dldnx,dldny];
%                 hh=feval(lh,xyfixo(ip,:))*h;
%                 np=4;
%                 try
%                     for ii=1:np-1
%                     if abs(xyfixo(ip,:)*zxy-xyfixo(ip+ii,:)*zxy)<(ii+1)*hh
%                         xyfixo(ip+ii,:)=xyfixo(ip+ii,:)-(np-ii)*dxyo/np;
%                     end
%                     end
%                 catch
%                 end
%                 try
%                     for ii=1:np-1
%                     if abs(xyfixo(ip,:)*zxy-xyfixo(ip-ii,:)*zxy)<(ii+1)*hh
%                         xyfixo(ip-ii,:)=xyfixo(ip-ii,:)-(np-ii)*dxyo/np;
%                     end
%                     end
%                 catch
%                 end
            end
            
        end
            out=feval(ld,xyfixo)>=deps;
            xyfixn=xyfixo;
        xyfixo=xyfixo(~out,:);
%            xyfixn=xyfixo;
    end
end
while ~meshok&&ntry<1
yo=roi(3):h*sqrt(3)/2:roi(4);
yo=(yo-roi(3)-1)*(roi(4)-1-roi(3)-1)/(yo(end)-roi(3)-1)+roi(3)+1;
xo=roi(1):h:roi(2);
xo=(xo-roi(1)-1)*(roi(2)-1-roi(1)-1)/(xo(end)-roi(1)-1)+roi(1)+1;
[yo,xo]=meshgrid(yo,xo);
xo(:,2:2:end)=xo(:,2:2:end)+h/2;
xo(:,2)=roi(1);
xo(:,end)=roi(2);
xyo=[xo(:),yo(:)];
xyo=xyo(feval(ld,xyo)<deps,:);
if 1
p=0.5*feval(lh,xyo).^(-2);
%p=exp(-feval(lh,xyo).^(2));
xyo=xyo(rand(size(xyo,1),1)<p/max(p),:);
%conn=delaunay(xyo);
% figure
%         triplot(conn,xyo(:,1),xyo(:,2))
%         axis equal
%         xlim(roi(1:2))
%         ylim(roi(3:4))
%         pause(0.2)

else
figure
    conn=delaunay(xyo);
        xyc=[reshape(xyo(conn,1),size(conn))*tiers,reshape(xyo(conn,2),size(conn))*tiers];
        refine=1
        while refine
        triplot(conn,xyo(:,1),xyo(:,2))
        axis equal
        xlim(roi(1:2))
        ylim(roi(3:4))
        pause
        keep=feval(lh,xyc)<1/(3^(refine-1));
        if any(keep)
       xyo=[xyo;xyc(keep,:)];
       refine=refine+1;
    conn=delaunay(xyo);
        xyc=[reshape(xyo(conn,1),size(conn))*tiers,reshape(xyo(conn,2),size(conn))*tiers];
        else
            refine=0;
        end
        end
       
end
[in,on]=inpolygon(xyo(:,1),xyo(:,2),roi([1,2,2,1,1]),roi([3,3,4,4,3]));
%on=[];
xyfix=[xyfixo;xyo(on,:)];
%[~,ind,~]=unique(round(2*xyfix/h)*h/2,'rows','stable');
%[~,ind,~]=unique(round(2*xyfix/h)*h/2,'rows','first','legacy');
%xyfix=xyfix(ind,:);
xyfix=xyfix(feval(ld,xyfix)<deps,:);
if ~isempty(xyfix)
    xyo=setdiff(xyo,xyfix,'rows');
end
%xyfix=unique(xyfix,'rows');
nc=size(xyfix,1);
xyo=[xyfix;xyo];
Nnodes=size(xyo,1);

ii=0;
xyp=Inf;
res=Inf;
if check
    figure
end
while (res>(dutol*h)&&ii<maxiter)||ii<10
    ii=ii+1;
    if max(abs((xyo-xyp)*zxy))>utol*h||isinf(xyp(1))
        xyp=xyo;
        try conn=delaunay(xyo);catch, end
        xyc=[reshape(xyo(conn,1),size(conn))*tiers,reshape(xyo(conn,2),size(conn))*tiers];
        keep=feval(ld,xyc)<-deps;
        conn=conn(keep,:);
        seg=reshape(conn(:,[1,2,1,3,2,3])',2,3*size(conn,1))';
        seg=unique(sort(seg,2),'rows');
        Nseg=size(seg,1);
    end
    if check&&(ii/50==round(ii/50))
        triplot(conn,xyo(:,1),xyo(:,2))
        title(sprintf('iteration %d',ii'))
        axis equal
        xlim(roi(1:2))
        ylim(roi(3:4))
        pause(0.2)
    end
    xys=xyo(seg(:,1),:)-xyo(seg(:,2),:);
    ls=abs(xys*zxy);
    hs=feval(lh,xyo(seg(:,2),:)+0.5*xys)*h;
    lp=hs*(aa*norm(ls)/norm(hs));
    
    if any(lp>1.5*ls)&&mod(ii,maxiter)==0
        reject=setdiff(reshape(seg(lp>2*ls,:),[],1),1:nc);
        xyo(reject,:)=[];
        Nnodes=size(xyo,1);
        xyp=Inf;
        continue;
    end
    
    dL=((max(lp-ls,0)./ls)*[1,1]).*xys;
    n1=sparse(seg(:,1),1:Nseg,1,Nnodes,Nseg);
    n2=sparse(seg(:,2),1:Nseg,-1,Nnodes,Nseg);
    U=(n1+n2)*dL;
    U(1:nc,:)=0;
%     U(xyo(:,1)==roi(1),:)=0;
%     U(xyo(:,2)==roi(3),:)=0;
%     U(xyo(:,1)==roi(2),:)=0;
%     U(xyo(:,2)==roi(4),:)=0;
    xyo=xyo+alpha*U;
    
    ldn=feval(ld,xyo);
    ind=ldn>0;
    if any(ind)
    dldnx=(feval(ld,[xyo(ind,1)+dxy,xyo(ind,2)])-ldn(ind))/dxy;
    dldny=(feval(ld,[xyo(ind,1),xyo(ind,2)+dxy])-ldn(ind))/dxy;
    ndld=dldnx.^2+dldny.^2;
    dxyo=diag(sparse(ldn(ind)./ndld))*[dldnx,dldny];
    xyo(ind,:)=xyo(ind,:)-dxyo;
    end
   U=U(ldn<-deps,:);
    res=alpha*max(abs(U*zxy))/h;
    
end
if res>(dutol*h)
%display('Convergence failure in mesh optimization')
else
%display(sprintf('Convergence achieved in mesh optimization after %d iterations',ii))
meshok=1;
end
ntry=ntry+1;
end
% figure
%         triplot(conn,xyo(:,1),xyo(:,2))
%         axis equal
if max(conn(:))>size(xyo,1)
conn=delaunay(xyo);
        xyc=[reshape(xyo(conn,1),size(conn))*tiers,reshape(xyo(conn,2),size(conn))*tiers];
        keep=feval(ld,xyc)<-deps;
        conn=conn(keep,:);
end
xo=xyo(:,1);yo=xyo(:,2);
[xo,yo,conn]=CleanTriMesh(xyo(:,1),xyo(:,2),conn);

end