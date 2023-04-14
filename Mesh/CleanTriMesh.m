function [xo,yo,conn,zo]=CleanTriMesh(xo,yo,conn,zo)
if nargin<4
    ptol=1024*eps;
    ptol=1.e-3;
    jtol=1.e-1;
    d=max(max(xo)-min(xo),max(yo)-min(yo))*ptol;
        seg=reshape(conn(:,[1,2,1,3,2,3])',2,3*size(conn,1))';
        seg=unique(sort(seg,2),'rows');
    d=0.1*min(abs(diff(xo(seg)+1i*yo(seg),1,2)));
%    d=.1;
    %[~,ind,jnd]=unique(round([xo,yo]/d)*d,'rows','first','legacy');
    [~,ind,jnd]=unique(round([xo,yo]/d)*d,'rows','stable');
    %[ind,ss]=sort(ind);
    try
        %jnd=ss(jnd);
        conn=reshape(jnd(conn),size(conn));
        xo=xo(ind);
        yo=yo(ind);
    catch
        return
    end
    
    if size(conn,2)==3
        detJ=((xo(conn(:,2))-xo(conn(:,1))).*(yo(conn(:,3))-yo(conn(:,1)))-(yo(conn(:,2))-yo(conn(:,1))).*(xo(conn(:,3))-xo(conn(:,1))))/2;
        flip=detJ<0;
        conn(flip,[1,2])=conn(flip,[2,1]);
        remove=abs(detJ)<100*0.5*0.5*d*d*jtol;
        conn(remove,:)=[];
    end
    [pind,~,j1]=unique(conn);
    conn=reshape(j1,size(conn));
    xo=xo(pind);
    yo=yo(pind);
else
    d=0.1;
    [~,ind,jnd]=unique(round([xo,yo,zo]/d)*d,'rows','stable');
    %[ind,ss]=sort(ind);
    try
        %jnd=ss(jnd);
        conn=reshape(jnd(conn),size(conn));
        xo=xo(ind);
        yo=yo(ind);
        zo=zo(ind);
    catch
        return
    end
    
    [pind,~,j1]=unique(conn);
    conn=reshape(j1,size(conn));
    xo=xo(pind);
    yo=yo(pind);
    zo=zo(pind);
end
end