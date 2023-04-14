function GenerateBezierTriangles(nmod)
 load(fullfile('TMP',sprintf('%d_params',nmod)),'param');

 iscale=1;
 load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'elt');
 assert(~any(elt>3),'The initial mesh must contains trianlge only !!');
 load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'conn','xo','yo','selected','ns');
 
 zo=0*xo;
 selected=selected(:);
 
    eltn=zeros(length(elt),1)+10;
    connn=repmat(0,length(eltn),10);
    xa=[];ya=[];za=[];

     for ie=1:length(elt)
        inods=conn(ie,1:elt(ie));
        xn=xo(inods);
        yn=yo(inods);
        zn=zo(inods);

             nf=zeros(1,3*2+1);
% ie
% inods
% xn
% yn

            for i1=1:elt(ie)
                iface=conn(ie,mod(i1+(0:1)-1,elt(ie))+1);
                xface=xo(iface(1))+diff(xo(iface))*[1/3;2/3];
                yface=yo(iface(1))+diff(yo(iface))*[1/3;2/3];
                zface=zo(iface(1))+diff(zo(iface))*[1/3;2/3];
%                 i1
%                 iface
%                 xface
%                 yface
%                 figure
%                 plot(xn,yn,'x')
%                 hold on
%                 plot(xface,yface,'rx')
                for i2=1:2
                found=find(abs(xa-xface(i2))<1&abs(ya-yface(i2))<1&abs(za-zface(i2))<1);
                if isempty(found)
                    xa=[xa;xface(i2)];
                    ya=[ya;yface(i2)];
                    za=[za;zface(i2)];
                    found=length(xa);
                 end
                nf(i2+(i1-1)*2)=found;
                 if ~any(selected(iface))
                    selected(found+length(xo))=0;
                else
                    selected(found+length(xo))=1;
                 end
               end
            end
                    xa=[xa;mean(xn)];
                    ya=[ya;mean(yn)];
                    za=[za;mean(zn)];
                    found=length(xa);
            nf(7)=found;
            connn(ie,:)=[conn(ie,1:elt(ie)),nf+length(xo)];
     end
     
         elt=eltn;
    conbt=connn;
%     figure
%     plot(xo,yo,'bo')
%     hold on
%     plot(xa,ya,'rx')
    xo=[xo;xa];
    yo=[yo;ya];
    zo=[zo;za];
Nbsnodes=[length(xo),1,1];
Nbselems=[length(elt),1,1];
Px=xo;
Py=yo;

elt=3*ones(length(eltn)*9,1);
conni=conn;
conn=zeros(length(eltn)*9,3);

    iconn=[1,4,9;4,10,9;4,5,10;5,6,10;5,2,6;9,10,8;10,7,8;10,6,7;8,7,3];
for ie=1:size(conbt,1)
    inods=conbt(ie,:);
    ies=(1:9)+(ie-1)*9;
    conn(ies,:)=inods(iconn);

    
    
end
cons=conn;
%     figure
%     hold on
%     axis equal
% for ie=1:size(conn,1)
%  plot(xo(conn(ie,[1:3,1])),yo(conn(ie,[1:3,1])),'b')
% end


[elt,conn,xo,yo,Nnodes,Nelems,nselected,zo]=BuiltRefinedConnectivity(elt,conn,xo,yo,selected,2,zo);

%     figure
%     hold on
%     axis equal
% for ie=1:size(conn,1)
%  plot(xo(conn(ie,[1:3,1])),yo(conn(ie,[1:3,1])),'b')
% end
ng=1;
ns=round(0.5*param.mesh_size);
mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
save(mesh_file,'ns','ng','Px','Py','conbt','cons','Nbsnodes','Nbselems','Nnodes','Nelems','conn','elt','-append');
[phio,xo,yo]=CreateBezierTriangleBasis(mesh_file,'nodes');
save(mesh_file,'xo','yo','zo','-append');
save(fullfile('TMP',sprintf('%d_phio_%d',nmod,iscale-1)),'phio');
 end