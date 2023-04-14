function [xo,yo,conn,indc]=GenTriMesh(fem_model,roi)

h=fem_model.mesh_size;
roi=roi+0.5*[h(1),-h(1),h(2),-h(2)];

xe=roi(1):h:roi(2);
ye=roi(3):h:roi(4);
[ye,xe]=meshgrid(ye,xe);


[d,xyfixo,indc]=GetMeshDensity(fem_model,xe(:),ye(:));

    ld=@(xy) GetSignedDistanceToZone(fem_model,roi,xy(:,1),xy(:,2));

lh=@(xy) (GetMeshDensity(fem_model,xy(:,1),xy(:,2)));
[xo,yo,conn,xyfix]=GenMeshFromLVL7(roi,mean(h)*min(d),ld,lh,xyfixo);
if ~isempty(indc)
    out=feval(ld,xyfix)>=(0.001*mean(h)*min(d));
    new_id=0*out;
    new_id(~out)=1:sum(~out);
    
    
    duplicated=[];
    for iz=1:size(fem_model.zone,2)
        if fem_model.zone{4,iz}==5
            for id=1:size(indc,1)
                nodes=indc{id,iz}';
                if ~isempty(nodes)
                    for in=1:numel(nodes)
                        dist=abs((xyfix(nodes(in),1)-xo)+1i*(xyfix(nodes(in),2)-yo));
                        [dmin,idmin]=min(dist);
                        if dmin>(0.001*mean(h)*min(d))
                            nodes(in)=0;
                        else
                            nodes(in)=idmin;
                        end
                    end
                    nodes(nodes==0)=[];
                    if id==1,tipin=[nodes(1),nodes(end)]>0;end
                    segc=[nodes(1:end-1),nodes(2:end)];
                    eneighboor=GetEltsFromNodes(conn,3*ones(size(conn,1),1),nodes);
                    neighboor=conn(eneighboor,:);
                    segs=reshape(neighboor(:,[1,2,1,3,2,3])',2,3*size(neighboor,1))';
                    segs=unique(sort(segs,2),'rows');
                    
                    [~,isc,~]=intersect(segc,segs,'rows');
                    segc(isc,:)=[];
                    toinsert=zeros(size(segc,1),1);
                    for ic=1:size(segc,1)
                        has1=neighboor((sum(neighboor==segc(ic,1),2)>0),:);
                        has2=neighboor((sum(neighboor==segc(ic,2),2)>0),:);
                        nodes12=intersect(has1(:),has2(:))';
                        switch numel(nodes12)
                            case 2
                                for ii=1:2
                                    elti=sort([segc(ic,:),nodes12(ii)],2);
                                    try
                                        [~,ie,~]=intersect(sort(neighboor,2),sort([nodes12,segc(ic,ii)],2),'rows');
                                    catch
                                        keyboard
                                    end
                                    detJ=((xo(elti(2))-xo(elti(1)))*(yo(elti(3))-yo(elti(1)))-(yo(elti(2))-yo(elti(1)))*(xo(elti(3))-xo(elti(1))))/2;
                                    if detJ<0
                                        elti([1,2])=elti([2,1]);
                                    end
                                    conn(eneighboor(ie),:)=elti;
                                end
                            case 1
                                toinsert(ic)=nodes12;
                        end
                    end
                    if any(toinsert)
                        nnodes=zeros(length(nodes)+sum(toinsert>0),1);
                        for ic=1:size(segc,1)
                            if toinsert(ic)>0
                                in=find(nodes==segc(ic,1));
                                nnodes(in+1)=toinsert(ic);
                            end
                        end
                        in=1;
                        for ii=1:length(nnodes)
                            if nnodes(ii)==0
                                nnodes(ii)=nodes(in);
                                in=in+1;
                            end
                        end
                        nodes=nnodes;
                    end
                    
                    
                    if id==1
                        nodes=repmat(nodes,1,2);
                        toadd=(1+tipin(1)):(size(nodes,1)-tipin(2)*fem_model.zone{9,iz});
                        if ~isempty(duplicated)
                            ldup=[];
                            inn=numel(toadd)-1;
                            inod=nodes(toadd(inn),1);
                            elts=sum(conn==inod,2)>0;inods=conn(elts,:);
                            
                            for ii=1:numel(inods)
                                if any(duplicated(:,2)==inods(ii))
                                    try
                                        nodes(toadd(end),:)=duplicated(duplicated(:,1)==nodes(toadd(end),1),2);
                                    catch
                                    end
                                    break
                                end
                            end
                        end
                        nodes(toadd,2)=length(xo)+(1:length(toadd));
                        xo=[xo;xo(nodes(toadd,1))];
                        yo=[yo;yo(nodes(toadd,1))];
                        duplicated=[duplicated;nodes(toadd,:)];
                        for ip=1:length(toadd)
                            
                            
                            inods=nodes(toadd(ip),1);
                            zn=(xo(inods)+1i*yo(inods));
                            elts=sum(conn==inods,2)>0;
                            elts=find(elts);
                            for ii=1:length(elts)
                                xc=sum(xo(conn(elts(ii),:)))/3;
                                yc=sum(yo(conn(elts(ii),:)))/3;
                                
                                if toadd(ip)==1
                                    behind=0;
                                    infront=1;
                                elseif toadd(ip)==size(nodes,1)
                                    behind=1;
                                    infront=0;
                                else
                                    t=diff(xo(nodes(toadd(ip)+(0:1),1))+1i*yo(nodes(toadd(ip)+(0:1),1)));
                                    infront=(real(((xc+1i*yc)-zn)*(t)'))>0;
                                    t=diff(xo(nodes(toadd(ip)+(-1:0),1))+1i*yo(nodes(toadd(ip)+(-1:0),1)));
                                    behind=(real(((xc+1i*yc)-zn)*(t)'))<0;
                                end
                                if infront
                                    t=diff(xo(nodes(toadd(ip)+(0:1),1))+1i*yo(nodes(toadd(ip)+(0:1),1)));
                                elseif behind
                                    t=diff(xo(nodes(toadd(ip)+(-1:0),1))+1i*yo(nodes(toadd(ip)+(-1:0),1)));
                                else
                                    t=diff(xo(nodes(toadd(ip)+(0:1),1))+1i*yo(nodes(toadd(ip)+(0:1),1)));
                                    t=t+diff(xo(nodes(toadd(ip)+(-1:0),1))+1i*yo(nodes(toadd(ip)+(-1:0),1)));
                                end
                                n=-t*exp(1i*pi/2);
                                side=real(((xc+1i*yc)-zn)*(n)');
                                if side<0
                                    nelt=conn(elts(ii),:);
                                    nelt(nelt==inods)=nodes(toadd(ip),2);
                                    conn(elts(ii),:)=nelt;
                                end
                                
                            end
                            
                        end
                    else
                        if ~isempty(nodes)
                            if tipin(mod(id,2)+1)
                                face_nodes=indc{1,iz};
                                tipid=(mod(id,2)==0)+size(face_nodes,1)*(mod(id,2)==1);
                                ztip=xo(face_nodes(tipid,1))+1i*yo(face_nodes(tipid,1));
                                zzone=xo(nodes)+1i*yo(nodes);
                                zface=xo(face_nodes(:,1))+1i*yo(face_nodes(:,1));
                                rtip=mean(abs(zzone-ztip));
                                [dmin,idmin]=min(abs(abs(zface-ztip)-rtip));
                                nodes(end+1)=face_nodes(idmin,1);
                                nodes(end+1)=face_nodes(idmin,2);
                            end
                        end
                    end
                    indc{id,iz}=nodes;
                end
            end
        end
    end
    
end

