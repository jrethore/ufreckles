function model=Remeshing(param,model,morf)
if nargin<3,morf=0;end
if abs(morf)>0
    check=1;
    clear functions
    nmod=0;
    nuo=model.material_parameters.nu;
    model.material_parameters.nu=0.;
    LoadParameters(param);
    LoadParameters(model,nmod);
%    ReferenceImage(nmod);
    LoadMask(nmod);
    LoadMeshes(nmod);
    LoadMat(nmod);
    
    CreateGradBasisFunction(1,nmod);
    AssembleMechanicalOperator(1,nmod);
    roi=param.roi;
    load(fullfile('TMP','0_mesh_0.mat'))
    xo=xo+roi(1)-1;
    yo=yo+roi(3)-1;
    zo=xo+1i*yo;
    
    conn=conn(:,1:3);
    TR=triangulation(conn,xo,yo);
    load(fullfile('TMP','0_k_operator_0.mat'),'K')
    
    be=freeBoundary(TR);
    
    bnodes=zeros(prod(Nnodes),1);
    nx=zeros(prod(Nnodes),1);
    ny=zeros(prod(Nnodes),1);
    
    for ib=1:size(be,1)
        inods=be(ib,:);
        bnodes(inods)=1;
        ny(inods)=ny(inods)-diff(xo(inods));
        nx(inods)=nx(inods)+diff(yo(inods));
    end
    nn=abs(nx+1i*ny);
    
    
    tipnodes=zeros(prod(Nnodes),1);
    facenodes=zeros(prod(Nnodes),1);
    tx=zeros(prod(Nnodes),1);
    ty=zeros(prod(Nnodes),1);
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        if zone{4}==5
            if (zone{8}>0)
                indc=zone{10};
                cnodes=indc{1};
                ztip=zo(cnodes(1,1));
                tt=-diff(zo(cnodes(1:2,1)));
                tt=tt/abs(tt);
                rtip=mean(abs(zo(indc{2})-ztip));
                tipnodes(abs(zo-ztip)<=rtip)=1;
                tipnodes(indc{2})=1;
                bnodes(abs(zo-ztip)<=rtip)=0;
                bnodes(indc{2})=0;
                bnodes(cnodes(:,2))=0;
                facenodes(cnodes(:,1))=cnodes(:,2);
                facenodes(abs(zo-ztip)<=rtip)=0;
                facenodes(indc{2})=0;
                
                tx(abs(zo-ztip)<=rtip)=real(tt);
                ty(abs(zo-ztip)<=rtip)=imag(tt);
                tx(indc{2})=real(tt);
                ty(indc{2})=imag(tt);
                if  (zone{9}>0)
                    ztip=zo(cnodes(end,1));
                    tt=diff(zo(cnodes(end+(-1:0),1)));
                    tt=tt/abs(tt);
                    rtip=mean(abs(zo(indc{3})-ztip));
                    tipnodes(abs(zo-ztip)<=rtip)=1;
                    tipnodes(indc{3})=1;
                    bnodes(abs(zo-ztip)<=rtip)=0;
                    bnodes(indc{3})=0;
                    facenodes(abs(zo-ztip)<=rtip)=0;
                    facenodes(indc{3})=0;
                    tx(abs(zo-ztip)<=rtip)=real(tt);
                    ty(abs(zo-ztip)<=rtip)=imag(tt);
                    tx(indc{3})=real(tt);
                    ty(indc{3})=imag(tt);
                end
                
            end
        end
    end
    bnodes(~selected)=0;
    indi=find(bnodes>0);
    nx(indi)=nx(indi)./nn(indi);
    ny(indi)=ny(indi)./nn(indi);
    Cpx=sparse(indi,1:length(indi),nx(indi),2*prod(Nnodes),length(indi));
    Cpy=sparse(indi+prod(Nnodes),1:length(indi),ny(indi),2*prod(Nnodes),length(indi));
    Cp=Cpx+Cpy;
    Un=zeros(length(indi),1);
    
    facenodes(~selected)=0;
    indi=find(facenodes>0);
    Cfx=sparse(indi,1:length(indi),1,2*prod(Nnodes),length(indi))-sparse(facenodes(indi),1:length(indi),1,2*prod(Nnodes),length(indi));
    Cfy=sparse(indi+prod(Nnodes),1:length(indi),1,2*prod(Nnodes),length(indi))-sparse(facenodes(indi)+prod(Nnodes),1:length(indi),1,2*prod(Nnodes),length(indi));
    Cf=[Cfx,Cfy];
    Uf=zeros(2*length(indi),1);
    
    
    indi=find(tipnodes>0);
    Ctx=sparse(indi,1:length(indi),tx(indi),2*prod(Nnodes),length(indi));
    Cty=sparse(indi+prod(Nnodes),1:length(indi),ty(indi),2*prod(Nnodes),length(indi));
    Cnx=sparse(indi,1:length(indi),-ty(indi),2*prod(Nnodes),length(indi));
    Cny=sparse(indi+prod(Nnodes),1:length(indi),tx(indi),2*prod(Nnodes),length(indi));
    Ct=[Ctx+Cty,Cnx+Cny];
    
    Up=[ones(length(indi),1);zeros(length(indi),1)];
    
    indi=find(~selected);
    indi=[indi;indi+prod(Nnodes)];
    Cb=sparse(indi,1:length(indi),1,2*prod(Nnodes),length(indi));
    Ub=zeros(length(indi),1);
    
    a=1e8*max(diag(K));
    
    C=[Ct,Cb,Cp,Cf];
    Kc=[K,C;C',zeros(size(C,2),size(C,2))];
    Fc=[zeros(2*prod(Nnodes),1);Up;Ub;Un;Uf];
    
    X=Kc\Fc;
    load(fullfile('TMP','0_mesh_0.mat'))
    if check
        U=morf*X(1:2*prod(Nnodes),:);
        save('morphing.mat','U','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
        postproVTK(['morphing'],0,0);
    end
    xo=xo+roi(1)-1+morf*X(1:prod(Nnodes));
    yo=yo+roi(3)-1+morf*X(prod(Nnodes)+(1:prod(Nnodes)));
    model.material_parameters.nu=nuo;
else
    keyboard
    load(fullfile('TMP','sample0'),'sizeim')
    display('Remeshing.....')
    check=1;
    ho=model.mesh_size;
    h=round(ho/2);
    model.mesh_size=h;
    roio=param.roi;
    roi=[1,sizeim(1),1,sizeim(2)];
    roi=roi+0.5*[h(1),-h(1),h(2),-h(2)];
    [ye,xe]=meshgrid(roi(3):h(2):roi(4),roi(1):h(1):roi(2));
    [d,xyfix,indc]=GetMeshDensity(model,xe(:),ye(:));
    hmin=min(d);
    ld=@(xy) GetSignedDistanceToZone(model,roi,xy(:,1),xy(:,2));
    lh=@(xy) (GetMeshDensity(model,xy(:,1),xy(:,2)));
    [xo,yo,conn,xyfix]=GenMeshFromLVL7(roi,2*mean(h)*hmin,ld,lh,xyfix,1*check);
    
    try
        if check
            close all
            % figure
            % scatter(xe(:),ye(:),d*0+10,d)
            % axis equal
            figure
            triplot(conn,xo,yo)
            axis equal
            hold on
        end
    catch
        keyboard
    end
    out=feval(ld,xyfix)>(0.001*2*mean(h)*hmin);
    new_id=0*out;
    new_id(~out)=1:sum(~out);
    for iz=1:size(model.zone,2)
        if model.zone{4,iz}==5
            if check
                plot(model.zone{2,iz}*[1;1i],'r-x')
                pause(0.1)
            end
            for id=1:size(indc,1)
                nodes=indc{id,iz};
                if ~isempty(nodes)
                    nodes=new_id(nodes);
                    if id==1,tipin=[nodes(1),nodes(end)]>0;end
                    nodes(nodes==0)=[];
                    
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
                                    if detJ<0;
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
                        toadd=(1+tipin(1)):(size(nodes,1)-tipin(2)*model.zone{9,iz});
                        nodes(toadd,2)=length(xo)+(1:length(toadd));
                        xo=[xo;xo(nodes(toadd,1))];
                        yo=[yo;yo(nodes(toadd,1))];
                        
                        for ip=1:length(toadd)
                            %                         inods=nodes(toadd(ip),1);
                            %                         seg=-t(toadd(ip)-1);
                            %                         mid=(xo(inods)+1i*yo(inods));
                            %                         elts=sum(conn==inods,2)>0;
                            %                         xc=sum(xo(conn(elts,:)),2)/3;
                            %                         yc=sum(yo(conn(elts,:)),2)/3;
                            %                         side=real(((xc+1i*yc)-mid)*(seg*exp(1i*pi/2))');
                            %                         elts=find(elts);
                            %
                            %
                            %                         for ii=1:length(side)
                            %                             if side(ii)<0
                            %                                 nelt=conn(elts(ii),:);
                            %                                 nelt(nelt==inods)=nodes(toadd(ip),2);
                            %                                 conn(elts(ii),:)=nelt;
                            %                             end
                            %
                            %                         end
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
    Nnodes=size(xo);
    Nelems=[size(conn,1),1,1];
    elt=3*ones(size(conn,1),1);
    selected=ones(Nnodes);
    for iz=1:size(model.zone,2)
        
        switch  model.zone{4,iz}
            case 5
                model.zone{10,iz}=indc(:,iz);
            case 6
                zone=model.zone{2,iz};
                in=inpolygon(xo,yo,zone(:,1),zone(:,2));
                in=find(in);
                xon=xo(in);
                yon=yo(in);
                model.zone{5,iz}=[xon,yon];
                model.zone{6,iz}=in;
                if any(sum(isnan(model.zone{7,iz})))
                    selected(in)=0;
                end
        end
    end
    model.mesh_size=ho;
    display('Remeshing.....done')
end
save(strrep(param.result_file,'.res','.dat'),'param','model','xo','yo','Nnodes','Nelems','conn','elt','selected')
writeVTKmesh(strrep(param.result_file,'.res','.dat'));

end