function run_fdfea_job(param,model)
%persistent Ko
Ko=[];
warning('on')
checkc=0;
nmod=0;
iscale=1;
param.onflight=1;
sizeim=param.roi([2,4])-param.roi([1,3])+1;
model.nscale=1;
LoadParameters(param);
LoadParameters(model,nmod);
save(fullfile('TMP','sample0_0.mat'),'sizeim')
save(fullfile('TMP','sample0.mat'),'sizeim')
LoadMask(nmod);
LoadMeshes(nmod);

filres=param.result_file;
%%
load(fullfile('TMP','0_mesh_0.mat'),'xo','yo','Nnodes','conn','elt')
conno=conn;
elto=elt;
conn=conn(:,1:3);
extract=0;
ii=[1,2,3];
jj=circshift(ii,[0,-1]);
kk=circshift(ii,[0,-2]);
tips=zeros(size(model.zone,2),2);
cracks=[];
ntip=0;

for iz=1:size(model.zone,2)
    zone=model.zone(:,iz);
    switch zone{4}
        case 5 %CRACK
            xo=xo+param.roi(1)-1;
            yo=yo+param.roi(3)-1;
            if (zone{8}>0)
                extract=1;
                tips(iz,1)=ntip+1;
                ntip=ntip+1;
                if  (zone{9}>0)
                    tips(iz,2)=ntip+1;
                    ntip=ntip+1;
                end
            end
            xyc=zone{2};
            hmax=2*mean(model.mesh_size);
            box=[floor(min(xyc(:,1))-hmax),ceil(max(xyc(:,1))+hmax),floor(min(xyc(:,2))-hmax),ceil(max(xyc(:,2))+hmax)];
            %            htip=real(zone{7});
            ftip=0.5;
            xg=mean(xo(conn),2);
            yg=mean(yo(conn),2);
            inbox=find(inpolygon(xg,yg,box([1,2,2,1,1,]),box([3,3,4,4,3])));
            seg=reshape(conn(inbox,[1,2,1,3,2,3])',2,3*size(inbox,1))';
            seg=unique(sort(seg,2),'rows');
            htip=min(abs(diff(xo(seg),[],2)+1i*diff(yo(seg),[],2)));
            [crack,front]=GetSignedDistanceToCrack(xyc,xg(inbox)+1i*yg(inbox));
            inbox=inbox(((abs(crack)<hmax)&(front<hmax)));
            inods=conn(inbox,:);
            [inods,~,ic]=unique(inods);
            [cn,fn,nx,ny]=GetSignedDistanceToCrack(xyc,xo(inods)+1i*yo(inods));
            tooclose=find((abs(cn)<ftip*htip)&(fn<0));
            nmove=numel(tooclose);
            
            if nmove>0
                intext=GetSignedDistanceToZone(model,param.roi,xo(inods(tooclose)),yo(inods(tooclose)));
                %             xo(inods(tooclose))=xo(inods(tooclose))+ftip*htip*(2*rand(nmove,1)-1);
                %             yo(inods(tooclose))=yo(inods(tooclose))+ftip*htip*(2*rand(nmove,1)-1);
                xo(inods(tooclose))=xo(inods(tooclose))+(intext<-ftip*htip).*abs(abs(cn(tooclose))-ftip*htip).*nx(tooclose).*sign(cn(tooclose));
                yo(inods(tooclose))=yo(inods(tooclose))+(intext<-ftip*htip).*abs(abs(cn(tooclose))-ftip*htip).*ny(tooclose).*sign(cn(tooclose));
            end
            [cn,fn]=GetSignedDistanceToCrack(xyc,xo(inods)+1i*yo(inods));
            cn=reshape(cn(ic),numel(inbox),3);
            fn=reshape(fn(ic),numel(inbox),3);
            edges=[1,2;2,3;3,1];
            face_nodes=[];
            face_elts=[];
            for ij=1:length(inbox)
                ie=inbox(ij);
                inods=conn(ie,:);
                cns=cn(ij,:);
                if any(prod(sign(cns(edges)),2)<0)&&(any(fn(ij,:)<0))
                    face_elts(end+1)=ie;
                    if any(fn(ij,:)>0)
                        xn=xo(inods);
                        yn=yo(inods);
                        A=0.5*(xn(2)*yn(3)-xn(3)*yn(2)+xn(3)*yn(1)-xn(1)*yn(3)+xn(1)*yn(2)-xn(2)*yn(1));
                        a=xn(jj).*yn(kk)-xn(kk).*yn(jj);
                        b=yn(jj)-yn(kk);
                        c=xn(kk)-xn(jj);
                        xg1=xyc(1,1);yg1=xyc(1,2);
                        L1=0.5*((1+0*xg1)*a'+xg1*b'+yg1*c')/A;
                        xg2=xyc(end,1);yg2=xyc(end,2);
                        L2=0.5*((1+0*xg2)*a'+xg2*b'+yg2*c')/A;
                        %                         figure
                        %                         plot(xyc(:,1),xyc(:,2),'k')
                        %                         hold on
                        %                         plot(xn,yn,'x')
                        %                         keyboard
                        if ~(~any(L1<0)||~any(L2<0))
                            if sum(fn(ij,:)>0)>1
                                skip=1;
                                for ed=1:elt(ie)
                                    if (prod(sign(cns(edges(ed,:))))<0)&&(~any(fn(ij,edges(ed,:))>0))
                                        face_nodes=unique([face_nodes,inods(edges(ed,:))]);
                                        skip=0;
                                    end
                                end
                                if skip, face_elts(end)=[];end
                            else
                                face_nodes=unique([face_nodes,inods]);
                            end
                        else
                            if any(L2<0)
                                id=1;
                            end
                            if any(L1<0)
                                id=size(xyc,1);
                            end
                            xfront=xyc(id,:);
                            for ed=1:elt(ie)
                                if (prod(sign(cns(edges(ed,:))))<0)&&(any(fn(ij,edges(ed,:))>0))
                                    ce=cns(edges(ed,:));
                                    xe=xn(edges(ed,:));
                                    ye=yn(edges(ed,:));
                                    a=ce(end)/diff(ce);
                                    xfront=xfront+0.99*([a*xe(1)+(1-a)*xe(2),a*ye(1)+(1-a)*ye(2)]-xfront);
                                    
                                end
                            end
                            xyc(id,:)=xfront;
                            
                        end
                        
                    else
                        face_nodes=unique([face_nodes,inods]);
                    end
                    
                end
                
            end
            if checkc
                figure
                triplot(conn(inbox,:),xo,yo)
                axis equal
                hold on
                plot(xyc(:,1),xyc(:,2),'k')
                plot(xo(face_nodes),yo(face_nodes),'om')
                plot(xg(face_elts),yg(face_elts),'xr')
            end
            %             figure
            %             trimesh(conn(:,1:3),xo,yo,[],crack)
            %             hold on
            %             plot(xyc(:,1),xyc(:,2),'k')
            %             axis equal
            %             figure
            %             trimesh(conn(:,1:3),xo,yo,[],front)
            %             hold on
            %             plot(xyc(:,1),xyc(:,2),'k')
            %             axis equal
            
            new_nodes=prod(Nnodes)+(1:length(face_nodes))';
            new_ids=(1:prod(Nnodes))';
            new_ids(face_nodes)=new_nodes;
            inods=conn(face_elts,:);
            nconn=new_ids(inods);
            [inods,~,ic]=unique(inods);
            [cn,~]=GetSignedDistanceToCrack(xyc,xo(inods)+1i*yo(inods));
            cn=reshape(cn(ic),numel(face_elts),3);
            nconn1=conn(face_elts,:).*(cn>0)+nconn.*(cn<0);
            nconn2=conn(face_elts,:).*(cn<0)+nconn.*(cn>0);
            xo=[xo;xo(face_nodes)];
            yo=[yo;yo(face_nodes)];
            Nnodes=[numel(xo),1,1];
            model.zone{10,iz}=face_elts';
            model.zone{11,iz}=[face_nodes',new_nodes];
            model.zone{12,iz}={nconn1,nconn2};
                        xo=xo-param.roi(1)+1;
            yo=yo-param.roi(3)+1;

    end
end
LoadParameters(model,nmod);
save(fullfile('TMP','0_mesh_0.mat'),'xo','yo','Nnodes','-append')
selected=[];
load(strrep(param.result_file,'.res','.dat'),'-mat','selected')
if isempty(selected)
    fromdic=0;
    
    indi=[];Up=[];
    Fo=zeros(2*prod(Nnodes),1);
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        switch zone{4}
            case 6 % BCS
                loads=zone{7};
                nodes=zone{6};
                for ii=1:2
                    Fo(nodes+(ii-1)*prod(Nnodes))=loads(ii,2)/length(nodes);
                    up=loads(ii,1);
                    if ~isnan(up)
                        indi=[indi;nodes+(ii-1)*prod(Nnodes)];
                        Up(end+(1:length(nodes)))=up;
                    end
                end
        end
    end
else
    fromdic=1;
    load(strrep(param.result_file,'.res','.dat'),'-mat','U')
    indi=find(~selected);
    indio=[indi(:);indi(:)+size(U,1)/2];
    indi=[indi(:);indi(:)+prod(Nnodes)];
    Up=U(indio,1);
    
    Fo=zeros(2*prod(Nnodes),1);
    clear U
    
end

C=sparse(indi,1:length(indi),1,2*prod(Nnodes),length(indi));
Up=Up(:);
Nddl=2*prod(Nnodes);
areg=0;
if isfield(model.material_parameters,'lc')
    if abs(model.material_parameters.lc)>0
        areg=model.material_parameters.lc;
    end
end
%%
LoadMat(nmod);
if isempty(Ko)
    CreateGradBasisFunction(1,nmod);
    clear AssembleStiffnessMatrix
    Ko=AssembleStiffnessMatrix(1,nmod);
    K=Ko;
    dn=1;
else
    indi=(1:size(Ko,1)/2);
    indi=[indi,indi+prod(Nnodes)];
    indj=1:size(Ko,1);
    dn=sparse(indi,indj,1,Nddl,size(Ko,1));
    K=(dn*Ko)*dn';
end

mesh_file=fullfile('TMP','0_mesh_0');
if areg>0
    [phi]=CreateFiniteElementBasis(mesh_file,sizeim,1,[],'Gauss_points');
    load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'Uxx','wdetJ');
    load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'Uyy');
    load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uyx','Uxy');
    phio=sparse(size(phi,1),size(phi,2));
    ux=[phi,phio];
    uy=[phio,phi];
    Mgo=(ux'*(wdetJ*ux)+uy'*(wdetJ*uy));
    Kgo=Uxx'*(wdetJ*Uxx)+Uxy'*(wdetJ*Uxy)+Uyx'*(wdetJ*Uyx)+Uyy'*(wdetJ*Uyy);
end

for iz=1:size(model.zone,2)
    zone=model.zone(:,iz);
    switch zone{4}
        case 5 %CRACK
            xyc=zone{2};
            face_elts=model.zone{10,iz};
            nodes=model.zone{11,iz};
            nconn=model.zone{12,iz};
            conn=conno(face_elts,:);
            elt=elto(face_elts);
            Nelems=[numel(elt),1,1];
            save(fullfile('TMP','0_mesh_0.mat'),'conn','elt','Nelems','-append')
            
            maskg=1;
            save(fullfile('TMP','0_mask_0'),'maskg','-append')
            CreateGradBasisFunction(1,nmod);
            clear AssembleStiffnessMatrix
            Ke=AssembleStiffnessMatrix(1,nmod);
            
            K=K-Ke;
            
            
            if areg>0
                [phi]=CreateFiniteElementBasis(mesh_file,sizeim,1,[],'Gauss_points');
                load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'Uxx','wdetJ');
                load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'Uyy');
                load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uyx','Uxy');
                phio=sparse(size(phi,1),size(phi,2));
                ux=[phi,phio];
                uy=[phio,phi];
                Mge=(ux'*(wdetJ*ux)+uy'*(wdetJ*uy));
                Kge=Uxx'*(wdetJ*Uxx)+Uxy'*(wdetJ*Uxy)+Uyx'*(wdetJ*Uyx)+Uyy'*(wdetJ*Uyy);
                Mg=Mgo-Mge;
                Kg=Kgo-Kge;
            end
            
            
            
            
            ns=10*[1,1];
            for ii=1:2
                conn=nconn{ii};
                save(fullfile('TMP','0_mesh_0.mat'),'conn','ns','-append')
                
                phi=CreateFiniteElementBasis(mesh_file,sizeim,1,[],'sub_cells');
                xg=phi*xo+param.roi(1)-1;
                yg=phi*yo+param.roi(3)-1;
                [maskg,~]=GetSignedDistanceToCrack(xyc,xg+1i*yg);
                switch ii
                    case 1
                        maskg=maskg>0;
                    case 2
                        maskg=maskg<0;
                end
                %                 figure
                %                 triplot(conn,xo,yo)
                %                 axis equal
                %                 hold on
                %                 %             plot(xg,yg,'xr')
                %                 plot(xg(maskg),yg(maskg),'xr')
                
                maskg=diag(sparse(maskg));
                save(fullfile('TMP','0_mask_0'),'maskg','-append')
                [dphidx,dphidy]=CreateGradFiniteElementBasis(mesh_file,sizeim,1,[],'sub_cells');
                [wdetJ,inde]=GetWeigthDetJ(mesh_file,sizeim,1,'sub_cells');
                phi0=sparse(size(dphidx,1),size(dphidx,2));
                epsxx=[dphidx,phi0];
                Uxx=[dphidx,phi0];
                Uxy=[dphidy,phi0];
                Uyx=[phi0,dphidx];
                Uyy=[phi0,dphidy];
                save(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'Uxx','epsxx','wdetJ','sizeim');
                epsyy=[phi0,dphidy];
                save(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'Uyy','epsyy','sizeim');
                epsxy=[dphidy,dphidx];
                save(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'Uyx','Uxy','epsxy','sizeim');
                
                clear AssembleStiffnessMatrix
                Ke=AssembleStiffnessMatrix(1,nmod);
                K=K+Ke;
                if areg>0
                    wdetJ=maskg*wdetJ;
                    ux=[phi,phi0];
                    uy=[phi0,phi];
                    Mge=(ux'*(wdetJ*ux)+uy'*(wdetJ*uy));
                    Kge=Uxx'*(wdetJ*Uxx)+Uxy'*(wdetJ*Uxy)+Uyx'*(wdetJ*Uyx)+Uyy'*(wdetJ*Uyy);
                    Mg=Mg+Mge;
                    Kg=Kg+Kge;
                end
                
            end
            
    end
end


Kc=[K,C;C',sparse(size(C,2),size(C,2))];
Fc=[Fo;Up];
X=Kc\Fc;
U=X(1:Nddl,:);
Fl=K*U;

if areg>0
    Kg=Kg*model.lc+Mg;
    
    Kc=[Kg,C;C',sparse(size(C,2),size(C,2))];
    Fc=[Mg*U;Up];
    X=Kc\Fc;
    U=X(1:Nddl,:);
    Fl=Kg*U;
    
    
end

if areg>0
    Ks=[];
else
    [Ks]=SIF(nmod,model,U,xo+param.roi(1)-1,yo+param.roi(3)-1,conno(:,1:3));
end
Ks=Ks;
%%
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
model=param;
%%
for iz=1:size(model.zone,2)
    uu=zeros(2,size(U,2));
    ff=zeros(2,size(U,2));
    zone=model.zone(:,iz);
    switch zone{4}
        case 6
            nodes=zone{6};
            for ii=1:2
                indo=nodes+(ii-1)*prod(Nnodes);
                ff(ii,:)=sum(Fl(indo,:),1);
                uu(ii,:)=mean(U(indo,:),1);
            end
            model.zone{7,iz}=[uu;ff];
    end
end
load(fullfile('TMP','params'),'param');
conn=conno;
elt=elto;
for iz=1:size(model.zone,2)
    zone=model.zone(:,iz);
    switch zone{4}
        case 5 %CRACK
            face_elts=model.zone{10,iz};
            nconn=model.zone{12,iz};
            conn(face_elts,:)=[nconn{1},zeros(numel(face_elts),1)];
            conn=[conn;[nconn{2},zeros(numel(face_elts),1)]];
            elt=[elt;3*ones(numel(face_elts),1)];
    end
end

Nelems=[numel(elt),1,1];
load(fullfile('TMP','0_mesh_0'),'Nnodes','xo','yo','ng','rflag','rint');
try
    tips(~(cell2mat(model.zone(4,:))==5),:)=[];
    cracks=zeros(size(model.zone,2),1);
    iscrack=find(cell2mat(model.zone(4,:))==5);
    cracks(iscrack)=iscrack;
    cracks(~((cell2mat(model.zone(4,:))==5)|(cell2mat(model.zone(4,:))==6)))=[];
    model.zone(:,~((cell2mat(model.zone(4,:))==5)|(cell2mat(model.zone(4,:))==6)))=[];
catch
end
save(filres,'U','Fl','Ks','tips','cracks','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
postproVTK([filres,''],0,0);

end