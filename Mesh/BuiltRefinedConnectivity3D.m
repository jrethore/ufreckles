function [elt,conn,xo,yo,zo,Nnodes,Nelems,nselected,Uf]=BuiltRefinedConnectivity3D(elt1,conn1,xo,yo,zo,selected,fac,Ug)

error('not finished')
if nargin<7, fac=1;end
Uf=[];
elt=elt1;conn=conn1;
if nargin==8
    Ux=Ug((1:length(xo)),:,:);
    Uy=Ug(length(xo)+(1:length(xo)),:,:);
    Uz=Ug(2*length(xo)+(1:length(xo)),:,:);
end
for ir=1:fac
    foundq4=find(elt==8);
    foundt3=find(elt==6);
    if ~isempty(foundt3)
        error('not coded yet')
    end
    eltn=zeros(8*length(foundt3)+8*length(foundq4),1);
    connn=repmat(0,length(eltn),8);
    xa=[];ya=[];za=[];
    Uxa=[];Uya=[];Uza=[];
    nelt=0;
    for ie=1:length(elt)
        inods=conn(ie,1:elt(ie));
        xn=xo(inods);
        yn=yo(inods);
        zn=zo(inods);
        
        if elt(ie)==8
            nf=zeros(12,1);
            
            xa=[xa;mean(xn)];
            ya=[ya;mean(yn)];
            za=[za;mean(zn)];
            if nargin==8
                Uxa=[Uxa;mean(Ux(inods,:,:),1)];
                Uya=[Uya;mean(Uy(inods,:,:),1)];
                Uza=[Uza;mean(Uz(inods,:,:),1)];
            end
            ng=length(xa)+length(xo);
            selected(ng)=1;
            for i1=1:4
                iface=conn(ie,mod(i1+(0:1)-1,4)+1);
                xface=mean(xo(iface));
                yface=mean(yo(iface));
                zface=mean(zo(iface));
                found=find(xa==xface&ya==yface&za==zface);
                if isempty(found)
                    xa=[xa;xface];
                    ya=[ya;yface];
                    za=[za;zface];
                    if nargin==8
                        Uxa=[Uxa;mean(Ux(iface,:,:),1)];
                        Uya=[Uya;mean(Uy(iface,:,:),1)];
                        Uza=[Uza;mean(Uz(iface,:,:),1)];
                    end
                    found=length(xa);
                end
                nf(i1)=found+length(xo);
                if ~any(selected(iface))
                    selected(found+length(xo))=0;
                else
                    selected(found+length(xo))=1;
                end
                iface=conn(ie,4+mod(i1+(0:1)-1,4)+1);
                xface=mean(xo(iface));
                yface=mean(yo(iface));
                zface=mean(zo(iface));
                found=find(xa==xface&ya==yface&za==zface);
                if isempty(found)
                    xa=[xa;xface];
                    ya=[ya;yface];
                    za=[za;zface];
                    if nargin==8
                        Uxa=[Uxa;mean(Ux(iface,:,:),1)];
                        Uya=[Uya;mean(Uy(iface,:,:),1)];
                        Uza=[Uza;mean(Uz(iface,:,:),1)];
                    end
                    found=length(xa);
                end
                nf(i1+4)=found+length(xo);
                if ~any(selected(iface))
                    selected(found+length(xo))=0;
                else
                    selected(found+length(xo))=1;
                end
                
                iface=conn(ie,i1+[0,4]);
                xface=mean(xo(iface));
                yface=mean(yo(iface));
                zface=mean(zo(iface));
                found=find(xa==xface&ya==yface&za==zface);
                if isempty(found)
                    xa=[xa;xface];
                    ya=[ya;yface];
                    za=[za;zface];
                    if nargin==8
                        Uxa=[Uxa;mean(Ux(iface,:,:),1)];
                        Uya=[Uya;mean(Uy(iface,:,:),1)];
                        Uza=[Uza;mean(Uz(iface,:,:),1)];
                    end
                    found=length(xa);
                end
                nf(i1+8)=found+length(xo);
                if ~any(selected(iface))
                    selected(found+length(xo))=0;
                else
                    selected(found+length(xo))=1;
                end
            end
            
                connn(nelt+1,1:elt(ie))=[ng,nf(i1+8),conn(ie,mod(i1-1+1,8)+1),nf(mod(i1+1-1,elt(ie))+1)];
                connn(nelt+1,1:elt(ie))=[ng,nf(i1+8),conn(ie,mod(i1-1+1,8)+1),nf(mod(i1+1-1,elt(ie))+1)];
                eltn(nelt+[1:8])=elt(ie);
                nelt=nelt+8;
           
        elseif elt(ie)==6
            
        end
        
    end
    
    elt=eltn;
    conn=connn;
    xo=[xo;xa];
    yo=[yo;ya];
    zo=[zo;za];
    if nargin==8
        Ux=[Ux;Uxa];
        Uy=[Uy;Uya];
        Uz=[Uz;Uza];
    end
end
Nnodes=[length(xo),1,1];
Nelems=[length(elt),1,1];
nselected=selected;
if nargin==8
    if dflag
        Uf=[Ux;Uy;Uz];
    else
        Uf=[Ux;Uy];
    end
end
end