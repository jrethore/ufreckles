function [elt,conn,xo,yo,Nnodes,Nelems,nselected,zo,Uf]=BuiltRefinedConnectivity(elt1,conn1,xo,yo,selected,fac,zo,Ug)
if nargin<6, fac=1;end
if nargin<7, zo=0*xo;end
Uf=[];
elt=elt1;conn=conn1;
if nargin==8
    Ux=Ug((1:length(xo)),:,:);
    Uy=Ug(length(xo)+(1:length(xo)),:,:);
    dflag=size(Ug,1)>2*length(xo);
    if dflag
        Uz=Ug(2*length(xo)+(1:length(xo)),:,:);
    else
        Uz=0*Ux;
    end
end
for ir=1:fac
    foundq4=find(elt==4);
    foundt3=find(elt==3);
    %    assert(length(foundt3)*length(foundq4)==0)
    eltn=zeros(4*length(foundt3)+4*length(foundq4),1);
    connn=repmat(0,length(eltn),4);
    xa=[];ya=[];za=[];
    Uxa=[];Uya=[];Uza=[];
    nelt=0;
    for ie=1:length(elt)
        inods=conn(ie,1:elt(ie));
        xn=xo(inods);
        yn=yo(inods);
        zn=zo(inods);
        
        if elt(ie)==4
            nf=zeros(4,1);
            
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
            for i1=1:elt(ie)
                iface=conn(ie,mod(i1+(0:1)-1,elt(ie))+1);
                xface=mean(xo(iface));
                yface=mean(yo(iface));
                zface=mean(zo(iface));
                found=find(xa==xface&ya==yface&za==zface);
                if isempty(found)||(~selected(iface(1))&&~selected(iface(2)))
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
            end
            
            for i1=1:elt(ie)
                %                 ie
                %                 i1
                %                 ng
                %                 nf'
                %                 iface=conn(ie,mod(i1+(0:1)-1,elt(ie))+1)
                %                 inods=[ng,nf(i1),conn(ie,mod(i1-1+1,4)+1),nf(mod(i1+1-1,elt(ie))+1)]
                %                 xt=[xo;xa];
                %                 yt=[yo;ya];
                %                 xn'
                %                 yn'
                %                 xt(inods)'
                %                 yt(inods)'
                %                 pause
                connn(nelt+1,1:elt(ie))=[ng,nf(i1),conn(ie,mod(i1-1+1,4)+1),nf(mod(i1+1-1,elt(ie))+1)];
                eltn(nelt+1)=elt(ie);
                nelt=nelt+1;
            end
        elseif elt(ie)==3
            nf=zeros(3,1);
            for i1=1:elt(ie)
                iface=conn(ie,mod(i1+(0:1)-1,elt(ie))+1);
                xface=mean(xo(iface));
                yface=mean(yo(iface));
                zface=mean(zo(iface));
                found=find(xa==xface&ya==yface&za==zface);
                if isempty(found)||(~selected(iface(1))&&~selected(iface(2)))
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
            end
            % figure
            % triplot(conn(:,1:3),xo,yo)
            % hold on
            % plot(xa,ya,'r+')
            % pause
            for i1=1:elt(ie)
                %                                 ie
                %                                 i1
                %                                 nf'
                %                                 iface=conn(ie,mod(i1+(0:1)-1,elt(ie))+1)
                %                                 inods=[nf(i1),conn(ie,mod(i1-1+1,3)+1),nf(mod(i1+1-1,elt(ie))+1)]
                %
                %
                %                                 xt=[xo;xa];
                %                                 yt=[yo;ya];
                %                                 xn'
                %                                 yn'
                %                                 xt(inods)'
                %                                 yt(inods)'
                %                                  triplot([1:3],xt(inods),yt(inods),'r')
                %                                pause
                
                connn(nelt+1,1:elt(ie))=[nf(i1),conn(ie,mod(i1-1+1,3)+1),nf(mod(i1+1-1,elt(ie))+1)];
                eltn(nelt+1)=elt(ie);
                nelt=nelt+1;
            end
            
            connn(nelt+1,1:elt(ie))=nf';
            eltn(nelt+1)=elt(ie);
            nelt=nelt+1;
            %         inods=nf;
            %                                   triplot([1:3],xt(inods),yt(inods),'k')
            %                                pause
            
            
            %             for i1=1:elt(ie)
            %                 connn(nelt+1,1:elt(ie))=[ng,conn(ie,i1),conn(ie,mod(i1,elt(ie))+1)];
            %                 eltn(nelt+1)=elt(ie);
            %                 nelt=nelt+1;
            %             end
            
            
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