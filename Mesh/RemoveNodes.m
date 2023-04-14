function U2=RemoveNodes(nmod,toremove,U1)
iscale=1;
check=1;
load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','elt','conn','Nnodes','Nelems','selected');
if check
    figure
    plot(xo,yo,'k+')
    hold on
     axis equal
   plot(xo(~selected),yo(~selected),'ro')
end
found=find(~toremove);
    nconn=zeros(size(conn));
    nelt=zeros(size(elt));
    nselected=ones(length(found),1);
    ielt=0;
    for ie=1:length(elt)
        inods=conn(ie,1:elt(ie));
        xn=xo(inods);
        yn=yo(inods);
        founde=find(~toremove(inods));
        if length(founde)==elt(ie)
            for in=1:elt(ie)
                nconn(ielt+1,in)=find(found==inods(in));
            end
            nelt(ielt+1)=elt(ie);
            ielt=ielt+1;
            if any(~selected(inods))
                founds=find(~selected(inods));
                for in=1:length(founds)
                    nselected(find(found==inods(founds(in))))=0;
                end
            end
       elseif any(~selected(inods))
            for in=1:elt(ie)
                if any(founde==in)
                    nselected(find(found==inods(in)))=0;
                end
            end

        end
    end
    if nargin<2
        U2=0;
    else
Ux=U1(found,:);
Uy=U1(found+prod(Nnodes),:);
U2=[Ux;Uy];
    end
    xo=xo(found);
    yo=yo(found);
    conn=nconn(1:ielt,:);
    elt=nelt(1:ielt);
    selected=nselected;
    Smesh=floor([max(xo)-min(xo),max(yo)-min(yo)]);
Vmesh=floor([max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)]);
Nnodes=[length(xo),1,1];
Nelems=[size(conn,1),1,1];

if check
    plot(xo,yo,'bx')
    plot(xo(~selected),yo(~selected),'gs')
end
    
    
save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','-append');

end