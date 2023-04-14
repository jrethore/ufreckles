function [elts,S]=GetEltsFromNodes(conn,elt,nodes,crit)
if any(nodes==0),nodes=find(nodes);end
if nargin<4, crit=0;end
elts=[];
S=cell(4,1);
if crit
    for i1=1:numel(elt)
        go=1;
        for in=1:elt(i1)
            if ~any(nodes==conn(i1,in))
                go=0;
            end
        end
        if go
            elts=[elts,i1];
        end
    end
    
else
    for i1=1:numel(elt)
        go=0;
        for in=1:elt(i1)
            if any(nodes==conn(i1,in))
                go=1;
            end
        end
        if go
            elts=[elts,i1];
        end
    end
end
if nargout>1
    for i1=1:length(elts)
        ie=elts(i1);
        inods=conn(ie,1:elt(ie));
        is=zeros(elt(ie),1);
        for in=1:elt(ie)
            is(in)=any(nodes==inods(in));
        end
        if sum(is)>1
        if diff(find(is))==1
            id=min(find(is));
        else
            id=max(find(is));
        end
        S{id}=[S{id},ie];
        end
    end
    
end


end

