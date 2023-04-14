function EE=MedianFilterCell(conn,E)
persistent cond
EE=zeros(size(E));
% if size(conn,2)>3
%     assert(~any(conn(:,4)>0));
% end
if size(conn,2)>3
    if any(conn(:,4)>0)
        
        elt=4*ones(size(conn,1),1);
    else
        conn=conn(:,1:3);
        elt=3*ones(size(conn,1),1);
    end
else
    conn=conn(:,1:3);
    elt=3*ones(size(conn,1),1);
end
ne=numel(elt);
nn=max(conn(:));
if isempty(cond)
    ie=repmat((1:ne)',1,size(conn,2));
    cond=sparse(conn(:),ie(:),1,nn,ne);
    cond=cond'*cond;
end
comp=size(E,1)/numel(elt);
for ie=1:size(conn,1)
    nelt=find(cond(ie,:)>0);
    for icomp=0:(comp-1)
        for in=1:size(E,2)
            EE(ie+icomp*numel(elt),in)=median(E(nelt+icomp*numel(elt),in));
        end
    end
    
end



end