function nod2=AddOneNodeLayer(conn,nod1)
selected=zeros(max(conn(:)),1);
for in=1:length(nod1)
    [ie,ine]=find(conn==nod1(in));
    neighboors=conn(ie,:);
    keep=find(neighboors(:)>0);
    selected(neighboors(keep))=1;
end
nod2=find(selected);

% outside=zeros(Nnodes);
% outside(nod1)=1;
% outside_nodes=find(outside(:));
% [i1s,i2s]=ind2sub(Nnodes,outside_nodes);
% i1s=min(Nnodes(1),max(1,[i1s-1,i1s,i1s+1,i1s-1,i1s,i1s+1,i1s-1,i1s,i1s+1]));
% i2s=min(Nnodes(2),max(1,[i2s-1,i2s-1,i2s-1,i2s,i2s,i2s,i2s+1,i2s+1,i2s+1]));
% outside_nodes=sub2ind(Nnodes,i1s,i2s);
% outside(outside_nodes)=1;
% nod2=find(outside(:));
end
