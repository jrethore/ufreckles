function [nodes,conn]=GetNodesFromElts(conn,elt,sel)
if nargin<3,sel=1:numel(elt);end
selected=zeros(max(conn(:)),1);
for i1=1:numel(sel)
    selected(conn(sel(i1),1:elt(i1)))=1;
end
nodes=find(selected);
conn(conn==0)=length(selected)+1;
selected=[selected;0];
selected(selected>0)=1:length(nodes);
conn=selected(conn);
end
