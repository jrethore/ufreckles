function [Uini]=InitializeSolution25D(Uo,iscale,nmod)

load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nnodes');
Nddl=prod(Nnodes);
if size(Uo,1)==3*prod(Nnodes)
    Uini=Uo;
else

    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale)),'unmasked_nodes');
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale)),'xo','yo','Nnodes','conn','elt');
    xg=2*(xo-1)+1;
    yg=2*(yo-1)+1;
    meshg.conn=conn;
    meshg.elt=elt;
    meshg.xo=xg;
    meshg.yo=yg;

    xg=reshape(xg,Nnodes);
    yg=reshape(yg,Nnodes);
    if isempty(unmasked_nodes)
        Ug=reshape(Uo((1:Nddl)),Nnodes(1:2));
        Vg=reshape(Uo(Nddl+(1:Nddl)),Nnodes(1:2));
        Wg=reshape(Uo(2*Nddl+(1:Nddl)),Nnodes(1:2));
    else
        Nddl=length(unmasked_nodes);
        Ug=repmat(0,Nnodes(1:2));
        Vg=repmat(0,Nnodes(1:2));
        Wg=repmat(0,Nnodes(1:2));
        Ug(unmasked_nodes)= Uo((1:Nddl));
        Vg(unmasked_nodes)= Uo(Nddl+(1:Nddl));
        Wg(unmasked_nodes)= Uo(2*Nddl+(1:Nddl));
    end

    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','Nnodes');
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');




    xo=min(xo,max(xg(:)));xo=max(xo,min(xg(:)));
    yo=min(yo,max(yg(:)));yo=max(yo,min(yg(:)));
 
%     Uf=interp2(yg,xg,Ug,yo(:),xo(:),'linear');
%     Vf=interp2(yg,xg,Vg,yo(:),xo(:),'linear');
%     Wf=interp2(yg,xg,Wg,yo(:),xo(:),'linear');

    coords.xi=xo;
    coords.yi=yo;
    [UVWf]=interpMesh(meshg,[Ug,Vg,Wg],coords);
    Uf=UVWf(:,1);
    Vf=UVWf(:,2);
    Wf=UVWf(:,3);

    if isempty(unmasked_nodes)
        Uini=[Uf(:);Vf(:);Wf(:)];
    else
        Uini=[Uf(unmasked_nodes);Vf(unmasked_nodes);Wf(unmasked_nodes)];
    end

end
end