function UpdateFEMVICMeshes(dU,nmods,iter)
iscale=1;
nel=0;
for im=1:length(nmods)
    nmod=nmods(im);
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
    lc=param.transition_length;
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'n','xo','yo','zo','Nnodes','Nelems','conn','elt');
    Ui=dU(nel+(1:prod(Nnodes)));
%     figure
%     scatter3(xo,yo,zo,10+0*xo,Ui)
%    hold on 
%     quiver3(xo,yo,zo,n(:,1).*Ui,n(:,2).*Ui,n(:,3).*Ui)
%     figure
%     plot3(xo,yo,zo,'bx')
%     hold on
    xo=xo+n(:,1).*Ui;
    yo=yo+n(:,2).*Ui;
    zo=zo+n(:,3).*Ui;
%     plot3(xo,yo,zo,'r+')
%    return
    Xo=xo;Yo=yo;Zo=zo;
    save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Xo','Yo','Zo','xo','yo','zo','-append');
    nel=nel+prod(Nnodes);    
    n=zeros(prod(Nnodes),3);
    for in=1:prod(Nnodes)
        nconn=find(conn==in);
        nn=0;
        for ic=1:length(nconn)
            [ie,ien]=ind2sub(size(conn),nconn(ic));
            inods=conn(ie,1:elt(ie));
            xn=xo(inods(1+mod(ien+[-1,0,1]-1,elt(ie))));
            yn=yo(inods(1+mod(ien+[-1,0,1]-1,elt(ie))));
            zn=zo(inods(1+mod(ien+[-1,0,1]-1,elt(ie))));
            nx=(yn(1)-yn(2))*(zn(3)-zn(2))-(zn(1)-zn(2))*(yn(3)-yn(2));
            ny=(zn(1)-zn(2))*(xn(3)-xn(2))-(xn(1)-xn(2))*(zn(3)-zn(2));
            nz=(xn(1)-xn(2))*(yn(3)-yn(2))-(yn(1)-yn(2))*(xn(3)-xn(2));
            nn=nn+[nx,ny,nz];
        end
        nn=nn/length(nconn);
        nn=nn/norm(nn);
        n(in,:)=-nn;
    end
    xo=[xo;xo+lc*n(:,1)];
    yo=[yo;yo+lc*n(:,2)];
    zo=[zo;zo+lc*n(:,3)];
    Nnodes=[length(xo),1,1];
    
    
    save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'n','-append');
    save(fullfile('TMP',sprintf('%d_3d_vicmesh_%d',nmod,iscale-1)),'xo','yo','zo','-append');

    unix(sprintf( 'cp %s.mat %s.mat',fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),fullfile('TMP',sprintf('%d_mesh_%d_%03d',nmod,iscale-1,iter))));
        writeVTKmesh25D(fullfile('TMP',sprintf('%d_mesh_%d_%03d',nmod,iscale-1,iter)));

end

end