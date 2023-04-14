function UpdateVirtualImage(U,iscale,nmod,refine)
if nargin<4,refine=0;end
if iscale>0
    load(fullfile('TMP','params'),'param');
    tau=param.transition_length;
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
thickness=0;
if isfield(param,'line_thickness')
    thickness=round(0.5*(param.line_thickness));
end

    load(fullfile('TMP',sprintf('sample0_%d',1-1)),'im0');


    load(fullfile('TMP',sprintf('sample0_%d',1-1)),'lso','ls1','sizeim','nband','on','nx','ny');
    load(fullfile('TMP',sprintf('%d_phi_%d',nmod,(iscale-1))),'phii','dphii');
    [xon,yon]=ind2sub(sizeim,nband(on));
    nxon=-nx(nband(on));
    nyon=-ny(nband(on));
    unx=nxon.*(phii*U);
    uny=nyon.*(phii*U);
    
        figure
        plot(xon,yon,'bx')
        hold on
        plot(xon+unx,yon+uny,'ro')
    
    xon=xon+unx;
    yon=yon+uny;

    M=phii'*phii;
    L=phii'*xon;
    xs=M\L;
    L=phii'*yon;
    ys=M\L;
    xon=phii*xs;
    yon=phii*ys;
    
    
    plot(xon,yon,'kx')
    
    dxon=dphii*xs;
    dyon=dphii*ys;
    nxon=dyon;
    nyon=-dxon;
    nnorm=abs(nxon+1i*nyon);
    nxon=nxon./nnorm;
    nyon=nyon./nnorm;
    
    %     txon=dxon;
    %     tyon=dyon;
    %     txon=txon./nnorm;
    %     tyon=tyon./nnorm;
         quiver(xon,yon,nxon,nyon,'k')
    %     quiver(xon,yon,txon,tyon,'g')
         axis equal
[lso,ls1,nmesh]=ComputeLevelSetFromPoints(nmod,iscale,xon,yon,nxon,nyon,refine);
    nx=FDgradient(lso,1);
    ny=FDgradient(lso,2);
    switch param.contour_type
        case 'edge'
            im0=255*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
            nband=find((lso(:)>=0)&(lso(:)<=tau)&nmesh);
        case 'line'
            im0=255*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
            nband=find((abs(lso(:))<=(tau+thickness))&nmesh);
    end

    on=find((lso(nband)<1)&(lso(nband)>=0));
    % figure
    % imagesc(im0)
%     [xon,yon]=ind2sub(sizeim,nband);
%     figure
%     plot(xon,yon,'x')
%     hold on
%     quiver(xon,yon,nx(nband),ny(nband),'r')

    save(fullfile('TMP','sample0_0'),'im0','nx','ny','lso','ls1','nband','on','-append');
    %      if iscale>2
    %         tau=2^(iscale-2)*tau;
    for iscale=2:param.nscale
        tau=2*tau;

        switch param.contour_type
            case 'edge'
                im0=255*0.5*(1-cos(pi*min(lso,tau)/tau)).*double(lso>=0);
                nband=find((lso(:)>=0)&(lso(:)<=tau)&nmesh);
            case 'line'
            im0=255*(0.5*(1+cos(pi*max(0,min(abs(lso)-thickness,tau))/tau)));
            nband=find((abs(lso(:))<=(tau+thickness))&nmesh);
        end

        on=find((lso(nband)<1)&(lso(nband)>=0));
%     [xon,yon]=ind2sub(sizeim,nband);
%     figure
%     plot(xon,yon,'x')
%     hold on
%     quiver(xon,yon,nx(nband),ny(nband),'r')

        save(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim','nband','on','-append');
    end
end
end