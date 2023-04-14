function VI = interpMesh(mesh, V, coords,keepnan)
if nargin<4, keepnan=0;end
    %
    % VI = interpMesh(mesh, V, coords)
    %
    % mesh.xo
    %     .yo
    %     .zo
    %     .conn
    %     .elt
    %
    % V
    % 
    % coords.xi
    %       .yi
    %       .zi
    zflag=isfield(mesh, 'zo');
if zflag
    mesh.elt(mesh.elt==4)=5;
end
VI = NaN(length(coords.xi), size(V,2));
% tic
% for in=1:length(coords.xi)
%     
%     dist=(mesh.xo-coords.xi(in)).^2+(mesh.yo-coords.yi(in)).^2;
%     if zflag
%         dist=dist+(mesh.zo-coords.zi(in)).^2;
%     end
%     dist=sqrt(dist);
%     [dmin,id]=min(dist);
%     
%     
% end
% for i1=1:size(mesh.conn,1)
% 
%     inods=mesh.conn(i1,1:(mesh.elt(i1)-(mesh.elt(i1)==5)));
%     xn=mesh.xo(inods);
%     yn=mesh.yo(inods);
%     if zflag
%         zn=mesh.zo(inods);
%     end
% 
%     Vn=V(inods,:);
% 
%     cxi=find(min(xn) <= coords.xi & coords.xi <= max(xn));
%     cyi=find(min(yn) <= coords.yi & coords.yi <= max(yn));
%     if isfield(mesh, 'zo')
%         czi=find(min(zn) <= coords.zi & coords.zi <= max(zn));
%     end
%     
%     if isfield(mesh, 'zo')
%         ci = intersect(intersect(cxi,cyi),czi);
%     else
%         ci = intersect(cxi,cyi);
%     end
%     
%     if ~isempty(ci)
%         switch mesh.elt(i1)
%             case 8
%             [xg,yg,zg]=GetGaussPointsHexaedron(0,[1,1,1],xn,yn,zn,coords.xi(ci),coords.yi(ci),coords.zi(ci));
%             N=[0.125*(1-xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1+yg).*(1-zg),0.125*(1-xg).*(1+yg).*(1-zg),...
%                 0.125*(1-xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1+yg).*(1+zg),0.125*(1-xg).*(1+yg).*(1+zg)];
%             case 5
%             [xg,yg,zg]=GetGaussPointsTetrahedron(0,1,xn,yn,zn,coords.xi(ci),coords.yi(ci),coords.zi(ci));
%             N=[1-xg-yg-zg,xg,yg,zg];
%             case 4
%             [xg,yg]=GetGaussPointsQuadrangle(0,1,xn,yn,coords.xi(ci),coords.yi(ci));
%             N=[0.25*(1-xg).*(1-yg),0.25*(1+xg).*(1-yg),0.25*(1+xg).*(1+yg),0.25*(1-xg).*(1+yg)];
%             case 3
%             [xg,yg]=GetGaussPointsTriangle(0,1,xn,yn,coords.xi(ci),coords.yi(ci));
%             N=[1-xg-yg,xg,yg];
%         end
%         for di=1:size(V,2)
%             Vin=N*Vn(:,di);
%             found=~isnan(Vin);
%             VI(ci(found),di)=Vin(found);
%         end
%     end
% end
% 
% toc

%tic


for i1=1:size(mesh.conn,1)

    inods=mesh.conn(i1,1:(mesh.elt(i1)-(mesh.elt(i1)==5)));
    xn=mesh.xo(inods);
    yn=mesh.yo(inods);
    if zflag
        zn=mesh.zo(inods);
    end

    Vn=V(inods,:);
if 1
    cxi=find(min(xn) <= coords.xi & coords.xi <= max(xn));
    cyi=find(min(yn) <= coords.yi & coords.yi <= max(yn));
    if zflag
        czi=find(min(zn) <= coords.zi & coords.zi <= max(zn));
    end
    
    if zflag
        ci = intersect(intersect(cxi,cyi),czi);
    else
        ci = intersect(cxi,cyi);
    end
else
    
%    ci=(min(xn) <= coords.xi & coords.xi <= max(xn))&(min(yn) <= coords.yi & coords.yi <= max(yn));
%    if zflag
%        ci=ci&(min(zn) <= coords.zi & coords.zi <= max(zn));
%    end
 ci=inpolygon(coords.xi,coords.yi,[xn;xn(1)],[yn;yn(1)]);   
    ci=find(ci);
end
    
    
    
    if ~isempty(ci)
%    if any(ci)
        switch mesh.elt(i1)
            case 8
            [xg,yg,zg,wg]=GetGaussPointsVoxels(8,xn,yn,zn,coords.xi(ci),coords.yi(ci),coords.zi(ci));
            N=[0.125*(1-xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1-yg).*(1-zg),0.125*(1+xg).*(1+yg).*(1-zg),0.125*(1-xg).*(1+yg).*(1-zg),...
                0.125*(1-xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1-yg).*(1+zg),0.125*(1+xg).*(1+yg).*(1+zg),0.125*(1-xg).*(1+yg).*(1+zg)];
            case 5
            [xg,yg,zg,wg]=GetGaussPointsVoxels(4,xn,yn,zn,coords.xi(ci),coords.yi(ci),coords.zi(ci));
            N=[1-xg-yg-zg,xg,yg,zg];
            case 4
            [xg,yg,wg]=GetGaussPointsPixels(4,xn,yn,coords.xi(ci),coords.yi(ci));
            N=[0.25*(1-xg).*(1-yg),0.25*(1+xg).*(1-yg),0.25*(1+xg).*(1+yg),0.25*(1-xg).*(1+yg)];
            case 3
            [xg,yg,wg]=GetGaussPointsPixels(3,xn,yn,coords.xi(ci),coords.yi(ci));
            N=[1-xg-yg,xg,yg];
        end
        for di=1:size(V,2)
            Vin=N*Vn(:,di);
            found=wg>0;
            VI(ci(found),di)=Vin(found);
        end
    end
end
if any(isnan(VI(:,1)))&&~keepnan
    found=find(isnan(VI(:,1)));
    for in=1:length(found)
        inod=found(in);
        xn=coords.xi(inod);
        yn=coords.yi(inod);
        if zflag
            zn=coords.zi(inod);
            [dmin,~]=min((mesh.xo-xn).^2+(mesh.yo-yn).^2+(mesh.zo-zn).^2);
            idmin=find(((mesh.xo-xn).^2+(mesh.yo-yn).^2+(mesh.zo-zn).^2)<dmin+2);
        else
            [dmin,~]=min(((mesh.xo-xn).^2+(mesh.yo-yn).^2));
            idmin=find(((mesh.xo-xn).^2+(mesh.yo-yn).^2)<dmin+2);
            
        end
        VI(inod,:)=mean(V(idmin,:),1);
        
    end
end
end