function [phi,Xi,Yi,wdetJ]=CreateBezierTriangleBasis(mesh_file,GaussPts)

load(mesh_file,'Px','Py','Nbsnodes','Nbselems','conbt','cons','ng','ns');
if nargin < 2 , GaussPts='pixels';end

assert(ng>0);
Po=[[0.;1;0],[0;0;1]];
Po=[3*Po;2*Po(1,:)+Po(2,:);Po(1,:)+2*Po(2,:);2*Po(2,:)+Po(3,:);Po(2,:)+2*Po(3,:);2*Po(3,:)+Po(1,:);Po(3,:)+2*Po(1,:);Po(1,:)+Po(2,:)+Po(3,:)]/3;
xo=Po(:,1);
yo=Po(:,2);
wo=1-xo-yo;
iconn=[1,4,9;4,10,9;4,5,10;5,6,10;5,2,6;9,10,8;10,7,8;10,6,7;8,7,3];
switch GaussPts
    case {'pixels','Gauss_points'}
        if strcmp(GaussPts,'pixels')
            ngt=ng*prod(ns)*2;
            [xgt,ygt,wgt]=GetGaussPointsTriangle(max(1,ngt/(2*prod(ns))),ns);
            Nt_r=[-1+0*xgt,1+0*xgt,0*ygt];
            Nt_s=[-1+0*ygt,0*xgt,1+0*ygt];
            Selt=length(xgt);
            Nt=[1-xgt-ygt,xgt,ygt];
            xg=zeros(length(xgt)*size(iconn,1),1);
            yg=zeros(length(xgt)*size(iconn,1),1);
            wgs=zeros(length(xgt)*size(iconn,1),1);
            for ie=1:size(iconn,1)
                inods=iconn(ie,:);
                xn=xo(inods);
                yn=yo(inods);
                dxdr=Nt_r*xo(inods);
                dydr=Nt_r*yo(inods);
                dxds=Nt_s*xo(inods);
                dyds=Nt_s*yo(inods);
                wge=wgt.*(dxdr.*dyds-dydr.*dxds);
                ig=(1:Selt)+(ie-1)*Selt;
                xg(ig)=Nt*xn;
                yg(ig)=Nt*yn;
                wgs(ig)=wge;
            end
            
        elseif strcmp(GaussPts,'Gauss_points')
            [xg,yg,wgs]=GetGaussPointsTriangle(25);
        end
        
        wg=1-xg-yg;
        N=[wg.^3,xg.^3,yg.^3,...
            3*xg.*wg.^2,3*wg.*xg.^2,...
            3*yg.*xg.^2,3*xg.*yg.^2,...
            3*wg.*yg.^2,3*yg.*wg.^2,...
            6*wg.*xg.*yg];
        N_r=[-3*wg.^2,3*xg.^2,0*yg,...
            3*wg.^2-6*xg.*wg,-3*xg.^2+6*wg.*xg,...
            6*yg.*xg,3*yg.^2,...
            -3*yg.^2,-6*yg.*wg,...
            6*wg.*yg-6*xg.*yg];
        N_s=[-3*wg.^2,0*xg,3*yg.^2,...
            -6*xg.*wg,-3*xg.^2,...
            3*xg.^2,6*xg.*yg,...
            -3*yg.^2+6*wg.*yg,3*wg.^2-6*yg.*wg,...
            6*wg.*xg-6*xg.*yg];
        
        Selt=length(xg);
        npts=10*length(xg)*prod(Nbselems);
        ny=prod(Nbsnodes);
        nx=Selt*prod(Nbselems);
        indn=zeros(npts,1);
        indp=zeros(npts,1);
        val=zeros(npts,1);
        wdetJ=zeros(nx,1);
        nel=0;
        npix=0;
        for i1=1:prod(Nbselems)
            inods=conbt(i1,1:10);
            dxdr=N_r*Px(inods);
            dydr=N_r*Py(inods);
            dxds=N_s*Px(inods);
            dyds=N_s*Py(inods);
            detJ=(dxdr.*dyds-dydr.*dxds);
            
            for in=1:10
                id=inods(in);
                indn(nel+(1:Selt))=id;
                indp(nel+(1:Selt))=npix+(1:Selt);
                val(nel+(1:Selt))=N(:,in);
                nel=nel+Selt;
            end
            wdetJ(npix+(1:Selt))=wgs.*detJ;
            npix=npix+Selt;
        end
        
        phi=sparse(indp,indn,val,nx,ny);
        wdetJ=diag(sparse(wdetJ));
    case 'nodes'
        wdetJ=1;
        load(mesh_file,'conn');
        [eltr,iconnr,xg,yg,nNnodes,nNelems,nselected]=BuiltRefinedConnectivity(3*ones(size(iconn,1),1),iconn,xo,yo,1+0*xo,2);
        wg=1-xg-yg;
        
        N=[wg.^3,xg.^3,yg.^3,...
            3*xg.*wg.^2,3*wg.*xg.^2,...
            3*yg.*xg.^2,3*xg.*yg.^2,...
            3*wg.*yg.^2,3*yg.*wg.^2,...
            6*wg.*xg.*yg];
        
        npts=10*3*size(iconnr,1)*prod(Nbselems);
        ny=prod(Nbsnodes);
        nx=max(conn(:));
        indn=zeros(npts,1);
        indp=zeros(npts,1);
        val=zeros(npts,1);
        nn=zeros(nx,1);
        nel=0;
        for i1=1:prod(Nbselems)
            inods=conbt(i1,1:10);
            for isub=1:size(iconnr,1)
                insub=iconnr(isub,1:3);
                Nsub=N(insub,:);
                isubg=isub+(i1-1)*size(iconnr,1);
                ipix=conn(isubg,1:3);
                for in=1:10
                    id=inods(in);
                    indn(nel+(1:3))=id;
                    indp(nel+(1:3))=ipix;
                    val(nel+(1:3))=Nsub(:,in);
                    nel=nel+3;
                end
                nn(ipix)=nn(ipix)+(abs(sum(Nsub,2))>0);
            end
        end
        nn=1./max(1,nn);
        phi=sparse(indp,indn,val,nx,ny);
        
        phi=diag(sparse(nn))*phi;
        
end
Xi=phi*Px(:);
Yi=phi*Py(:);
end
