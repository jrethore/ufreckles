function [Uini3D,M,P]=ExtrudeNURBSModel(nmod,Uini,dooperator)
P=0;
M=0;
Uini3D=[];
check=0;
iscale=1;
tic;
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
params=param.extrusion_parameters;

if ~isfield(param0,'calibration_data')
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rint','ns','ng','conn','elt','ui','vi','xo','yo','zo','uo','vo','wo','Px','Py','Nbsnodes','Nnodes','Nbselems','Nelems','Smesh','Vmesh','selected');
    
    wo=0:params.thickness/(params.nlayers-1):(params.thickness);
    Pz=[0,-params.thickness];
    wd=wo([1,1,length(wo),length(wo)]);
    if params.degree>1
        [Pz,wd] = bspdegelev(1,Pz,wd,params.degree-1);
    end
    if length(wo)>2
        [Pz,wnew] = bspkntins(params.degree,Pz,wd,wo(2:length(wo)-1));
        wd=wnew;
    end
    wo=wd;
    pnlayers=length(Pz);
    Px=repmat(Px,[1,1,pnlayers]);
    Py=repmat(Py,[1,1,pnlayers]);
    Pz=repmat(reshape(Pz,[1,1,length(Pz)]),[size(Px,1),size(Px,2),1]);
    Nbsnodes=[size(Px)];
    Nbselems(3)=params.nlayers-1;

    
    
    to=0:params.thickness/((params.nlayers-1)*2^((params.degree)*(params.degree>1))):(params.thickness);
    wi=repmat(to,length(xo),1);
    enlayers=length(to);
    wi=wi(:);
    ui=repmat(ui,enlayers,1);
    vi=repmat(vi,enlayers,1);
    foundt3=find(elt==3);
    foundq4=find(elt==4);
    conno=conn;
    conn=repmat(0,size(conno,1),8);
    conn(foundt3,1:3)=conno(foundt3,1:3)+prod(Nnodes);
    conn(foundt3,4:6)=conno(foundt3,1:3);
    conn(foundq4,1:4)=conno(foundq4,1:4)+prod(Nnodes);
    conn(foundq4,5:8)=conno(foundq4,1:4);
    elt=2*elt;
    if enlayers>2
        connad=repmat(0:(enlayers-2),length(elt),1);
        connad=repmat(connad(:)*prod(Nnodes),1,size(conn,2));
    else
        connad=0;
    end
    elt=repmat(elt,enlayers-1,1);
    conn=repmat(conn,enlayers-1,1);
    conn=conn+connad;
    selected=repmat(selected,enlayers,1);
    
    
    degree=[param.degree,params.degree];
    save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'degree','rint','ng','conn','elt','ui','vi','wi','uo','vo','wo','Px','Py','Pz','Nbsnodes','Nbselems','selected');
    
    [phio,xo,yo,zo]=CreateNURBSBasis3D(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),degree,'nodes');
    Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
    Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
    Nnodes=[length(xo),1,1];
    Nelems=[size(conn,1),1,1];
    ns=ones(1,3);
    save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'xo','yo','zo','ns','Nnodes','Nelems','Smesh','Vmesh','-append');
    save(fullfile('TMP',sprintf('%d_3d_phio_%d',nmod,iscale-1)),'phio');

    if ~isempty(Uini)
        Uxini=repmat(Uini((1:size(Uini,1)/2),:),pnlayers,1);
        Uyini=repmat(Uini(size(Uini,1)/2+(1:size(Uini,1)/2),:),pnlayers,1);
        Uini3D=[Uxini;Uyini;0*Uxini];
    end
    
    
    
    
    
    
    
    
    if dooperator
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'phix');
        load(fullfile('TMP',sprintf('%d_phiy_%d',nmod,10*(iscale-1))),'phiy');
        phix=phix(:,(1:prod(Nbsnodes(1:2))));
        phiy=phiy(:,prod(Nbsnodes(1:2))+(1:prod(Nbsnodes(1:2))));
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'Xi','Yi','wdetJ');
        phi0=sparse(size(phix,1),(pnlayers-1)*prod(Nbsnodes(1:2)));
        phix=[phix,phi0;...
            phi0,phix];
        phiy=[phiy,phi0;...
            phi0,phiy];
        
        phi0=sparse(size(phix,1),(pnlayers)*prod(Nbsnodes(1:2)));
        phix=[phix,phi0,phi0];
        phiy=[phi0,phiy,phi0];
        Xi=repmat(Xi,2,1);
        Yi=repmat(Yi,2,1);
        if numel(wdetJ)>1
            wdetJ=repmat(diag(wdetJ),2,1);
            wdetJ=diag(wdetJ);
        end
        Nddl_tot=size(phix,2);
        save(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'phix','Nddl_tot','Xi','Yi','wdetJ','-append');
        save(fullfile('TMP',sprintf('%d_phiy_%d',nmod,10*(iscale-1))),'phiy','-append');
        load(fullfile('TMP',sprintf('%d_phidf_%d',nmod,iscale-1)),'phidf');
        
        phidfx=phidf(:,(1:prod(Nbsnodes(1:2))));
        phidfy=phidf(:,prod(Nbsnodes(1:2))+(1:prod(Nbsnodes(1:2))));
        phi0=sparse(size(phidfx,1),(pnlayers-1)*prod(Nbsnodes(1:2)));
        phidfx=[phidfx,phi0;...
            phi0,phidfx];
        phidfy=[phidfy,phi0;...
            phi0,phidfy];
        phi0=sparse(size(phix,1),pnlayers*prod(Nbsnodes(1:2)));
        phidf=[phidfx,phidfy,phi0];
        save(fullfile('TMP',sprintf('%d_phidf_%d',nmod,iscale-1)),'phidf');
        AssembleCorrelationOperator(iscale,nmod,1);
        
    end
    if nargout>1
        Mo=diag(sparse(ones(prod(Nbsnodes(1:2)),1)));
        M0=sparse((Nbsnodes(3)-2)*prod(Nbsnodes(1:2)),(Nbsnodes(3)-2)*prod(Nbsnodes(1:2)));
        Mx=blkdiag(Mo,blkdiag(M0,Mo));
        My=blkdiag(Mo,blkdiag(M0,Mo));
        Mz=sparse(size(Mx,1),size(Mx,1));
        M=blkdiag(Mx,blkdiag(My,Mz));
    end
else
    error('NOT YET READY');
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rint','Xo','Yo','Zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
    switch params.type
        case 'shell'
            assert(params.nlayers==2);
            if dooperator
                ncam=size(param0.deformed_image,1);
                for icam=1:ncam
                    load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),10*(iscale-1))),'phix');
                    load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icam-1),10*(iscale-1))),'phiy');
                    phixx=phix(:,(1:prod(Nnodes)));
                    phixy=phix(:,prod(Nnodes)+(1:prod(Nnodes)));
                    phixz=phix(:,2*prod(Nnodes)+(1:prod(Nnodes)));
                    phiyx=phiy(:,(1:prod(Nnodes)));
                    phiyy=phiy(:,prod(Nnodes)+(1:prod(Nnodes)));
                    phiyz=phiy(:,2*prod(Nnodes)+(1:prod(Nnodes)));
                    
                    %                     phi0=sparse(size(phix,1),(params.nlayers-1)*prod(Nnodes));
                    % %                     phixx=0.5*[phixx,phixx];
                    % %                     phixy=0.5*[phixy,phixy];
                    % %                     phixz=0.5*[phixz,-phixz];
                    % %                     phiyx=0.5*[phiyx,phiyx];
                    % %                     phiyy=0.5*[phiyy,phiyy];
                    % %                     phiyz=0.5*[phiyz,-phiyz];
                    %                     phixx=[phixx,0*phixx];
                    %                     phixy=[phixy,0*phixy];
                    %                     phixz=[phixz,-0*phixz];
                    %                     phiyx=[phiyx,0*phiyx];
                    %                     phiyy=[phiyy,0*phiyy];
                    %                     phiyz=[phiyz,-0*phiyz];
                    %
                    %                     phi0=sparse(size(phix,1),(params.nlayers)*prod(Nnodes));
                    %                     phix=[phixx,phixy,phixz];
                    %                     phiy=[phiyx,phiyy,phiyz];
                    %                     Nddl_tot=size(phix,2);
                    %                     Nddl_tot=size(phix,2);
                    %                     save(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),10*(iscale-1))),'phix','Nddl_tot','-append');
                    %                     save(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icam-1),10*(iscale-1))),'phiy','-append');
                    load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),10*(iscale-1))),'Xi','Yi','wdetJ');
                    phi0=sparse(size(phix,1),(params.nlayers-1)*prod(Nnodes));
                    phixx=[phixx,phi0;...
                        phi0,phixx];
                    phixy=[phixy,phi0;...
                        phi0,phixy];
                    phixz=[phixz,phi0;...
                        phi0,phixz];
                    phiyx=[phiyx,phi0;...
                        phi0,phiyx];
                    phiyy=[phiyy,phi0;...
                        phi0,phiyy];
                    phiyz=[phiyz,phi0;...
                        phi0,phiyz];
                    
                    phix=[phixx,phixy,phixz];
                    phiy=[phiyx,phiyy,phiyz];
                    Nddl_tot=size(phix,2);
                    Xi=repmat(Xi,2,1);
                    Yi=repmat(Yi,2,1);
                    if numel(wdetJ)>1
                        wdetJ=repmat(diag(wdetJ),2,1);
                        wdetJ=diag(wdetJ);
                    end
                    save(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),10*(iscale-1))),'phix','Nddl_tot','Xi','Yi','wdetJ','-append');
                    save(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icam-1),10*(iscale-1))),'phiy','-append');
                end
                ComputeGradFPhi25D(iscale,nmod);
                AssembleCorrelationOperator25D(iscale,nmod);
            end
            
            if numel(params.thickness)==length(Xo)
                t=params.thickness;
            else
                t=repmat(params.thickness,length(Xo),1);
            end
            xo=zeros(length(Xo),1);
            yo=zeros(length(Xo),1);
            zo=zeros(length(Xo),1);
            for in=1:prod(Nnodes)
                found=find(sum(conn==in,2));
                n=0;
                %          figure
                %          hold on
                for ie=1:length(found)
                    ielt=found(ie);
                    inods=conn(ielt,1:elt(ielt));
                    Xn=Xo(inods);
                    Yn=Yo(inods);
                    Zn=Zo(inods);
                    M=[Xn,Yn,Zn];
                    F=1+0*Xn;
                    a=M\F;
                    n=n+sign(a(3))*a(1:3)./norm(a(1:3));
                    %             plot3([Xn;Xn(1)],[Yn;Yn(1)],[Zn;Zn(1)])
                end
                n=n/length(found);
                n=n/norm(n);
                %           quiver3(Xo(in),Yo(in),Zo(in),0.01*n(1),0.01*n(2),0.01*n(3))
                %           pause
                xo(in)=Xo(in)+0.5*t(in)*n(1);
                yo(in)=Yo(in)+0.5*t(in)*n(2);
                zo(in)=Zo(in)+0.5*t(in)*n(3);
                xo(in+prod(Nnodes))=Xo(in)-0.5*t(in)*n(1);
                yo(in+prod(Nnodes))=Yo(in)-0.5*t(in)*n(2);
                zo(in+prod(Nnodes))=Zo(in)-0.5*t(in)*n(3);
            end
            
            foundt3=find(elt==3);
            foundq4=find(elt==4);
            conno=conn;
            conn=repmat(0,size(conno,1),8);
            conn(foundt3,1:3)=conno(foundt3,1:3)+prod(Nnodes);
            conn(foundt3,4:6)=conno(foundt3,1:3);
            conn(foundq4,1:4)=conno(foundq4,1:4)+prod(Nnodes);
            conn(foundq4,5:8)=conno(foundq4,1:4);
            elt=2*elt;
            selected=repmat(selected,params.nlayers,1);
            Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
            Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
            Nnodes=[length(xo),1,1];
            Nelems=[size(conn,1),1,1];
            save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
            if nargin>1
                Uxini=repmat(Uini((1:size(Uini,1)/3),:),params.nlayers,1);
                Uyini=repmat(Uini(size(Uini,1)/3+(1:size(Uini,1)/3),:),params.nlayers,1);
                Uzini=repmat(Uini(2*size(Uini,1)/3+(1:size(Uini,1)/3),:),params.nlayers,1);
                Uini3D=[Uxini;Uyini;Uzini];
            end
            if nargout>1
                Mo=diag(sparse(ones(size(Uini,1),1)));
                Nnods=length(Xo);
                inods=(1:length(Xo))';
                indj=1:(3*prod(Nnodes));
                indi=[inods;inods;inods+Nnods;inods+Nnods;inods+2*Nnods;inods+2*Nnods];
                val=0.5;
                P=sparse(indi,indj,val,3*Nnods,3*prod(Nnodes));
                M=P'*Mo*P;
            end
        case {'solid','solid-from-topo'}
            if dooperator
                ncam=size(param0.deformed_image,1);
                for icam=1:ncam
                    load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),10*(iscale-1))),'phix','Xi','Yi','wdetJ');
                    load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icam-1),10*(iscale-1))),'phiy');
                    phixx=phix(:,(1:prod(Nnodes)));
                    phixy=phix(:,prod(Nnodes)+(1:prod(Nnodes)));
                    phixz=phix(:,2*prod(Nnodes)+(1:prod(Nnodes)));
                    phiyx=phiy(:,(1:prod(Nnodes)));
                    phiyy=phiy(:,prod(Nnodes)+(1:prod(Nnodes)));
                    phiyz=phiy(:,2*prod(Nnodes)+(1:prod(Nnodes)));
                    
                    phi0=sparse(size(phix,1),(params.nlayers-1)*prod(Nnodes));
                    phixx=[phixx,phi0;phi0,phixx];
                    phixy=[phixy,phi0;phi0,phixy];
                    phixz=[phixz,phi0;phi0,-phixz];
                    phiyx=[phiyx,phi0;phi0,phiyx];
                    phiyy=[phiyy,phi0;phi0,phiyy];
                    phiyz=[phiyz,phi0;phi0,-phiyz];
                    
                    phix=[phixx,phixy,phixz];
                    phiy=[phiyx,phiyy,phiyz];
                    Nddl_tot=size(phix,2);
                    Xi=repmat(Xi,2,1);
                    Yi=repmat(Yi,2,1);
                    if numel(wdetJ)>1
                        wdetJ=repmat(diag(wdetJ),2,1);
                        wdetJ=diag(wdetJ);
                    end
                    save(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),10*(iscale-1))),'phix','Nddl_tot','Xi','Yi','wdetJ','-append');
                    save(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icam-1),10*(iscale-1))),'phiy','-append');
                end
                ComputeGradFPhi25D(iscale,nmod);
                AssembleCorrelationOperator25D(iscale,nmod);
            end
            if nargout>1
                Mo=diag(sparse(ones(size(Uini,1)/3,1)));
                M0=sparse(size(Uini,1)/3,(params.nlayers-1)*prod(Nnodes));
                Mi=[Mo,M0];
                M=blkdiag(Mi,blkdiag(Mi,Mi));
            end
            L=[Xo,Yo,Zo];
            F=1+0*Xo;
            a=L\F;
            n=sign(a(3))*a(1:3)./norm(a(1:3));
            
            if numel(params.thickness)==1
                if strcmp(params.type,'solid-from-topo')
                    Xoo=mean(Xo);Yoo=mean(Yo);Zoo=mean(Zo);
                    z=(Xo-Xoo)*n(1)+(Yo-Yoo)*n(2)+(Zo-Zoo)*n(3);
                    zmin=min(z);
                    to=0:-1/(params.nlayers-1):-1;
                    t=(2*(z-zmin)+params.thickness)*to;
                    warning('check the expression');
                else
                    to=params.thickness*(0:-1/(params.nlayers-1):-1);
                    t=repmat(to,length(Xo),1);
                    
                end
            elseif numel(params.thickness)==length(Xo)
                to=0:-1/(params.nlayers-1):-1;
                t=params.thickness*to;
            end
            xo=n(1)*t(:)+repmat(Xo,params.nlayers,1);
            yo=n(2)*t(:)+repmat(Yo,params.nlayers,1);
            zo=n(3)*t(:)+repmat(Zo,params.nlayers,1);
            foundt3=find(elt==3);
            foundq4=find(elt==4);
            conno=conn;
            conn=repmat(0,size(conno,1),8);
            conn(foundt3,1:3)=conno(foundt3,1:3)+prod(Nnodes);
            conn(foundt3,4:6)=conno(foundt3,1:3);
            conn(foundq4,1:4)=conno(foundq4,1:4)+prod(Nnodes);
            conn(foundq4,5:8)=conno(foundq4,1:4);
            elt=2*elt;
            if params.nlayers>2
                connad=repmat(0:(params.nlayers-2),length(elt),1);
                connad=repmat(connad(:)*prod(Nnodes),1,size(conn,2));
            else
                connad=0;
            end
            elt=repmat(elt,params.nlayers-1,1);
            conn=repmat(conn,params.nlayers-1,1);
            conn=conn+connad;
            selected=repmat(selected,params.nlayers,1);
            Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
            Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
            Nnodes=[length(xo),1,1];
            Nelems=[size(conn,1),1,1];
            ns=ones(1,3);
            save(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');
            if nargin>1
                phi=[ 1+0*Xo,0*Yo,0*Zo,0*Xo,Zo,-Yo;...
                    0*Xo,1+0*Yo,0*Zo,-Zo,0*Yo,Xo;...
                    0*Xo,0*Yo,1+0*Zo,Yo,-Xo,0*Zo];
                rbt=phi\Uini;
                Urbt=phi*rbt;
                Uxini=Uini((1:size(Uini,1)/3),:);
                Uyini=Uini(size(Uini,1)/3+(1:size(Uini,1)/3),:);
                Uzini=Uini(2*size(Uini,1)/3+(1:size(Uini,1)/3),:);
                Uxrbt=Urbt((1:size(Urbt,1)/3),:);
                Uyrbt=Urbt(size(Urbt,1)/3+(1:size(Urbt,1)/3),:);
                Uzrbt=Urbt(2*size(Urbt,1)/3+(1:size(Urbt,1)/3),:);
                
                
                Un=(Uxini-Uxrbt)*n(1)+(Uyini-Uyrbt)*n(2)+(Uzini-Uzrbt)*n(3);
                Un=Un-repmat(mean(Un,1),size(Un,1),1);
                Uxini3D=repmat(Uxini-n(1)*Un,params.nlayers,1);
                Uyini3D=repmat(Uyini-n(2)*Un,params.nlayers,1);
                Uzini3D=repmat(Uzini-n(3)*Un,params.nlayers,1);
                to=1+(0:-1/(params.nlayers-1):-1)*2;
                %                to=1:-1/(params.nlayers-1):0;
                for iim=1:size(Uini,2)
                    Un3D=Un(:,iim)*to;
                    Uxini3D(:,iim)=Uxini3D(:,iim)+Un3D(:)*n(1);
                    Uyini3D(:,iim)=Uyini3D(:,iim)+Un3D(:)*n(2);
                    Uzini3D(:,iim)=Uzini3D(:,iim)+Un3D(:)*n(3);
                end
                
                
                Uini3D=[Uxini3D;Uyini3D;Uzini3D];
                
            end
    end
end

end

