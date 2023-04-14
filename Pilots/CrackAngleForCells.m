function Oc=CrackAngleForCells(E,EE,phio)
check=0;
load('cell-hom','model');
matmod=model.material_parameters;
clear model
load('cell-mesh','n','zp','indc')
[~,indz]=sort(angle(zp));
lbox=max(abs(zp))+10;
zp=zp(indz);
nmod=5;
filres='cell-failure';
param.analysis='mechanics';
param.result_file=filres;
param.pixel_size=1e-6;
pix2m=param.pixel_size;
param.roi=[1,2*lbox+10,1,2*lbox+10];

EE=EE*(param.pixel_size);


model.basis='fem';
model.mesh_file='cell-mesh.vtk';
xc=lbox;
yc=lbox;
model.gluing_parameters{1}={'scale',[0;0;0],1};
model.gluing_parameters{2}={'translate',[xc-1;yc-1;0]+1};

mu=1/2*matmod.young/(1+matmod.nu);
lambda=matmod.young*matmod.nu/(1+matmod.nu)/(1-2*matmod.nu);


model.material_model='elastic_homogeneous_isotropic';
model.material_parameters=matmod;
model.vtk_export={'S'};
model.material_parameters=matmod;
model.vtk_export={'S'};

iscale=1;
LoadParameters(param);
LoadParameters(model,nmod);
sizeim=[diff(param.roi(1:2)),diff(param.roi(3:4))]+1;
save(fullfile('TMP','sample0_0'),'sizeim')
save(fullfile('TMP','sample0'),'sizeim')

LoadMask(nmod);   % ???
LoadMeshes(nmod); % generate the meshes
LoadMat(nmod);    % define some materials constants

%%

meshfile=fullfile('TMP',sprintf('%d_mesh_%d',nmod,10*(iscale-1)));
load(meshfile);
phig=CreateFiniteElementBasis(meshfile,sizeim,1,[],'Gauss_points');
zp=[zp(:);zp(1)];
selected=GetSignedDistanceToZone([real(zp),imag(zp)],[min(xo-xc)-10,max(xo-xc)+10,min(yo-xc)-10,max(yo-xc)+10],xo-xc,yo-xc);
selected=selected<-1;

if 0
    zo=(phig*xo-lbox)+1i*(phig*yo-lbox);
    xyc=[0,0;2*lbox*cos(phio+pi),2*lbox*sin(phio+pi)];
    [crack,front]=GetSignedDistanceToCrack(xyc,zo);
    inbox=(((abs(crack)<0.25*lbox)&(front<0)));
    figure
    triplot(conn(inbox,:),xo-lbox,yo-lbox)
    axis equal
    hold on
    plot(xyc(:,1),xyc(:,2),'k')
    face_elts=[];
    inods=conn(inbox,1:3);
    [inods,~,ic]=unique(inods);
    [cn,fn]=GetSignedDistanceToCrack(xyc,xo(inods)-lbox+1i*(yo(inods)-lbox));
    cn=reshape(cn(ic),sum(inbox),3);
    fn=reshape(fn(ic),sum(inbox),3);
    edges=[1,2;2,3;3,1];
    ij=1;
    for ie=1:numel(inbox)
        if inbox(ie)
            inods=conn(ie,:);
            cns=cn(ij,:);
            if any(prod(sign(cns(edges)),2)<0)
                face_elts(end+1)=ie;
            end
            ij=ij+1;
        end
    end
    plot(real(zo(face_elts)),imag(zo(face_elts)),'xr')
    
    conn(face_elts,:)=[];
    elt(face_elts)=[];
    Nelems=[numel(elt),1,1];
    save(meshfile,'conn','elt','Nelems','-append')
    
    zo=(xo-lbox)+1i*(yo-lbox);
    anglp=angle(zp*exp(-1i*phio));
    anglo=angle(zo*exp(-1i*phio));
    selected(anglo>max(anglp))=1;
    selected(anglo<min(anglp))=1;
    selected(abs(anglo)>pi-2*pi/n)=1;
    plot(xo(~selected)-lbox,yo(~selected)-lbox,'ro')
end

CreateGradBasisFunction(1,nmod);%calcul le gradient des fonctions de formes
AssembleMechanicalOperator(1,nmod);%assemblage de la matrice de raideur elements finis
%%
load(fullfile('TMP',sprintf('%d_k_operator_%d.mat',nmod,10*(iscale-1))),'K');
load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'wdetJ');
phigo=CreateFiniteElementBasis(meshfile,sizeim,1,[],'Gauss_points');





Ub=[(xo-xc),(yo-yc),0*xo,0*xo,1/2*(xo-xc).*(xo-xc),1/2*(xo-xc).*(yo-yc),1/2*(yo-yc).*(yo-yc),0*xo                ,0*xo                ,0*xo;...
    0*yo,0*xo,(xo-xc),(yo-yc),0*yo                ,0*yo                ,0*yo                ,1/2*(xo-xc).*(xo-xc),1/2*(xo-xc).*(yo-yc),1/2*(yo-yc).*(yo-yc)];
indi=[find(~selected)];
indi=[indi;indi+prod(Nnodes)];
Up=Ub(indi,:);

C=sparse(indi,1:length(indi),1,2*prod(Nnodes),length(indi));
Fo=zeros(2*prod(Nnodes),size(Up,2));

Kc=[K,C;C',sparse(size(C,2),size(C,2))];
Fc=[Fo;Up];
X=Kc\Fc;
U=X(1:2*prod(Nnodes),:);
U=[U*diag([E;EE]),U*[E;EE]];


for ii=1:size(U,2)
    UpdateInternalState(nmod,U(:,ii),ii);
end
save(filres,'U','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
postproVTK(filres,0,0);
U=U(:,end);

%%
load('cell-mesh','zpp','homo','anglep')
ome=pi+2*pi/n;
np=size(homo,1);
alps=SigValueElasticity(ome);
xoo=xo;
yoo=yo;
gm1=lambda+2*mu; % in plane strain
gm2=lambda;
conno=conn(:,1:3);

%%% Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mode I
H1_I=@(aa)(lambda+3*mu-aa*(lambda+mu))/(lambda+mu)/(1-aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);
H2_I=@(aa)-(lambda+3*mu+aa*(lambda+mu))/(lambda+mu)/(1-aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);


Hs1_I=@(aa)(3-aa)/(1-aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);
Hs2_I=@(aa)(1+aa)/(1-aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);
Hs3_I=@(aa)sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);


ctt_I=@(aa)2*mu*aa*((1+aa)*sin(ome*(1+aa)/2)/((1-aa)*sin(ome*(1-aa)/2))-1);


ctts_I=@(aa)((1+aa)*sin(ome*(1+aa)/2)/((1-aa)*sin(ome*(1-aa)/2))-1);


% Mode II
H1_II=@(aa)(lambda+3*mu-aa*(lambda+mu))/(lambda+mu)/(1+aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);
H2_II=@(aa)(lambda+3*mu+aa*(lambda+mu))/(lambda+mu)/(1+aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);


Hs1_II=@(aa)(3-aa)/(1+aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);
Hs2_II=@(aa)sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);
Hs3_II=@(aa)-(1-aa)/(1+aa)*sin(ome*(1+aa)/2)/sin(ome*(1-aa)/2);


ctt_II=@(aa)2*mu*aa*(1-(1-aa)*sin(ome*(1+aa)/2)/((1+aa)*sin(ome*(1-aa)/2)));


ctts_II=@(aa)(1-(1-aa)*sin(ome*(1+aa)/2)/((1+aa)*sin(ome*(1-aa)/2)));

%%% Displacement and stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mode I

UU_I=@(aa,cctt,cH1,cH2,tt,r,aome)r.^(aa).*(cos((1+aa)*(tt-aome))+cH1*cos((1-aa)*(tt-aome)))/cctt.*cos((tt-aome))...
    -r.^(aa).*(-sin((1+aa)*(tt-aome))+cH2*sin((1-aa)*(tt-aome)))/cctt.*sin((tt-aome));

VV_I=@(aa,cctt,cH1,cH2,tt,r,aome)r.^(aa).*(cos((1+aa)*(tt-aome))+cH1*cos((1-aa)*(tt-aome)))/cctt.*sin((tt-aome))...
    +(r).^(aa).*(-sin((1+aa)*((tt-aome)))+cH2*sin((1-aa)*((tt-aome))))/cctt.*cos((tt-aome));

Ur_I=@(aa,cctt,cH1,cH2,tt,r,aome)(r).^(aa).*(cos((1+aa)*((tt-aome)))+cH1*cos((1-aa)*((tt-aome))))/cctt;

Ut_I=@(aa,cctt,cH1,cH2,tt,r,aome) (r).^(aa).*(-sin((1+aa)*((tt-aome)))+cH2*sin((1-aa)*((tt-aome))))/cctt;

Srr_I=@(aa,cctts,cHs1,cHs2,cHs3,tt,r,aome)(r).^(aa-1).*(cos((1+aa)*((tt-aome)))+cHs1*cos((1-aa)*((tt-aome))))/cctts;

Stt_I=@(aa,cctts,Hs1,Hs2,Hs3,thetha,r,aome)(r).^(aa-1).*(-cos((1+aa)*(thetha+aome))+Hs2*cos((1-aa)*(thetha+aome)))/cctts;

Srt_I=@(aa,cctts,cHs1,cHs2,cHs3,tt,r,aome)(r).^(aa-1).*(-sin((1+aa)*((tt-aome)))+cHs3*sin((1-aa)*((tt-aome))))/cctts;

% Mode II


UU_II=@(aa,cctt,cH1,cH2,tt,r,aome)r.^(aa).*(sin((1+aa)*(tt-aome))+cH1*sin((1-aa)*(tt-aome)))/cctt.*cos((tt-aome))...
    -r.^(aa).*(-sin((1+aa)*(tt-aome))+cH2*sin((1-aa)*(tt-aome)))/cctt.*sin((tt-aome));

VV_II=@(aa,cctt,cH1,cH2,tt,r,aome)r.^(aa).*(cos((1+aa)*(tt-aome))+cH1*cos((1-aa)*(tt-aome)))/cctt.*sin((tt-aome))...
    +(r).^(aa).*(-sin((1+aa)*((tt-aome)))+cH2*sin((1-aa)*((tt-aome))))/cctt.*cos((tt-aome));

Ur_II=@(aa,cctt,cH1,cH2,tt,r,aome)(r).^(aa).*(cos((1+aa)*((tt-aome)))+cH1*cos((1-aa)*((tt-aome))))/cctt;

Ut_II=@(aa,cctt,cH1,cH2,tt,r,aome) (r).^(aa).*(-sin((1+aa)*((tt-aome)))+cH2*sin((1-aa)*((tt-aome))))/cctt;

Srr_II=@(aa,cctts,cHs1,cHs2,cHs3,tt,r,aome)(r).^(aa-1).*(cos((1+aa)*((tt-aome)))+cHs1*cos((1-aa)*((tt-aome))))/cctts;

Stt_II=@(aa,cctts,Hs1,Hs2,Hs3,thetha,r,aome)(r).^(aa-1).*(-sin((1+aa)*(thetha+aome))+Hs2*sin((1-aa)*(thetha+aome)))/cctts;

Srt_II=@(aa,cctts,cHs1,cHs2,cHs3,tt,r,aome)(r).^(aa-1).*(-sin((1+aa)*((tt-aome)))+cHs3*sin((1-aa)*((tt-aome))))/cctts;

anglc=zeros(np,1);

if 1
    crit=3;
    fics=zeros(np,3);
    for ic=1:np
%        if ic==5,check=1;end
        xc=xoo(indc(:,ic));
        yc=yoo(indc(:,ic));
        Uz=U(1:numel(xoo),end)+1i*U(numel(xoo)+(1:numel(xoo)),end);
        
        Px=fminsearch(@(P)(norm((P(1)^2-abs(xc-P(2)+1i*(yc-P(3))).^2))),[abs(max(xc)-min(xc)),mean(xc),mean(yc)]);
        if check
        figure
        plot(xc,yc,'ro')
        hold on
        plot(Px(2),Px(3),'rx')
        plot(Px(2)+Px(1)*cos(0:2*pi/100:2*pi),Px(3)+Px(1)*sin(0:2*pi/100:2*pi),'b-')

        end
        anglc(ic)=angle(Px(2)-lbox+1i*(Px(3)-lbox));
        [~,idmin]=min(abs(xoo-Px(2)+1i*(yoo-Px(3))));
        Uz=(Uz-Uz(idmin))*exp(-1i*anglep((ic)));
        zc=((xoo-Px(2))+1i*(yoo-Px(3)))*exp(-1i*anglep((ic)));
        zg=phigo*zc;
        inside=abs(zg)<Px(1);
        conn=conno(inside,:);
        xo=real(zc);
        yo=imag(zc);
        zo=1;
        elt=3*ones(size(conn,1),1);
        ng=0;ns=[1,1];rint=0;
        Nnodes=[numel(xo),1,1];
        Nelems=[size(conn,1),1,1];
        save('contour','xo','yo','zo','conn','elt','Nnodes','Nelems','rint','ng')
        wg=GetWeigthDetJ('contour',[1,1],1,'Gauss_points');
        phig=CreateFiniteElementBasis('contour',[1,1],1,[],'Gauss_points');
        [dphigdx,dphigdy]=CreateGradFiniteElementBasis('contour',[1,1],1,[],'Gauss_points');
        zg=phig*zc;
        a1=0;
    aome=0;
    thp=angle(zg);
        rp=abs(zg);
        
        uxp=phig*real(Uz);
        uyp=phig*imag(Uz);
        dfac=0.25;
        rmax=max(Px(1));
        dqpdr=(abs(zg)>dfac*rmax)/rmax/(1-dfac);
        
        urp=uxp.*cos(thp-aome)+uyp.*sin(thp-aome);
        utp=uxp.*(-sin(thp-aome))+uyp.*cos(thp-aome);
        
        exxp=dphigdx*real(Uz);
        eyyp=dphigdy*imag(Uz);
        exyp=(dphigdy*real(Uz)+dphigdx*imag(Uz))/2;
        
        
        tre=exxp+eyyp;
        sxxp=lambda*tre+2*mu*exxp;
        syyp=lambda*tre+2*mu*eyyp;
        sxyp=2*mu*exyp;
        srrp=(cos(thp+aome).*sxxp+sin(thp+aome).*sxyp).*cos(thp+aome)+(cos(thp+aome).*sxyp+sin(thp+aome).*syyp).*sin(thp+aome);
        srtp=(cos(thp+aome).*sxxp+sin(thp+aome).*sxyp).*(-sin(thp+aome))+(cos(thp+aome).*sxyp+sin(thp+aome).*syyp).*cos(thp+aome);
        for ia=1:sum(alps<1)
            
            alp=alps(ia);
            
                alp1=alps(1);
            if abs(a1)>0
                

                cctt=ctt_I(alp1);
                cH1=H1_I(alp1);
                cH2=H2_I(alp1);
                cHs1=Hs1_I(alp1);
                cHs2=Hs2_I(alp1);
                cHs3=Hs3_I(alp1);
                cctts=ctts_I(alp1);

                uxp=uxp-a1*UU_I(alp1,cctt,cH1,cH2,thp,Px(1),aome);
                uyp=uyp-a1*VV_I(alp1,cctt,cH1,cH2,thp,Px(1),aome);      


                cSrr=Srr_I(alp1,cctts,cHs1,cHs2,cHs3,thp,Px(1),aome);
                cSrt=Srt_I(alp1,cctts,cHs1,cHs2,cHs3,thp,Px(1),aome);
                cStt=Stt_I(alp1,cctts,cHs1,cHs2,cHs3,thp,Px(1),aome);


                sxxp=sxxp-a1*((cos(thp).*cSrr-sin(thp).*cSrt).*cos(thp)-(cos(thp).*cSrt-sin(thp).*cStt).*sin(thp));
                syyp=syyp-a1*((sin(thp).*cSrr+cos(thp).*cSrt).*cos(thp)-(sin(thp).*cSrt+cos(thp).*cStt).*sin(thp));
                sxyp=sxyp-a1*((sin(thp).*cSrr+cos(thp).*cSrt).*sin(thp)+(sin(thp).*cSrt+cos(thp).*cStt).*cos(thp));

            end
            switch ia
            case 1 
                cctt=ctt_I(alp);
                cH1=H1_I(alp);
                cH2=H2_I(alp);

                cHs1=Hs1_I(alp);
                cHs2=Hs2_I(alp);
                cHs3=Hs3_I(alp);
                cctts=ctts_I(alp);

                ncH1=H1_I(-alp);
                ncH2=H2_I(-alp);
                ncctt=ctt_I(-alp);

                ncctts=ctts_I(-alp);
                ncHs1=Hs1_I(-alp);
                ncHs2=Hs2_I(-alp);
                ncHs3=Hs3_I(-alp);

                usx=UU_I(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);
                usy=VV_I(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);

        
                usx=((usx));
                usy=((usy));
                nx=((cos(thp)));
                ny=((sin(thp)));


        %%%%

                 ssrr=Srr_I(-alp,ncctts,ncHs1,ncHs2,ncHs3,thp,Px(1),aome);
        %        
                 ssrt=Srt_I(-alp,ncctts,ncHs1,ncHs2,ncHs3,thp,Px(1),aome);

                urs=Ur_I(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);

                uts=Ut_I(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);

                I3=dqpdr.*(urs.*srrp+uts.*srtp);
                I4=dqpdr.*(urp.*ssrr+utp.*ssrt);


                icontour=sum(wg*I3,1)-sum(wg*I4,1);


                srr=Srr_I(alp,cctts,cHs1,cHs2,cHs3,thp,Px(1),aome);
                srt=Srt_I(alp,cctts,cHs1,cHs2,cHs3,thp,Px(1),aome);


                nUUr=Ur_I(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);
                UUr=Ur_I(alp,cctt,cH1,cH2,thp,Px(1),aome);

                nVVt=Ut_I(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);
                VVt=Ut_I(alp,cctt,cH1,cH2,thp,Px(1),aome);

                I5=dqpdr.*(urs.*srr+uts.*srt);
                I6=dqpdr.*(UUr.*ssrr+VVt.*ssrt);


                Aalp=sum(wg*I5,1)-sum(wg*I6,1);

                fics(ic,ia)=icontour/Aalp;

                a1=fics(ic,ia);
                if check
                figure
                plot(thp-aome,uxp,'rx')
                hold on
                plot(thp-aome,uyp,'bx')
                title(sprintf('mode %d uef',ia))
                figure


                plot(thp-aome,a1*UU_I(alp,cctt,cH1,cH2,thp,rp,aome),'rx')
                hold on
                plot(thp-aome,a1*VV_I(alp,cctt,cH1,cH2,thp,rp,aome),'bx')
                title(sprintf('mode %d uana',ia))

                figure
                 plot(thp-aome,uxp,'rx')
                 hold on
                 plot(thp-aome,a1*UU_I(alp,cctt,cH1,cH2,thp,rp,aome),'gx')
                 title(sprintf('u fem vs u analytique mode %d uana',ia))
                 
                 
                figure
                  plot(thp-aome,uyp,'rx')
                 hold on
                 plot(thp-aome,a1*VV_I(alp,cctt,cH1,cH2,thp,rp,aome),'bx')
                 title(sprintf('v fem vs v analytique mode %d uana',ia))
                 
            end
            case 2
                
                
                cctt=ctt_II(alp);
                cH1=H1_II(alp);
                cH2=H2_II(alp);

                cHs1=Hs1_II(alp);
                cHs2=Hs2_II(alp);
                cHs3=Hs3_II(alp);
                cctts=ctts_II(alp);

                ncH1=H1_II(-alp);
                ncH2=H2_II(-alp);
                ncctt=ctt_II(-alp);

                ncctts=ctts_II(-alp);
                ncHs1=Hs1_II(-alp);
                ncHs2=Hs2_II(-alp);
                ncHs3=Hs3_II(-alp);



%                 uxp=uxp+mean(UU_II(alp,cctt,cH1,cH2,thp,Px(1),aome)-uxp);
%                 uyp=uyp+mean(VV_II(alp,cctt,cH1,cH2,thp,Px(1),aome)-uyp);

                usx=UU_II(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);
                usy=VV_II(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);

       
                usx=((usx));
                usy=((usy));
                nx=((cos(thp)));
                ny=((sin(thp)));


        %%%%

                ssrr=Srr_II(-alp,ncctts,ncHs1,ncHs2,ncHs3,thp,Px(1),aome);
                ssrt=Srt_II(-alp,ncctts,ncHs1,ncHs2,ncHs3,thp,Px(1),aome);

                urs=Ur_II(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);

                uts=Ut_II(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);

                I3=dqpdr.*(urs.*srrp+uts.*srtp);
                I4=dqpdr.*(urp.*ssrr+utp.*ssrt);


                icontour=sum(wg*I3,1)-sum(wg*I4,1);


                icontour=sum(wg*I3,1)-sum(wg*I4,1);


                srr=Srr_II(alp,cctts,cHs1,cHs2,cHs3,thp,Px(1),aome);
                srt=Srt_II(alp,cctts,cHs1,cHs2,cHs3,thp,Px(1),aome);


                nUUr=Ur_II(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);
                UUr=Ur_II(alp,cctt,cH1,cH2,thp,Px(1),aome);

                nVVt=Ut_II(-alp,ncctt,ncH1,ncH2,thp,Px(1),aome);
                VVt=Ut_II(alp,cctt,cH1,cH2,thp,Px(1),aome);

                I5=dqpdr.*(urs.*srr+uts.*srt);
                I6=dqpdr.*(UUr.*ssrr+VVt.*ssrt);


                Aalp=sum(wg*I5,1)-sum(wg*I6,1);

                fics(ic,ia)=icontour/Aalp;

                a1=fics(ic,ia);
                if check
                 figure
                plot(thp-aome,uxp,'rx')
                hold on
                plot(thp-aome,uyp,'bx')
                title(sprintf('mode %d uef',ia))
                figure


                plot(thp-aome,fics(ic,1)*UU_I(alp1,ctt_I(alp1),H1_I(alp1),H2_I(alp1),thp,rp,aome)+a1*UU_II(alp,cctt,cH1,cH2,thp,rp,aome),'rx')
                hold on
                plot(thp-aome,fics(ic,1)*VV_I(alp1,ctt_I(alp1),H1_I(alp1),H2_I(alp1),thp,rp,aome)+a1*VV_II(alp,cctt,cH1,cH2,thp,rp,aome),'bx')
                title(sprintf('mode %d uana',ia))

                figure
                 plot(thp-aome,uxp,'rx')
                 hold on
                 plot(thp-aome,fics(ic,1)*UU_I(alp1,ctt_I(alp1),H1_I(alp1),H2_I(alp1),thp,rp,aome)+a1*UU_II(alp,cctt,cH1,cH2,thp,rp,aome),'gx')
                 title(sprintf('u fem vs u analytique mode %d uana',ia))
                 
                 
                figure
                  plot(thp-aome,uyp,'rx')
                 hold on
                 plot(thp-aome,fics(ic,1)*VV_I(alp1,ctt_I(alp1),H1_I(alp1),H2_I(alp1),thp,rp,aome)+a1*VV_II(alp,cctt,cH1,cH2,thp,rp,aome),'bx')
                 title(sprintf('v fem vs v analytique mode %d uana',ia))
                end
            end % of switch
            
%        
        end

%        keyboard
    end
    fics(:,3)=(fics(:,1)>=0).*sqrt(fics(:,1).^2+fics(:,2).^2);
    if any(isnan(fics(:)))
%        keyboard
    else
        close all
        pause(0.1)
    end
else
    xo=xoo;
    yo=yoo;
    fics=zeros(np,4);
    load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
    load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'epsxy');
    load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
    [sxx,syy,sxy,sxz,syz,szz]=ComputeStress(epsxx*U,epsyy*U,epsxy*U,nmod);
    S=[sxx,syy,szz,sxy,syz,sxz];
    J2=sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2 ...
        -S(:,1).*S(:,3)-S(:,2).*S(:,3)-S(:,1).*S(:,2) ...
        +3*(S(:,4).^2+S(:,5).^2+S(:,6).^2));
    %                         figure
    % scatter(phig*xo,phig*yo,[],J2)
    phig=CreateFiniteElementBasis(meshfile,sizeim,1,[],'Gauss_points');
    %  figure
    %  hold on
    for ic=1:np
        
        xc=xoo(indc(:,ic));
        yc=yoo(indc(:,ic));
        
        Px=fminsearch(@(P)(norm((P(1)^2-abs(xc-P(2)+1i*(yc-P(3))).^2))),[abs(max(xc)-min(xc)),mean(xc),mean(yc)]);
        anglc(ic)=angle(Px(2)-lbox+1i*(Px(3)-lbox));
        zc=(phig*xo-Px(2)+1i*(phig*yo-Px(3)));
        mask=abs(zc)<Px(1);
        nn=xc(end)-xc(1)+1i*(yc(end)-yc(1));
        nn=nn/abs(nn);
        %    nn=exp(1i*anglc(ic))*1i;
        nn=[real(nn);imag(nn)];
        npdf=exp(-0.5*(abs(zc)/(Px(1)/1)).^2)*0+1;
        %      figure
        % %     scatter(phig(mask,:)*xo,phig(mask,:)*yo,[],syy(mask))
        
        %    scatter(phig(mask,:)*xo,phig(mask,:)*yo,[],J2(mask))
        %    plot(real(zpp(ic))+lbox,imag(zpp(ic))+lbox,'or','MarkerSize',20,'LineWidth',2)
        
        %      scatter(phig(mask,:)*xo,phig(mask,:)*yo,[],npdf(mask))
        %      axis equal
        %      hold on
        %         plot(xc,yc,'ro')
        %         hold on
        %         plot(Px(2),Px(3),'rx','MarkerSize',20)
        %         plot(Px(2)+Px(1)*cos(0:2*pi/100:2*pi),Px(3)+Px(1)*sin(0:2*pi/100:2*pi),'b-')
        % quiver(Px(2),Px(3),Px(1)*nn(1),Px(1)*nn(2),'r','LineWidth',2)
        mask= diag(sparse(double(mask.*npdf)));
        Sc=sum(diag(mask*wdetJ));
        Sm=zeros(2);
        SM=zeros(2);
        Sm(1,1)=sum(mask*(wdetJ*sxx))/Sc;
        Sm(2,2)=sum(mask*(wdetJ*syy))/Sc;
        Sm(1,2)=sum(mask*(wdetJ*sxy))/Sc;
        Sm(2,1)=Sm(1,2);
        eS=0;%max(abs(diff(eigs(Sm))),max(abs(eigs(Sm))));
        
        SM(1,1)=sum((wdetJ*sxx))/Sc;
        SM(2,2)=sum((wdetJ*syy))/Sc;
        SM(1,2)=sum((wdetJ*sxy))/Sc;
        SM(2,1)=SM(1,2);
        fics(ic,1)=nn'*(Sm*nn);
        fics(ic,2)=Sm(1,1)+Sm(2,2);
        fics(ic,3)=eS;
        fics(ic,4)=sum(mask*(wdetJ*J2))/Sc;
        %    pause
        %    keyboard
    end
    crit=2;
end
anglc
fics
if 1
    %[~,id]=max(fics.*(abs(angle(zp(:)*exp(-1i*phio)))>=pi/2));
    %angl=-angle(-zp(id));
    %[~,id]=max((fics(:,1)).*(abs(angle(zp(:)*exp(-1i*phio)))<=pi/2))
    %angl=angle(zp(id));
    [~,id]=max(fics(:,crit))
    %keyboard
    [~,sid1]=sort(fics(mod(homo(id,:)-1,np)+1,crit))
    id1=homo(id,sid1(1))
    angl=angle(diff(zpp([id,id1])))
    figure
    triplot(conno,xoo-lbox,yoo-lbox)
    hold on
    plot(zpp([id,id1]),'r')
    plot(real(zpp([id])),imag(zpp([id])),'or','MarkerSize',20,'LineWidth',2)
    Oc=mod(angl-phio,pi)
    
    
    %keyboard
    if numel(sid1)>1&&abs(abs(Oc)-pi/2)<0.1*pi
        id1=homo(id,sid1(2))
        angl=angle(diff(zpp([id,id1])))
        Oc=mod(angl-phio,pi)
    end
    if abs(Oc)>pi/2
        Oc=Oc-sign(Oc)*pi
    end
    
    % [~,id]=max(fics(:,1))
    % angl=anglc(id);
    % if abs(angl)>pi/2
    %     angl=sign(angl)*pi-angl;
    % end
    
    
else
    Oc=sum(fics.*((angle(zp(:)*exp(-1i*phio)))).*(abs(angle(zp(:)*exp(-1i*phio)))<=pi/2))/sum(fics.*(abs(angle(zp(:)*exp(-1i*phio)))<=pi/2))
end
%keyboard
