function run_williams_curved(id,zone,cparam,filres)

check=cparam.check;
Step=cparam.steps;
Rmax=cparam.radius;% taille de la zone d'extraction
mask_radius=cparam.mask_radius;% taille de la zone a exclure au niveau de la pointe
mask_width=cparam.mask_width;% taille de la zone a exclure le long des faces de la fissure
km_indices=-3:7;
modes=[1,2];
E=cparam.young;% module de young
nu=cparam.nu;% coefficient de poisson
mu=0.5*E/(1+nu);
kappa=(3-4*nu);% kolossov constant en deformation plane (echantillon massif)
pix2m=cparam.pix2m;% correspondance pixel/metre
maxiter=500;
filk=sprintf('%s-crack-%02d-sif.res',strrep(filres,'.res',''),id);
load(filres,'-mat','U','xo','yo','param','Nnodes','Nelems','model','conn','elt');
etype=min(elt);
model.zone=zone;
U=U(:,Step);
U=U(1:prod(Nnodes),:)+1i*U(prod(Nnodes)+(1:prod(Nnodes)),:);
sizeim=[param.roi(2)-param.roi(1),param.roi(4)-param.roi(3)]+1;
xoo=xo;yoo=yo;
zo=xo+param.roi(1)-1+1i*(yo+param.roi(3)-1);
conno=conn;
elto=elt;
conn=cell(1,numel(Step));
elt=cell(1,numel(Step));
xo=cell(1,numel(Step));
yo=cell(1,numel(Step));
Ut=cell(numel(Step),1);
Uw=cell(numel(Step),1);
Urbt=zeros(1,numel(Step));
checka=zeros(numel(Step),1);

xyc=zone{2};
zcini=xyc*[1;1i];
tend=xyc(end,:)-xyc(end-1,:);
xyc(end,:)=xyc(end,:)+1e0*tend;
a=0;
scal=2*mu*sqrt(2*pi);
scalamp=scal*(pix2m.^(1-(km_indices)*.5));

Uo=U;
Iterationfin=ones(1,size(U,2));
xytips=ones(size(U,2),1);
Tfin=ones(1,size(U,2));
Bfin=ones(1,size(U,2));
K1fin=ones(1,size(U,2));
K2fin=zeros(1,size(U,2));
Afin=ones(1,size(U,2));
Widthfin=ones(1,size(U,2));

xc=xyc(1,1);yc=xyc(1,2);
xf=xyc(2,1);yf=xyc(2,2);
zci=(xc+1i*yc);
zf=(xf+1i*yf);
tt=zci-zf;
tt=tt/abs(tt);
nn=tt*exp(1i*pi/2);
thetap=angle(tt);
xcp=xc;
ycp=yc;
zc=xyc*[1;1i];

dcum=abs(diff(zc-zc(1)));
if check
    hf= figure;
end
decc=0;
for iim=size(U,2):-1:1
    disp(sprintf('Step %2d / %2d',Step(iim),Step(size(U,2))));
    iteration=0;
    decx=0;dda=Inf;
    while ((iteration<maxiter) && (((abs(decx)>1)&&(abs(dda)>0.1))||iteration<2))
        iteration=iteration+1;
%    if Step(iim)==14
%       keyboard
%    end
        if abs(decx)>max(sizeim)
            if iteration==1
                dda=0;
            else
                dda=-1.9*abs(dda)*sign(decx);
            end
        elseif abs(decx)>max(model.mesh_size)
            dda=2*mean(model.mesh_size)*sign(decx);
        else
            dda=decx;
        end
        dda=0.5/(1+(iteration>40))/(1+(iteration>50))/(1+(iteration>100))/(1+(iteration>200))/(1+(iteration>20))/(1+(iteration>30))*dda;
        decc=decc+dda;
        if (-(decc)>dcum(1))&(size(xyc,1)>2)
%             keyboard
            a=a-(decc-dda+dcum(1));
            xyc=xyc(2:end,:);
            zc=xyc*[1;1i];
            dcum=abs(diff(zc-zc(1)));
            decc=0;
        else
            a=a+dda;
            xyc(1,1)=xyc(1,1)+dda*real(tt);
            xyc(1,2)=xyc(1,2)+dda*imag(tt);
            zc=xyc*[1;1i];
        end
        
        xc=xyc(1,1);yc=xyc(1,2);
        xf=xyc(2,1);yf=xyc(2,2);
        zci=(xc+1i*yc);
        zf=(xf+1i*yf);
        tt=zci-zf;
        tt=tt/abs(tt);
        nn=tt*exp(1i*pi/2);
        thetap=angle(tt);
        
        Uef=U(:,iim)*exp(-1i*thetap);
        
        
        
        [crack,front,~,~]=GetSignedDistanceToCrack(xyc,zo,0);
        %front=real((zo-zc).*conj(tt));
        %crack=real((zo-zc).*conj(nn));
        %            front=front-0.5*decx;
        z=front+1i*crack;
        dist=max(abs(z),1);
        angl=angle(z);
%                     hf1=figure;
%                     scatter(real(zo),imag(zo),[],front)
%                     hold on
%                     plot(zcini,'r-')
%                     plot(zc(1:end-1),'ro')
%                     hf2=figure;
%                     scatter(real(zo),imag(zo),[],crack)
%                     hold on
%                     plot(zcini,'r-')
%                     plot(zc(1:end-1),'ro')
%                      pause
%                     close(hf1)
%                     close(hf2)
        mask=double(~((dist>Rmax)|(dist<mask_radius)|((abs(crack)<mask_width)&(front<0))));
        
        mask=diag(sparse(mask(:)));
        
        icon=0;
        phic=zeros(length(dist),length(modes)*length(km_indices));
        for m=1:length(modes)
            for ii=1:length(km_indices)
                icon=icon+1;
                meq=modes(m);
                heq=km_indices(ii);
                [harm1,amp1]=KMPotFourrier(meq,heq,kappa);
                feq=0*dist(:);
                for kk=1:3
                    feq=feq+amp1(kk)*exp(harm1(kk)*angl(:)*1i/2);
                end
                feq=feq.*(dist(:).^(heq/2));
                %               feq=feq.*exp(1i*thetap);
                
                phic(:,icon)=feq(:);
            end
        end
        
        M=real(phic'*mask*phic);
        F=real(phic'*mask*Uef);
        [LT,UT]=lu(sparse(M));
        Ukm=UT\(LT\F);
        %        Ukm=M\F;
        found=find(km_indices==1);
        K1=scalamp(found)*Ukm(found);
        K2=scalamp(found)*Ukm(found+length(km_indices));
        found=find(km_indices==-1);
        SK1=scalamp(found)*Ukm(found);
        found=find(km_indices==-3);
        SK3=scalamp(found)*Ukm(found);
        decx=-2*SK1/K1/pix2m;
        width=(sqrt(-8*SK3/K1))/pix2m;
        found=find(km_indices==2);
        %        T=[scalamp(found)*Ks(found,:)];
        T=[scalamp(found)*Ukm(found)]*4/sqrt(2*pi);
        found=find(km_indices==3);
        B=[scalamp(found)*Ukm(found)];
        if dda==0&&iteration>1
            decx=0;
        end
        disp(sprintf('At iteration %2d decx %f applied dda %f',iteration,decx,dda));
        
        if check&&~(((iteration<maxiter) && (abs(decx)>1||iteration==0)))
            Xi=real(zo*exp(-1i*thetap));
            Yi=imag(zo*exp(-1i*thetap));
            maskn=diag(mask);
            maskn(maskn==0)=NaN;
            set(0,'CurrentFigure',hf);
            subplot(2,2,1)
            scatter(Xi,Yi,[],imag(maskn.*Uef))
            colorbar
            axis equal
            subplot(2,2,2)
            scatter(Xi,Yi,[],imag(maskn.*(phic*Ukm)));
            colorbar
            axis equal
            subplot(2,2,3:4)
            scatter(Xi,Yi,[],abs(imag(maskn.*Uef)-imag(maskn.*(phic*Ukm))))
            colorbar
            axis equal
            pause%(1)
        end
        
        
    end
   
    decc=decc+(abs(decx)<1)*decx;
    a=a+(abs(decx)<1)*decx;
            xyc(1,1)=xyc(1,1)+(abs(decx)<1)*decx*real(tt);
            xyc(1,2)=xyc(1,2)+(abs(decx)<1)*decx*imag(tt);

    Tfin(iim)=T;
    Bfin(iim)=B;
    K1fin(iim)=K1;
    if any(modes==2)
        K2fin(iim)=K2;
    end
    Afin(iim)=a*pix2m;
    Iterationfin(iim)=iteration;
    Widthfin(iim)=width*pix2m;
    xytips(iim)=xyc(1,:)*[1;1i];
    
    
    checka(iim)=decx;
    nodes=diag(mask)>0;
    [eltn,~]=GetEltsFromNodes(conno,elto,nodes,0);
    [pind,~,j1]=unique(conno(eltn,1:etype));
    conn{iim}=reshape(j1,[numel(eltn),etype]);
    xo{iim}=xoo(pind);
    yo{iim}=yoo(pind);
    elt{iim}=elto(eltn);
    Ut{iim}=[real(Uo(pind,iim));imag(Uo(pind,iim))];
    wil=(phic*Ukm)*exp(1i*thetap);
    Uw{iim}=[real(wil(pind));imag(wil(pind))];
            found=find(km_indices==0);
            Urbt(iim)=exp(1i*thetap)*(phic(1,found)*Ukm(found)+phic(1,found+length(km_indices))*Ukm(found+length(km_indices)));
end
K1=K1fin;
K2=K2fin;
T=Tfin;
B=Bfin;
da=Afin;
rho=Widthfin;
try
if numel(Step)>1
    param.deformed_image=param.deformed_image(:,Step);
else
    if iscell(param.deformed_image)
        param.deformed_image=param.deformed_image{:,Step};
    end
end
catch
end
param.result_file=filk;
param.detect=1;
rint=0;
U=Uw;
if isfield(model,'visu')
    model=rmfield(model,'visu');
end
save(filk,'Urbt','K1','K2','T','B','checka','xytips','Widthfin','da','rho','xo','yo','conn','U','Ut','Nnodes','Nelems','rint','elt','model','param');
%         param.result_file=strrep(filk,'.res','-w.res');
%         U=Uw;
%  save(param.result_file,'K1','K2','T','B','xytips','Widthfin','da','xo','yo','conn','U','Nnodes','Nelems','rint','elt','model','param');
%         param.result_file=strrep(filk,'.res','-dfit.res');
%         for iim=1:size(U)
%             U{iim}=Ut{iim}-Uw{iim};
%         end
%  save(param.result_file,'K1','K2','T','B','xytips','Widthfin','da','xo','yo','conn','U','Nnodes','Nelems','rint','elt','model','param');

end
