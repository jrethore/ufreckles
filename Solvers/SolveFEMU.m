function [U,P,Mt,Ft]=SolveFEMU(Uref,Pini,nmod)
ttic=cputime;
clear AdditionalConstrains
if nargout<3,Fe=0;end
check=0;
load(fullfile('TMP','params'),'param');
param0=param;
roi=param0.roi;
maxiter=param0.iter_max;
conv=param0.convergance_limit;
if isfield(param0,'search_convergance_limit')
    conv=param0.search_convergance_limit;
end
nim=size(Uref,2);
if isfield(param0,'under_relaxation')
    alpha=param0.under_relaxation;
else
    alpha=1;
end
if isfield(param0,'weighting')
    weighting=param0.weighting;
else
    weighting='uniform';
end
if strcmp(weighting,'user')
    w=param0.weight;
else
    w=GetWeight(weighting,nim);
end
if isfield(param0,'coupling_weight')
    l=param0.coupling_weight;
else
    l=1;
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
disp(sprintf('Starting resolution ...'));
iscale=1;
load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'phix','sizeim');
load(fullfile('TMP',sprintf('%d_phiy_%d',nmod,10*(iscale-1))),'phiy');
load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
load(fullfile('TMP',sprintf('%d_operator_%d',nmod,iscale-1)),'M');
Mc=M;
if length(Uref)==prod(sizeim)
    onpixel=1;
    uref=Uref;
else
    uref=Uref;
    onpixel=0;
end
[U,Fe]=UpdateUField(nmod,Pini);
if onpixel
    umap=mask*(phix*U+i*(phiy*U));
else
    umap=U;
end
if check
    dmap=mask*(phix*U(:,nim)+i*(phiy*U(:,nim)));
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,0)),'Nnodes','xo','yo');
    xn=xo;
    xn(length(xn))=xo(length(xo))-1;
    yn=yo;
    yn(length(yn))=yo(length(yo))-1;
    tmap=reshape(real(dmap),sizeim);
    unmap=tmap(xn,yn)-mean(tmap(:));
    tmap=reshape(imag(dmap),sizeim);
    vnmap=tmap(xn,yn)-mean(tmap(:));
    fac=0.2*max(sizeim)/max(abs(unmap(:)+i*vnmap(:)));
    limd=[min(abs(unmap(:)+i*vnmap(:))),max(abs(unmap(:)+i*vnmap(:)))];
    DeformedMesh(xn,yn,unmap,vnmap,limd,fac,['solve-femu'],0,1);
    pause
end
tmp=uref-umap;
Foini=sum(diag(real(tmp'*tmp)).*w);
Fo=Foini;
No=sum(diag(real(uref'*uref)).*w);
disp(sprintf('Relative gap =%6.2f %%',100*sqrt(Fo/No)));
    if l<1
    ferror=0;
        for iim=1:nim
            ferror=ferror+w(iim)*Fe(:,iim)'*Fe(:,iim);
        end
    else
        ferror=1;
    end
    ferrorini=ferror;
res=1;
ii=1;
P=Pini;
        nU0=norm(P);
tiko=0.e-1;
while ( res>conv && ii< maxiter)
    Pold=P;
    Ju=repmat(0,[size(uref,1),length(P),nim]);
    Jf=repmat(0,[size(Fe,1),length(P),nim]);
    for ip=1:length(P)
        dPi=0*P;
        dPi(ip)=0.01*Pini(ip);
        [U,Fi]=UpdateUField(nmod,P+dPi);
            for iim=1:nim
        if onpixel
            Ju(:,ip,iim)=mask*(phix*U(:,iim)+i*(phiy*U(:,iim))-umap(:,iim))/dPi(ip);
        else
            Ju(:,ip,iim)=(U(:,iim)-umap(:,iim))/dPi(ip);
        end
        Jf(:,ip,iim)=-(Fi(:,iim)-Fe(:,iim))/dPi(ip);
    end
    end
    M=sparse(length(P),length(P));
    for iim=1:nim
        M=M+w(iim)*real(Ju(:,:,iim)'*Ju(:,:,iim));
    end
    Mf=sparse(length(P),length(P));
    for iim=1:nim
        Mf=Mf+w(iim)*Jf(:,:,iim)'*Jf(:,:,iim);
    end
    jj=1;
    res2=Inf;conv2=conv;maxiter2=2;%maxiter;

    while ( res2>conv2 && jj< maxiter2)

    
    
    F=zeros(length(P),1);
    if l>0
        for iim=1:nim
            F=F+w(iim)*real(Ju(:,:,iim)'*(uref(:,iim)-umap(:,iim)));
        end
%         figure
%      load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,0)),'Nnodes','xo','yo');
%        scatter(xo,yo,10+0*xo,abs(uref(length(xo)+(1:length(xo)),iim)-umap(length(xo)+(1:length(xo)),iim)))
%         pause
    end
    Ff=zeros(length(P),1);
    ferror=0;
    if l<1
        for iim=1:nim
            ferror=ferror+w(iim)*Fe(:,iim)'*Fe(:,iim);
            Ff=Ff+w(iim)*Jf(:,:,iim)'*Fe(:,iim);
        end
    end
    if (ii==1)
        if l>0
            normRu=norm(F);
        else
            normRu=1;
        end
        if l<1
            normRf=norm(Ff);
        else
            normRf=1;
        end
    ertot=l*Fo/normRu+(1-l)*ferror/normRf;
        
    end
    Mt=l*M/normRu+(1-l)*Mf/normRf;
    Ft=l*F/normRu+(1-l)*Ff/normRf;
    dP=alpha*((Mt+tiko*diag(diag(Mt)))\Ft);
    [dP]=AdditionalConstrains(Mt+tiko*diag(diag(Mt)),Ft,P,dP);
    [dP]=AdditionalConstrains(Mt+tiko*diag(diag(Mt)),Ft,P,dP);

close all
pause(0.1)


    [U,Fe]=UpdateUField(nmod,P+dP);
    

    
        if onpixel
        umap=mask*(phix*U+i*(phiy*U));
    else
        umap=U;
    end
    tmp=umap-uref;
    Fo=sum(diag(real(tmp'*tmp)).*w);
    ferror=0;
    if l<1
        for iim=1:nim
            ferror=ferror+w(iim)*Fe(:,iim)'*Fe(:,iim);
        end
    end

ertotnew=l*Fo/normRu+(1-l)*ferror/normRf;

    disp(sprintf('Update # %d At iteration # %d Tiko %f',ii,jj,tiko));
    disp(sprintf('Relative gap |Ue-Ui|/|Ue| =%6.2f %% |dF|/|dFo|  %6.2f %%',100*sqrt(Fo/No),sqrt(ferror/ferrorini)*100));
        res2=norm(dP)/nU0;
        res2=max(abs(dP./P));
        disp(sprintf('|dP|/|Po|=%f',res2));

    if ertotnew>ertot
    tiko=tiko*10;
    else
        tiko=tiko/10;
    end
    tiko=max(tiko,1.e-3);
    ertot=ertotnew;
    P=P+dP
    if check
        dmap=mask*(phix*U(:,nim)+i*(phiy*U(:,nim)));
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,0)),'Nnodes','xo','yo');
        xn=xo;
        xn(length(xn))=xo(length(xo))-1;
        yn=yo;
        yn(length(yn))=yo(length(yo))-1;
        tmap=reshape(real(dmap),sizeim);
        unmap=tmap(xn,yn)-mean(tmap(:));
        tmap=reshape(imag(dmap),sizeim);
        vnmap=tmap(xn,yn)-mean(tmap(:));
        fac=0.2*max(sizeim)/max(abs(unmap(:)+i*vnmap(:)));
        limd=[min(abs(unmap(:)+i*vnmap(:))),max(abs(unmap(:)+i*vnmap(:)))];
        DeformedMesh(xn,yn,unmap,vnmap,limd,fac,['solve-femu'],0,1);
        pause
    end
    jj=jj+1;
    end
            res=norm(P-Pold)/nU0;
    disp(sprintf('After update # %d',ii));
        disp(sprintf('|P-Pold|/|Po|=%f',res));

    Pold=P;
    ii=ii+1;
end
disp(sprintf('Enlapsed time for resolution = %6.2f s',cputime-ttic));

if l>0
    invA=inv(M);
    tmp=0;
       for iim=1:nim
    if onpixel
        tmp=tmp+w(iim)*(phix+i*phiy)'*Ju(:,:,iim)*invA;
    else
        tmp=tmp+w(iim)*Ju(:,:,iim)*invA;
    end
       end
    C=Mc\tmp;
    sens=real(tmp'*C);
if length(P)==1
    display(sprintf('Sensitivity to image noise %5.3e / gray level',sqrt(2*full(sens))));
end
else
    sens=0;
end
pscale=2^(iscale-1);
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','sizeim');
load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'Xi','Yi');
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
    if strcmp(param.basis,'fem')
        mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
        load(mesh_file,'ng');
        if ng>0
load(fullfile('TMP','sample0'),'im0','sizeim');
             [im00]=mexInterpSpline(Xi+roi(1)-1,Yi+roi(3)-1,im0);
             im0=im00;
        end
    end

dynamic=max(im0(:))-min(im0(:));
mean0=mean(im0(:));
std0=std(im0(:));
im0=im0-mean0;


fid=fopen(fullfile('TMP',sprintf('%d_error_%d.mat',nmod,iscale-1)),'w');
fwrite(fid,min(255,dynamic));
        for jjm=1:nim
Ux=phix*U(:,jjm);
Uy=phiy*U(:,jjm);
NestedResidual(jjm);

disc=mask*disc;
    if dynamic>255
        disc=255*disc/dynamic;
    end

    fwrite(fid,abs(disc));
        end
fwrite(fid,sens);
        fclose(fid);

    function NestedResidual(kkk)
            if nim==1
                fildef=param0.deformed_image;
                im1=double(readim(fildef));
            else
                fildef=param0.deformed_image{kkk};
                im1=double(readim(fildef));
            end
            if length(size(im1))==3
                im1=mean(im1,3);
            end
            if reverse
                im1=im1';
            end

 
        sizeim1=size(im1);
        disc=(mexInterpSpline(min(sizeim1(1)-2,max(3,(Xi+Ux+roi(1)-1))),min(sizeim1(2)-2,max(3,(Yi+Uy+roi(3)-1))),im1));


        mean1=mean(disc(:));
        std1=std(disc(:));
        disc=disc-mean1;
        st=std0/std1;
        disc=(im0(:)-st*disc(:));

    end
end