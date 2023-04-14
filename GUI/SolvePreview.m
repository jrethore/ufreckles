function [U]=SolvePreview(Uini,iscale,nmod,ijm)

persistent im0 M mean0 std0 roip R L MM Nnodes Mo
ttic=cputime;
load(fullfile('TMP','params'),'param');
param0=param;
restart=0;
onflight=1;
disc=0;
preview=1;
pgd=1;
if onflight
    global phiy phix Xi Yi phidf
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if ~isfield(param0,'deformed_image')
    mreader=VideoReader(param0.reference_image);
    nbf=mreader.NumberOfFrames;
    if isfield(param0,'number_of_frames')
        nbf=param0.number_of_frames;
    end
    dim=1;
    if isfield(param0,'video_sampling')
        dim=param0.video_sampling;
    end
    frames=2:dim:nbf;
    nim=length(frames);
else
    if iscell(param0.deformed_image)
        nim=size(param0.deformed_image,2);
    else
        nim=1;
    end
end
maxiter=50;
conv=1e-3;

pscale=2^(iscale-1);

im1=0;
U=Uini;
icamr=1;
if isempty(roip)
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nnodes','xo','yo');
    load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0','roi');
    load(fullfile('TMP',sprintf('%d_operator_%d',nmod,iscale-1)),'M','R','P');
    mean0=mean(im0(:));
    std0=std(im0(:));
    roip=((roi-1)/pscale)+1;
    if strcmp(param.basis,'fem')
    L=[1+0*xo,0*xo,-yo,xo,0*xo,yo;...
        0*yo,1+0*yo,xo,0*yo,yo,xo];
    else
        L=1;
    end
    MM=L'*(M)*L;
    Mo=M;
    M=M+R+P;
    
end
iz=1;

if ~isfield(param0,'deformed_image')
    im1=double(read(mreader,frames(ijm)));
else
    if ~iscell(param0.deformed_image)
        fildef=param0.deformed_image;
        im1=double(readim(fildef));
    else
        fildef=param0.deformed_image{1,ijm};
        im1=double(readim(fildef));
    end
end
if length(size(im1))==3
    im1=mean(im1,3);
end
if iscale>1
    for ip=1:(iscale-1)
        scale=2;
        imsiz0=size(im1);
        imsiz1=floor(imsiz0/2);
        nn=2*imsiz1;
        im1=im1(1:nn(1),1:nn(2));
        
        im1=reshape(im1,scale,prod(nn)/scale);
        im1=mean(im1,1);
        nn(1)=nn(1)/scale;
        im1=reshape(im1,nn);
        
        im1=im1';
        im1=reshape(im1,scale,prod(nn)/scale);
        im1=mean(im1,1);
        nn(2)=nn(2)/scale;
        im1=reshape(im1,nn([2,1]));
        im1=im1';
        
    end
end






if restart||ijm<4
    [Un,Vn]=rbt(im0,im1(roip(1):roip(2),roip(3):roip(4)));
    switch param.basis
        case 'fem'
    Ui=[(Un*pscale)*ones(prod(Nnodes),1);(Vn*pscale)*ones(prod(Nnodes),1)];
        case 'uni'
            Ui=zeros(6,1);
            Ui(1)=(Un*pscale);
            Ui(2)=(Vn*pscale);
    end
    dL=Inf;
    ii=1;
    while norm(dL)>0.001&&ii<100
        Ux=phix*Ui;
        Uy=phiy*Ui;
        NestedResidual(ijm);
        FF=phidf'*(disc);
        FF=L'*FF;
        dL=MM\FF;
        Ui=Ui+L*dL;
        ii=ii+1;
    end
else
    DU=U(:,ijm-1,iz)-U(:,ijm-2,iz);
    Ui=U(:,ijm-1,iz)+DU;
    
    if pgd
        LL=U(:,max(2,ijm-4):(ijm-1),iz);
        dL=Inf;
        ii=1;
        LMML=LL'*(Mo)*LL;
        while norm(dL)>0.001&&ii<10
            Ux=phix*Ui;
            Uy=phiy*Ui;
            NestedResidual(ijm);
            FF=phidf'*(disc);
            FF=LL'*FF;
            dL=LMML\FF
            Ui=Ui+LL*dL;
            ii=ii+1;
        end
    end
end
res=Inf;
ii=1;

while ( res>conv && ii< maxiter)
    Ux=phix*Ui;
    Uy=phiy*Ui;
    NestedResidual(ijm);
    
    F=phidf'*(disc);
    
    F=F-R*Ui;
    
    if isfield(param0,'solver')
        solv=param0.solver;
        
        if strcmp(solv,'lu')
            if ii==1
                [LT,UT]=lu(sparse(M));
            end
            dU=UT\(LT\F);
        elseif strcmp(solv,'chol')
            if ii==1
                Rc=chol(sparse(M));
            end
            dU=Rc\(Rc'\F);
        end
    else
        dU=M\F;
    end
%                disp(sprintf('At iteration # %d',ii));
    if ii==1
        nU0=max(sqrt(length(Ui)),norm(Ui+dU(:)));
    else
        if length(dU)>10
            res=norm(dU(:))/nU0;
        else
            res=max(abs(dU)./Ui);
        end
%                       disp(sprintf('|dU|=%f',res));
        
    end
    Ui=Ui+dU;
    ii=ii+1;
end
U(:,ijm,iz)=Ui;
%disp(sprintf('Enlapsed time for resolution = %6.2f s',cputime -ttic));
    function NestedResidual(kkk)
        
        disc=(mexInterpLinear((Xi-1)+Ux/pscale+roip(1),(Yi-1)+Uy/pscale+roip(3),im1));
        
        maske=disc<0;
        mean1=mean(disc(:));
        std1=std(disc(:));
        disc=disc-mean1;
        st=std0/std1;
        disc=(im0(:)-mean0-st*disc(:));
        disc(maske)=0;
        
        
    end
end
