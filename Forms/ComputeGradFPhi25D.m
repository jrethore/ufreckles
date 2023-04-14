function ComputeGradFPhi25D(iscale,nmod)
tic();
load(fullfile('TMP','params'),'param');
param0=param;
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
if iscale==1
    if isfield(param0,'opti_grad')
        opti_grad=param0.opti_grad;
    else
        opti_grad=1;
    end
else
    opti_grad=1;
end
ncamr=length(param0.reference_image);
ncamd=size(param0.deformed_image,1);


switch param0.analysis
    case 'topography'
    filims=param0.deformed_image;
    imin=2;
ncam=ncamd;
pscale=2^(iscale-1);
phidft=sparse(0);
    filim=filims{1};
    im0=double(imread(filim));
    if length(size(im0))==3
        im0=mean(im0,3);
    end
    if reverse
        im0=im0';
    end

    if iscale>1
        NestedCoarseImage(pscale)
    end
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(1-1),(iscale-1))),'Xi','Yi');
    Xi=((Xi-1)*pscale)/pscale+1;
    Yi=((Yi-1)*pscale)/pscale+1;
    switch opti_grad
        case 1
                    gradxp=mexFDGradient(im0);
        im0=im0';
        gradyp=mexFDGradient(im0);

        gradyp=gradyp';

                                gradx=mexInterpLinear(Xi,Yi,gradxp);
                    grady=mexInterpLinear(Xi,Yi,gradyp);
%            [gradx,grady]=mexGradLinear(Xi,Yi,im0);
        otherwise
            [gradx,grady]=mexGradSpline(Xi,Yi,im0);
    end
for icam=imin:ncam
 
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),(iscale-1))),'phix','Xi','Yi');
    load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icam-1),(iscale-1))),'phiy');
    phidf=diag(sparse(gradx(:)))*phix+diag(sparse(grady(:)))*phiy;
    phidft=phidft+phidf;
    save(fullfile('TMP',sprintf('%d_phidf_%d',10^(icam)*nmod,iscale-1)),'phidf');
end
 phidf=phidft;
 save(fullfile('TMP',sprintf('%d_phidf_%d',nmod,iscale-1)),'phidf');
    case {'correlation','mechanics'}

    filims=param0.reference_image;
    imin=1;
ncam=ncamr;

pscale=2^(iscale-1);
for icam=imin:ncam
 
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(icam-1),(iscale-1))),'phix','Xi','Yi');
    load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(icam-1),(iscale-1))),'phiy');
    filim=filims{icam};
    im0=double(imread(filim));
    if length(size(im0))==3
        im0=mean(im0,3);
    end
    if reverse
        im0=im0';
    end

    if iscale>1
        NestedCoarseImage(pscale)
    end
    Xi=((Xi-1)*pscale)/pscale+1;
    Yi=((Yi-1)*pscale)/pscale+1;
    switch opti_grad
        case 1
                    gradxp=mexFDGradient(im0);
        im0=im0';
        gradyp=mexFDGradient(im0);

        gradyp=gradyp';

                                gradx=mexInterpLinear(Xi,Yi,gradxp);
                    grady=mexInterpLinear(Xi,Yi,gradyp);

            
            
 %           [gradx,grady]=mexGradLinear(Xi,Yi,im0);
        otherwise
            [gradx,grady]=mexGradSpline(Xi,Yi,im0);
    end
    phidf=diag(sparse(gradx(:)))*phix+diag(sparse(grady(:)))*phiy;
    save(fullfile('TMP',sprintf('%d_phidf_%d',10^(icam)*nmod,iscale-1)),'phidf');
end
end
disp(sprintf('Computing phidf for model %d...%6.2f s',nmod,toc()));
    function NestedCoarseImage(scale)

        imsiz0=size(im0);
        imsiz1=floor(imsiz0/2);
        nn=2*imsiz1;
        im0=im0(1:nn(1),1:nn(2));

        im0=reshape(im0,scale,prod(nn)/scale);
        im0=mean(im0,1);
        nn(1)=nn(1)/scale;
        im0=reshape(im0,nn);

        im0=im0';
        im0=reshape(im0,scale,prod(nn)/scale);
        im0=mean(im0,1);
        nn(2)=nn(2)/scale;
        im0=reshape(im0,nn([2,1]));
        im0=im0';

    end
end

