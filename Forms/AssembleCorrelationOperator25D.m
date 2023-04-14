function [gradx,grady]=ComputeGradF(iscale,nmod)

tic();
load(fullfile('TMP','params'),'param');
param0=param;
if isfield(param0,'sampling_factor')
    psample=param0.sampling_factor;
else
    psample=1;
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if iscale==1
    if isfield(param0,'opti_grad')
        opti_grad=param0.opti_grad;
    else
        opti_grad=1;
    end
    if isfield(param,'nb_gauss_points')
        ng=param.nb_gauss_points;
    else
        ng=0;
    end
else
    opti_grad=1;
        ng=0;
end
switch opti_grad
    case 1
if ng==0
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
            NestedFDgradient();
else
load(fullfile('TMP',sprintf('%d_phix_%d',nmod,(iscale-1))),'Xi','Yi');
 roi=param0.roi;
load(fullfile('TMP','sample0'),'im0');
        [gradx,grady]=mexGradLinear((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
       gradx=gradx*psample;grady=grady*psample;
   
end
    otherwise
load(fullfile('TMP',sprintf('%d_phix_%d',nmod,(iscale-1))),'Xi','Yi');
 roi=param0.roi;
load(fullfile('TMP','sample0'),'im0');
       [gradx,grady]=mexGradSpline((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
       gradx=gradx*psample;grady=grady*psample;
end
gradx=diag(sparse(gradx(:)));
grady=diag(sparse(grady(:)));

    function NestedFDgradient()

            gradx=mexFDGradient(im0);
            im0=im0';
            grady=mexFDGradient(im0);

                grady=grady';




    end

end











