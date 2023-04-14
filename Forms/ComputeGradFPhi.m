function ComputeGradFPhi(iscale,nmod)
tic();
load(fullfile('TMP','params'),'param');
onflight=0;
if isfield(param,'onflight')
    onflight=param.onflight;
end
if onflight
    global phiy phix Xi Yi phidf
end
param0=param;
if iscell(param0.reference_image)
    ncams=length(param0.reference_image);
else
    ncams=1;
end
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
    if isfield(param,'gp_file')
        ng=1;
    end
else
    opti_grad=1;
    ng=0;
end
for ncam=1:ncams
    if ncam==1
        roi=param0.roi;
    else
        load(fullfile('TMP',sprintf('%d_mesh_%d',nmod*10^(ncam-1),(iscale-1))),'roi');
    end
    
    if ~onflight
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(ncam-1),(iscale-1))),'phix');
        load(fullfile('TMP',sprintf('%d_phiy_%d',nmod*10^(ncam-1),(iscale-1))),'phiy');
    end
    switch opti_grad
        case 1
            %            if ng==0
            load(fullfile('TMP',sprintf('sample%d_%d',ncam-1,iscale-1)),'im0');
            NestedFDgradient();
            if (ng>0)||((strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri'))&&(iscale==1))
                if ~onflight
                    load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(ncam-1),(iscale-1))),'Xi','Yi');
                end
                gradxi=mexInterpLinear(Xi,Yi,gradx);
                gradyi=mexInterpLinear(Xi,Yi,grady);
                gradx=gradxi;
                grady=gradyi;
            end
            %             else
            %                 load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(ncam-1),(iscale-1))),'Xi','Yi');
            %                 load(fullfile('TMP','sample0'),'im0');
            %                 [gradx,grady]=mexGradLinear((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
            %                 gradx=gradx*psample;grady=grady*psample;
            %             end
        otherwise
            % load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,(iscale-1))),'xo','yo');
            % XY=[xo;yo];
            % Xi=phix*XY;
            % Yi=phiy*XY;
            if ~onflight
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod*10^(ncam-1),(iscale-1))),'Xi','Yi');
            end
            roi=param0.roi;
            load(fullfile('TMP',sprintf('sample%d',ncam-1)),'im0');
            [gradx,grady]=mexGradSpline((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
            %       [gradx,grady]=mexGradLinear((Xi-1)*psample+roi(1),(Yi-1)*psample+roi(3),im0);
            gradx=gradx*psample;grady=grady*psample;
    end
    phidf=diag(sparse(gradx(:)))*phix+diag(sparse(grady(:)))*phiy;
    
    if ~onflight
        save(fullfile('TMP',sprintf('%d_phidf_%d',nmod*10^(ncam-1),iscale-1)),'phidf','-v7.3');
    end
end
disp(sprintf('Computing phidf for model %d...%6.2f s',nmod,toc()));
    function NestedFDgradient()
        
        gradx=mexFDGradient(im0);
        im0=im0';
        grady=mexFDGradient(im0);
        
        grady=grady';
        
        
        
        
    end

end


