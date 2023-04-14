function [Uini]=InitializeVICSolution(U0,iscale,nmod)
%if nargout==1
load(fullfile('TMP','params'));
param0=param;
restart=1;
if isfield(param0,'restart')
    restart=param0.restart;
end
reverse=0;
if isfield(param0,'reverse_image')
    reverse=param0.reverse_image;
end
if iscell(param0.deformed_image)
    nim=length(param0.deformed_image);
else
    nim=1;
end
if isempty(U0)
    do_rbt=1;
    load(fullfile('TMP','sample0_0'),'im0');
else
    do_rbt=0;
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
tau=param0.transition_length;
nscale=param.nscale;
Uini=[];
im1=[];
if restart
    nmax=nim;
else
    if nim<2
        nmax=nim;
    else
        nmax=2;
    end
    if iscale==param.nscale
        nmax=1;
    else
        nmax=nim;
    end
end

for iim=1:nmax
    if do_rbt
        if nim==1
            fildef=param0.deformed_image;
        else
            fildef=param0.deformed_image{iim};
        end
        roi=param0.roi;
        im1=double(readim(fildef));
        if numel(size(im1))==3
            im1=mean(im1,3);
        end
        if reverse
            im1=im1';
        end
        [U,V]=rbt(im0,im1(roi(1):roi(2),roi(3):roi(4)));
        U0=[U;V];
    end
    sbasis=param.basis;

    switch sbasis
        case 'vic-nurbs'
            if do_rbt
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod,iscale-1)),'phix','wdetJ');
                load(fullfile('TMP',sprintf('%d_phiy_%d',nmod,iscale-1)),'phiy');
                load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
                Ux=repmat(U0(1),size(phix,1),1);
                Uy=repmat(U0(2),size(phix,1),1);
                M=phix'*mask*wdetJ*phix+phiy'*mask*wdetJ*phiy;
                F=phix'*mask*wdetJ*Ux+phiy'*mask*wdetJ*Uy;
            else
                load(fullfile('TMP',sprintf('%d_phi_%d',nmod,iscale)),'phi');
                Ug=phi*U0(:,iim);
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod,iscale-1)),'wdetJ');
                load(fullfile('TMP',sprintf('%d_phi_%d',nmod,iscale-1)),'phi');
                load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
                M=phi'*mask*wdetJ*phi;
                F=phi'*mask*wdetJ*Ug;
            end
                Ui=M\F;
            if isempty(Uini)
                Uini=zeros(length(Ui),nim);
            end
%             if iscale==nscale
%                 Ui=Ui+tau*2^(iscale-1);
%             else
%                 Ui=Ui-tau;
%             end
                Uini(:,iim)=Ui;

        otherwise
            error ('INVALID FUNCTIONAL BASIS');

    end


end

end