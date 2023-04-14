function LoadMask(nmod)

load(fullfile('TMP','params'));
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=0;
end
param0=param;
%irestart=param.restart_image;
irestart=1;
if irestart
    load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
    nscale=param.nscale;

    if ~isfield(param0,'mask_file')
        loadfile=0;
    else
        filemask=param0.mask_file;
        if iscell(filemask)
            filemask=filemask{nmod};
        end
        if strcmp(filemask,'')
            loadfile=0;
        else
            loadfile=1;
        end

    end


    if ~loadfile
        mask=1;
        maskg=1;
        masks=1;
        for iscale=1:nscale
            save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask','maskg','masks');
        end
        if strcmp(param.basis,'fem')||strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri')||strcmp(param.basis,'nurbs-beam')
            iscale=1;
            maskn=1;unmasked_nodes=[];cut_nodes=[];
            for iscale=1:nscale
                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'cut_nodes','unmasked_nodes','maskn','-append');
            end
        end
        disp(sprintf('Mask Loaded for model %d...',nmod));

    else
        tic;
        iscale=1;
        filref=param0.mask_file;
        roi=param0.roi;
        load(filemask,'mask');
        if psample
      [Yi,Xi]=meshgrid(roi(3):psample:roi(4),roi(1):psample:roi(2));
  mask0=interp2(mask,Yi,Xi,'*linear');
else
        mask0=(mask(roi(1):roi(2),roi(3):roi(4)));
        end
        sizeim=size(mask0);
        mask=diag(sparse(mask0(:)));
        save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
        if strcmp(param.basis,'fem')||strcmp(param.basis,'nurbs')
            mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
            load(mesh_file,'ng');
            if ng==0
                if ~isfield(param,'mesh_file')
                load(mesh_file,'Nnodes');
                phi=CreateFiniteElementBasis(mesh_file,sizeim,1);
                maskn=sum(mask*phi>0,1)./sum(phi>0,1);
                tol=0.025+0./sum(phi>0,1);
                maskn=reshape(full(maskn(:)),Nnodes);
                tol=reshape(full(tol(:)),Nnodes);
                cut_nodes=find((maskn(:)>tol(:))&(maskn(:)<1-tol(:)));
                maskn=double(maskn>tol);
                unmasked_nodes=find(maskn(:));
                maskg=mask;
                masks=1;
                else
                maskn=1;unmasked_nodes=[];cut_nodes=[];
                maskg=1;
                masks=1;                    
                end
                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'cut_nodes','unmasked_nodes','maskn','masks','-append');
            else
                maskn=1;unmasked_nodes=[];cut_nodes=[];
                mask=1;
                maskg=1;
                masks=1;
                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'cut_nodes','unmasked_nodes','maskn','mask','masks','maskg','-append');
            end
        end
        for iscale=2:nscale
            mask0=CoarseImage(mask0);
            mask0=mask0==1;
            sizeim=size(mask0);
            mask=diag(sparse(mask0(:)));
            save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask');
            if strcmp(param.basis,'fem')||strcmp(param.basis,'nurbs')
                mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
                load(mesh_file,'Nnodes');
                phi=CreateFiniteElementBasis(mesh_file,sizeim,1);
                maskn=sum(mask*phi>0,1)./sum(phi>0,1);
                maskn=reshape(full(maskn(:)),Nnodes);
                cut_nodes=find((maskn(:)>0.1)&(maskn(:)<1-0.1));
                maskn=double(maskn>0.1);
                unmasked_nodes=find(maskn(:));

                save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'cut_nodes','unmasked_nodes','maskn','-append');
            end
        end


        disp(sprintf('Loading mask for model %d...%6.2f s',nmod,toc));
    end
end


end
