function AssembleEquilibriumGapOperator(iscale,nmod,domatrix)
if nargin<3, domatrix=true;end
tic();
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
dflag=isfield(param,'extrusion_parameters')||(length(param0.roi)==6);
if isfield(param,'crack_id')
    ic=param.crack_id;
else
    ic=1;
end
if domatrix load(fullfile('TMP',sprintf('%d_k_operator_%d',nmod,iscale-1)),'K');end
if dflag 
   load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,(iscale-1))),'Nddl_tot'); 
else
   load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,(iscale-1))),'Nddl_tot'); 
end
select=ones(Nddl_tot,1);
if strcmp(param.basis,'fem')

    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'unmasked_nodes');
    if dflag
        load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,iscale-1)),'Nnodes','selected');
    else
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nnodes','selected');
    end
    Nnods=prod(Nnodes);
    outer_nodes=~selected;
    if isfield(param0,'coupling_width')
        load(fullfile('TMP',sprintf('alrequin_%d',2)),'coupling_nodes');
        display('WARNING nmod assigned by hand !!!!');
        outer_nodes(coupling_nodes)=1;
    end
    if ~isempty(unmasked_nodes)
        outer_nodes=outer_nodes(unmasked_nodes);
        Nnods=length(unmasked_nodes);
    end
    reject=find(outer_nodes(:));
    select(reject)=0;
    select(Nnods+reject)=0;
    if dflag
        select(2*Nnods+reject)=0;
    end
    if isfield(param,'enrichment')
        load(fullfile('TMP',sprintf('%d_enrichment_%d',nmod,ic)),'enriched_nodes');
        enriched=zeros(prod(Nnodes),1);
        enriched(enriched_nodes)=1;
        if ~isempty(unmasked_nodes)
            enriched=enriched(unmasked_nodes);
        end
        enriched_unmasked_nodes=find(enriched);
        select(enriched_unmasked_nodes)=1;
        select(Nnods+enriched_unmasked_nodes)=1;
        if dflag
            select(2*Nnods+enriched_unmasked_nodes)=1;
        end
    end


    select=diag(sparse(select));
if domatrix, K=select*K;end
end
save(fullfile('TMP',sprintf('%d_kk_operator_%d',nmod,iscale-1)),'select');
if domatrix
    K=K'*K;
save(fullfile('TMP',sprintf('%d_kk_operator_%d',nmod,iscale-1)),'K','-append');
end


disp(sprintf('Computing equilibrium gap operator...%6.2f s',toc()));



end
