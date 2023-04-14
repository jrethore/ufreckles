function [dummy]=UpdateOperators(nmod)
dummy=1;

load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
AssembleMechanicalOperator(1,nmod)
if isfield(param,'coupling_parameter')
load(fullfile('TMP',sprintf('%d_k_operator_%d',nmod,1-1)),'K');
load(fullfile('TMP',sprintf('%d_imodel',nmod)),'select');
K=select*K;
K=K'*K;
save(fullfile('TMP',sprintf('%d_kk_operator_%d',nmod,1-1)),'K','select');
end

end