function AssembleMechanicalOperator(iscale,nmod,step)
if nargin<3,step=[];end
    tic();
        K=AssembleStiffnessMatrix(iscale,nmod,step);
        save(fullfile('TMP',sprintf('%d_k_operator_%d',nmod,iscale-1)),'K','-v7.3');


%    disp(sprintf('Computing mechnical operator...%6.2f s',toc()));



end
