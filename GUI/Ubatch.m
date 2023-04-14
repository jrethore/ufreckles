close all
clear all
clc
fid=fopen('ubatch.ufr','r');
filnames=textscan(fid,'%s');
filnames=filnames{1};
for ii=1:length(filnames)/2
    cmd=filnames{2*ii-1};
    filname=filnames{2*ii};
    switch cmd
        case 'path'
            cd(filname);
        case 'jobname'
            [pp,jobname,ext]=fileparts(filname);
            if strcmp(ext,'.res')
                load(filname,'-mat','model')
                for iz=1:size(model.zone,2)
                    if model.zone{4,iz}==5
                        zone=model.zone(:,iz);
                        run_williams_curved(1,zone,zone{5},filname)
                    end
                end
                
            else
                if strcmp(ext,'.ufr')
                    [param,model]=readINPFile(filname);
                else
                    load(filname,'-mat','param','model');
                end
                if ~isfield(param,'ulive'), param.ulive=0;end
                if ~isfield(param,'thermo'), param.thermo=0;end
                if ~isfield(param,'psample'), param.psample=1;end
                switch model.basis
                    case 'fem'
                        if ~isfield(model,'phantom_nodes'), model.phantom_nodes=0;end
                        switch param.analysis
                            case 'correlation'
                                switch numel(model.mesh_size)
                                    case 2
                                        if isfield(param,'detect')
                                            if param.detect
                                                run_crack_detect_job(param,model);
                                            else
                                                run_fem_job(param,model);
                                            end
                                        else
                                            run_fem_job(param,model);
                                        end
                                    case 3
                                        run_fem_job_3D(param,model)
                                        
                                end
                            case 'mechanics'
                                if param.da>0
                                    run_crack_propa_job(param,model);
                                else
                                    run_fea_job(param,model);
                                end
                        end
                    case 'uni'
                        run_uni_job(param,model);
                    case {'beam','beam-nurbs'}
                        run_beam_job(param,model);
                        
                end
            end
    end
end