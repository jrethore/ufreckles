function []=LoadParameters(param,nmod)
if ~exist('TMP')
    mkdir('TMP');
end
if ~exist('FIG')
    mkdir('FIG');
end
if ~exist('VTK')
    mkdir('VTK');
end
if nargin < 2
    if ~isfield(param,'stack_format')
        param.stack_format='bin';
    end
    
    if ~isfield(param,'pixel_size')
        param.pixel_size=1;
    end
    if ~isfield(param,'detect')
        param.detect=0;
    end
     if ~isfield(param,'thermo')
        param.thermo=0;
     end
save(fullfile('TMP','params'),'param');
disp(sprintf('General parameters loaded...'));
else
if ~isfield(param,'nscale')
    param.nscale=1;
end
save(fullfile('TMP',sprintf('%d_params',nmod)),'param');
disp(sprintf('Parameters for model %d loaded...',nmod));
end
 
end
