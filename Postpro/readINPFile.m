function [param,model]=readINPFile(filename)

param=[];
model=[];
fid=fopen(filename,'r');
buf=textscan(fid,'%s');
fclose(fid);
buf=buf{1};
com={};
ended=1;
for ii=1:length(buf)
        tmp=buf{ii};
    if length(tmp)>5
       if  strcmp(tmp(1:6),'model.')||strcmp(tmp(1:6),'param.')||~ended
           if strcmp(tmp((end-2):end),'...')
               if ended
             com{end+1}=tmp(1:(end-3));
               else
             com{end}=[com{end},tmp(1:(end-3))];
               end
             ended=0;
           else
               if ended
           com{end+1}=  tmp;
               else
                   com{end}=[com{end},tmp];
               end
           ended=1;
           end
       end
    end
end
for ii=1:length(com)
%    com{ii}
eval(com{ii});
end
if isfield(param,'convergence_limit')
    param.convergance_limit=param.convergence_limit;
    param=rmfield(param,'convergence_limit');
end
model.zone={};
end