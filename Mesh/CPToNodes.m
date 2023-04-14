function N=CPToNodes(P,nmod,iscale)
if nargin<3,iscale=1;end
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
assert(strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri'),'YOU MUST USE NURBS NURBS BASIS TO DO THAT !');

dflag=isfield(param,'extrusion_parameters')||(length(param0.roi)==6);


if dflag
load(fullfile('TMP',sprintf('%d_3d_phio_%d',nmod,iscale-1)),'phio');
else
load(fullfile('TMP',sprintf('%d_phio_%d',nmod,iscale-1)),'phio');
end
dim=size(P,1)/size(phio,2);
N=zeros([size(phio,1)*dim,size(P,2),size(P,3)]);
phi=phio;
for id=1:dim-1
phi=blkdiag(phi,phio);
end
for i2=1:size(P,2)
    for i3=1:size(P,3)
        
        N(:,i2,i3)=phi*P(:,i2,i3);
    end
end

end
