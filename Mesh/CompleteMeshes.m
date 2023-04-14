function CompleteMeshes(Zo,iscale,nmod,Xo,Yo)
if nargin<4
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
roi=param0.roi;
nurbs=strcmp(param.basis,'nurbs');
load(fullfile('TMP','sample0'),'sizeim');
pscale=2^(iscale-1);
mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1));
C=param0.calibration_data{1};
if nurbs&&iscale==1
                    load(mesh_file,'Px','Py');
                xo=Px(:);yo=Py(:);

else
load(mesh_file,'xo','yo');
xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
end
xo=(xo-1)+roi(1)-0.5*(sizeim(1)+1);
yo=(yo-1)+roi(3)-0.5*(sizeim(2)+1);
%[Xo,Yo]=GetXYFromabZ(C,xo,yo,Zo);
[Xo,Yo,Zo]=GetXYZFromxyZ1(C,xo,yo,Zo);
end
if nurbs&&iscale==1
     PX=Xo;PY=Yo;PZ=Zo;
     Xo=CPToNodes(PX,nmod,iscale);
     Yo=CPToNodes(PY,nmod,iscale);
     Zo=CPToNodes(PZ,nmod,iscale);
     save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'PX','PY','PZ','-append');   
end


save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Xo','Yo','Zo','-append');

end