function [P]=FEMUSolve(nmod,Pini)

load(fullfile('TMP','params'),'param');
param0=param;
maxiter=param0.iter_max;
conv=param0.convergance_limit;
if isfield(param0,'search_convergance_limit')
    conv=param0.search_convergance_limit;
end
if isfield(param0,'update_max')
    maxiter=param0.update_max;
end
if isfield(param0,'under_relaxation')
    alpha=param0.under_relaxation;
else
    alpha=1;
end
load('TMP/0_mesh_0')
P=Pini;
res=Inf;
ii=1;
Phist=zeros(length(P),maxiter);
Ehist=zeros(1,maxiter-1);
Phist(:,1)=P;

while ( res>conv && ii< maxiter)
        err=costfun(nmod,P);
%         figure
%         subplot(1,2,1)
%         trimesh(conn,xo,yo,err((1:prod(Nnodes))))
%         colorbar
%         subplot(1,2,2)
%         trimesh(conn,xo,yo,err(prod(Nnodes)+(1:prod(Nnodes))))
%         colorbar
        if ii==1
            erro=err;
        end
    J=repmat(0,[size(err,1),length(P)]);
   for ip=1:length(P)
        dPi=0*P;
        dPi(ip)=0.01*Pini(ip);
        erri=costfun(nmod,P+dPi);
        J(:,ip)=(err-erri)/dPi(ip);
   end
   M=J'*J;
   F=J'*err;
   dP=M\F;
       res=max(abs(dP./(P)));
    disp(sprintf('At iteration # %d, |dP|/|P|=%f',ii,res));
    disp(sprintf('Discrepancy =%6.2f Er/Ero  %6.2f\n',norm(err),norm(err)/norm(erro)));
    P=P+alpha*dP;
    disp(sprintf('Actual parameters values %6.2f\n',P));
    Ehist(ii)=norm(err);
    ii=ii+1;

    Phist(:,ii)=P;
end
niter=ii;
save(fullfile('TMP',sprintf('%d_femusolve_%d',nmod,0)),'Ehist','Phist','niter');
sens=full(inv(M));
save(fullfile('TMP',sprintf('%d_sens_%d.mat',nmod,0)),'sens');

end