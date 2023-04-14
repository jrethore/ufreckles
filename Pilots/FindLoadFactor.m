function [lf,K1eq,Oc]=FindLoadFactor(K1s,K2s,K1t,K2t,K1c)


[K1eq,Oc]=MaxHoopStressCriterion(K1t,K2t);
 lf=K1c/K1eq;
 dlf=1;
 iter=0;
%  figure
%  hold on
while abs(dlf)>1.e-6&&iter<1000
[K1eq,Oc]=MaxHoopStressCriterion(K1s+lf*K1t,K2s+lf*K2t);
%plot(lf,K1eq,'o');
dlf=0.05*(1-K1eq/K1c)/lf;
lf=lf+dlf*lf;
iter=iter+1;
end
 [K1eq,Oc]=MaxHoopStressCriterion(K1s+lf*K1t,K2s+lf*K2t);
if iter==1000
    display('WARNING: FAILURE TO CONVERGE IN 1000 ITERATIONS IN FindLoadFactor.m');
end


% lfs=0:0.0001:0.2;
% 
% Kevo=0*lfs;
% tevo=0*lf;
% for ii=1:length(lfs)
%    lfi=lfs(ii); 
%    [K1eqi,Oci]=MaxHoopStressCriterion(K1s+lfi*K1t,K2s+lfi*K2t);
%    Kevo(ii)=K1eqi;
%    tevo(ii)=Oci;
% end
%   figure
% plot(lfs,Kevo)
% grid on
%   figure
% plot(lfs,tevo*180/pi)
% grid on


end