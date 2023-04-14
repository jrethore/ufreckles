function [gradf,dgradf,beta]=YieldSurfaceGradient(model,S,Seq,H)

  crit='Von-Mises';
if isfield(model,'equivalent_stress')
crit=model.equivalent_stress;
end

switch crit
case 'Von-Mises'
Seq=0.5./Seq;
gradf=[(2*S(:,1)-(S(:,2)+S(:,3))).*Seq,...
       (2*S(:,2)-(S(:,1)+S(:,3))).*Seq,...
       (2*S(:,3)-(S(:,1)+S(:,2))).*Seq,...
       6*S(:,4).*Seq,...
       6*S(:,5).*Seq,...
       6*S(:,6).*Seq];

otherwise
error(sprintf('PLASTICITY MODEL %s IS NOT AVAILABLE',crit));

end 
E=model.Eelast;
inde=[1,1,1,3,3,3,2,2,2,2,2,2];
gradtodgrad=sparse([1:6,2,1,2,3,3,1],[1:6,1,2,3,1,2,3],E(inde),6,6);
dgradf=gradf*gradtodgrad;
if nargout>2
beta=sparse(1./(sum(gradf.*dgradf,2)+H));
end


end
