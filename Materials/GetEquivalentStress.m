function Seq=GetEquivalentStress(model,S)
  crit='Von-Mises';
if isfield(model,'equivalent_stress')
crit=model.equivalent_stress;
end

switch crit
case 'Von-Mises'
Seq=sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2-S(:,1).*S(:,3)-S(:,2).*S(:,3)-S(:,1).*S(:,2)...
	 +3*(S(:,4).^2+S(:,5).^2+S(:,6).^2));
otherwise

error(sprintf('PLASTICITY MODEL %s IS NOT AVAILABLE',crit));
end 

end
