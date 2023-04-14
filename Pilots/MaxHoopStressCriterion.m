function [K1eq,Oc]=MaxHoopStressCriterion(K1,K2)


rat=K1./K2;
  Oc = 2.*atan(0.25*(rat + sqrt(rat.*rat + 8.)));
K1eq=K1.*(cos(Oc/2)).^3-1.5*K2.*cos(Oc/2).*sin(Oc);
if any(K1eq<0)
    found=find(K1eq<0);
  Oc(found) = 2.*atan(0.25*(rat(found) - sqrt(rat(found).*rat(found) + 8.)));
  K1eq(found)=K1(found).*(cos(Oc(found)/2)).^3-1.5*K2(found).*cos(Oc(found)/2).*sin(Oc(found));
end




end