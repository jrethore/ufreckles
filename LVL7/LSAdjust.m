function [vphi]=LSAdjust(vphi,vpsi,psi)

vphi=vphi.*psi./(vpsi);
vphi=vphi.*double(psi>=0);


end