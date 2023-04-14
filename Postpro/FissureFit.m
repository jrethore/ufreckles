function [dect,A]=FissureFit(U,V,W,crack,front,mask0,parsif)

Rcore=parsif.mask_radius;
crackedgemin=parsif.mask_width;
crackedgemax=crackedgemin;
if isfield(parsif,'mask_width_max')
    crackedgemax=parsif.mask_width_max;
end
kappa=parsif.kolosov;
nfit=parsif.km_indices;
decx=Inf;
maxiter=20;
iter=0;
dect=0;

while (abs(decx)>2)&&(iter<maxiter)

    if isinf(decx),decx=0;end
front=front+decx;
Zp=front+i*crack;
 mask=~((abs(Zp)<Rcore)|((abs(imag(Zp))<crackedgemin+(crackedgemax-crackedgemin)*abs(real(Zp))/abs(min(real(Zp(:)))))&(real(Zp)<0)));
mask=mask.*mask0;

            mask=diag(sparse(mask(:)));
[Uref,Vref,Wref]=CreateKMBasis(Zp,nfit,kappa);
            


iter=iter+1;
end

end
