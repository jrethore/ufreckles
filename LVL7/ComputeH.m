function H=ComputeH(commut,lset)
siz=size(lset);

[Nplus, Nminus]=ComputeNabla(lset);

H=max(commut,repmat(0,siz)).*Nplus+min(commut,repmat(0,siz)).*Nminus;


