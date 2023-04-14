function [dummy]=UpdateModel(nmod,P)
dummy=1;

[dummy]=UpdateMaterial(nmod,P);
[dummy]=UpdateOperators(nmod);


end