function [lim]=FindMinMaxDisp(dmin,dmax)
  dmin=sign(dmin)*0.5*((dmin<0)+floor(2*abs(dmin)));
dmax=sign(dmax)*0.5*(-1*(dmax<0)+ceil(2*abs(dmax)));

lim=[dmin,dmax];

end
