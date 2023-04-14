function [phi]=LSPropagate(F,phi,trupix)

ndt=ceil(max(abs(F(:)))/0.25);
dtp=1/ndt;

tt=0;
    while tt<=ndt
phi=mexLSPropagate(phi,F,trupix,dtp);
tt=tt+1;
    end

end