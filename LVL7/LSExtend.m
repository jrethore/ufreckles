function [vphi]=LSExtend(vphi,phi,dist,trupix)

dt=0.5*trupix;
maxiter=floor(dist/dt);

for iter=1:maxiter
    
    [vphi]=mexLSReOrtho(vphi,phi,trupix,dt);
end
    
end



