function [psi]=LSOrtho(psi,phi,dist,trupix)

dt=0.5*trupix;
maxiter=floor(dist/dt);

for iter=1:maxiter
    
if numel(size(psi))==2
    [psi]=mexLSReOrtho(psi,phi,trupix,dt);
elseif numel(size(psi))==3
    [psi]=mexLSReOrtho3D(psi,phi,trupix,dt);
end
end
    
end



