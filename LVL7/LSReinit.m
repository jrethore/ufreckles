function [ls]=LSReinit(ls,dist,trupix)
dt=0.5*trupix;
maxiter=floor(dist/dt);

if numel(size(ls))==2
    for iter=1:maxiter
        ls=mexLSReinit0(ls,trupix,dt);
    end
elseif numel(size(ls))==3
    for iter=1:maxiter
        ls=mexLSReinit3D(ls,trupix,dt);
    end
else
    error('NOT CODED YET');
end
end


