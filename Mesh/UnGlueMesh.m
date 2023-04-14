function [xo,yo,zo]=UnGlueMesh(nmod,xo,yo,zo)
if nargin<4,zo=0*xo;end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
operations=param.gluing_parameters;
for iop=length(operations):-1:1
    operation=operations{iop};
    optype=operation{1};
    switch optype
        case 'translate'
            vec=-operation{2};
            xo=xo+vec(1);
            yo=yo+vec(2);
            zo=zo+vec(3);
        case 'rotate'
            z=operation{2};
            assert(sum(z==[0;0;1])==3);
            cent=operation{3};
            angl=-operation{4};
            xo=xo-cent(1);
            yo=yo-cent(2);
            zo=zo-cent(3);
            xi=xo*cos(angl)+yo*sin(angl);
            yi=yo*cos(angl)-xo*sin(angl);
            zi=zo;
            xo=xi+cent(1);
            yo=yi+cent(2);
            zo=zi+cent(3);
        case 'scale'
            cent=operation{2};
            fac=1/operation{3};
            xo=xo-cent(1);
            yo=yo-cent(2);
            zo=zo-cent(3);
            xo=xo*fac;
            yo=yo*fac;
            zo=zo*fac;
            xo=xo+cent(1);
            yo=yo+cent(2);
            zo=zo+cent(3);

    end
end

end