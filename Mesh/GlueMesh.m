function [xo,yo,zo]=GlueMesh(operations,xo,yo,zo)
if nargin<4,zo=0*xo;end
for iop=1:length(operations)
    operation=operations{iop};
    optype=operation{1};
    switch optype
        case 'translate'
            vec=operation{2};
            xo=xo+vec(1);
            yo=yo+vec(2);
            zo=zo+vec(3);
        case 'rotate'
            z=operation{2};
            cent=operation{3};
            angl=operation{4};
            xo=xo-cent(1);
            yo=yo-cent(2);
            zo=zo-cent(3);
            
           if (sum(z==[0;0;1])==3)
            xi=xo*cos(angl)+yo*sin(angl);
            yi=yo*cos(angl)-xo*sin(angl);
            zi=zo;
           elseif (sum(z==[0;1;0])==3)
            xi=xo*cos(angl)-zo*sin(angl);
            yi=yo;
            zi=zo*cos(angl)+xo*sin(angl);
           elseif (sum(z==[1;0;0])==3)
            xi=xo;
            yi=yo*cos(angl)+zo*sin(angl);
            zi=zo*cos(angl)-yo*sin(angl);
           else
               error(sprintf('WRONG ROTATION AXIS'))
           end
            xo=xi+cent(1);
            yo=yi+cent(2);
            zo=zi+cent(3);
        case 'scale'
            cent=operation{2};
            fac=operation{3};
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