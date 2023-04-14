function [xg,yg,zg,wg]=GetGaussPointsWedge(nb_gauss_point,nb_sub_cell)

switch nb_gauss_point
    case {1,2}
        xg=[1/3;1/3];
        yg=[1/3;1/3];
        zg=[-0.577350269189626;0.577350269189626];
        wg=[0.5;0.5];

end
