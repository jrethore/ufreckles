function [I]=GetInvariant(E,dim,order)
switch dim
    case 2
        switch order
            case 1
                np=size(E,1)/2;
                I=sqrt(E(0*np+(1:np),:).^2+E(1*np+(1:np),:).^2);
            case 2
                np=size(E,1)/3;
                I=[E(0*np+(1:np),:)+E(1*np+(1:np),:);sqrt(E(0*np+(1:np),:).^2+E(1*np+(1:np),:).^2-...
E(1*np+(1:np),:).*E(0*np+(1:np),:)+4*E(2*np+(1:np),:).^2)];
            case 23
                np=size(E,1)/3;
                I=[E(0*np+(1:np),:)+E(1*np+(1:np),:);sqrt(E(0*np+(1:np),:).^2+E(1*np+(1:np),:).^2-...
E(1*np+(1:np),:).*E(0*np+(1:np),:)+3*E(2*np+(1:np),:).^2)];
                
            otherwise
        end
   otherwise
        error('NOT CODED YET')
end


end






















