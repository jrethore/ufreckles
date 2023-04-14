function [w]=GetWeight(weighting,n)


switch weighting
    case 'uniform'
        w=ones(n,1);
        w=w/sum(w);
    case 'linear'
        w=(1:n)';
        w=w/sum(w);

    otherwise
        error ('INVALID WEIGHTING');

end


end