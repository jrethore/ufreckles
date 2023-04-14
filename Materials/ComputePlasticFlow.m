function [Sy,H]=ComputePlasticFlow(model,Ep)

switch model.hardening
    case 'swift'
        B=model.B;
        C=model.C;
        d=model.delta;
        if isfield(model,'H')
            Sy=sparse(model.H*Ep+B*(C+Ep).^d);
            H=sparse(model.H+d*B*(C+Ep).^(d-1));
        else
            Sy=sparse(B*(C+Ep).^d);
            H=sparse(d*B*(C+Ep).^(d-1));
        end
    case 'exp'
        Syo=model.Sy;
        Q=model.Q;
        b=model.b;
        Sy=sparse(Syo+Q*(1-exp(-b*Ep)));
        H=sparse(b*Q*exp(-b*Ep));
    case 'lin'
        Syo=model.Sy;
        Ho=model.H;
        Sy=sparse(Syo+Ho*Ep);
        H=sparse(Ho+0*Ep);
    case 'poly'
        Syi=model.Sy;
        Epi=model.Ep;
        Ep=full(Ep);
        Sy=interp1(Epi,Syi,Ep,'linear','extrap');
        Sy=sparse(Sy);
        de=1.e-6;
        ei=de:de:max(Epi);
        si=interp1(Epi,Syi,ei,'linear','extrap');
        hi=gradient(si)./gradient(ei);
        H=interp1([0,ei],[model.young,hi],Ep,'linear','extrap');
        H=sparse(H);
    otherwise
        error(sprintf('HARDENING LAW %s IS NOT AVAILABLE',crit));

end





end
