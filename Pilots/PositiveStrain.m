function Pp=PositiveStrain(exx,eyy,exy)
    delta=sqrt((exx-eyy).^2+(exy).^2);
    m=(-exx+eyy+delta)./(exy);
    m(isnan(m))=0;
    delta(isinf(m))=-delta(isinf(m));
    m(isinf(m))=0;
    nn=1./(1+m.^2);
    
    
    l1=0.5*((exx+eyy)+(delta));
    Pp=diag(sparse(0.5*(l1+abs(l1))))*[nn,nn.*m.^2,2*nn.*m];
    
    l1=0.5*((exx+eyy)-(delta));
    Pp=Pp+diag(sparse(0.5*(l1+abs(l1))))*[nn.*m.^2,nn,-2*nn.*m];
end