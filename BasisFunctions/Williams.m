function phi=Williams(z,modes,ind,kappa)
            dist=abs(z);
            angl=angle(z);
            phi=zeros(length(z),length(modes)*length(ind));

icon=0;
for m=1:length(modes)
    for ii=1:length(ind)
        icon=icon+1;
        meq=modes(m);
        heq=ind(ii);
        [harm1,amp1]=KMPotFourrier(meq,heq,kappa);
        feq=0*dist(:);
        for kk=1:3
            feq=feq+amp1(kk)*exp(harm1(kk)*angl(:)*1i/2);
        end
        feq=feq.*(dist(:).^(heq/2));
        phi(:,icon)=feq;
    end
end


end