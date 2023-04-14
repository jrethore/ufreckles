function [dtdu,tu]=ComputeCohesiveTractionAndTangentStiffness(un2,tn2,dtndun,du,dt)
tu=interp1(un2,tn2,du,'cubic','extrap');

dtdu=interp1(un2,dtndun,du,'linear','extrap');

for ip=1:length(dtdu)
    if du(ip)<0
        dtdu(ip)=10000*max(abs(dtndun));
        tu(ip)=-10000*max(abs(tn2));
    elseif du(ip)>max(un2)
        dtdu(ip)=0000*max(abs(dtndun));
        tu(ip)=0;
%         elseif dt(ip)<0&&du(ip)>0
%             dtdu(ip)=10000*max(abs(dtndun));
    end
    
end

end