function [harm1,amp1]=KMPotFourrier(imode,id,kappa)

harm1=[id (4-id) (-id)];
amp1=zeros(1,3);
if imode==1
        amp1(1)=kappa;
        amp1(2)=(-id/2);
        amp1(3)=(id/2+(-1)^id);
elseif imode==2
        amp1(1)=-1i*kappa;
        amp1(2)=-1i*(id/2);
        amp1(3)=1i*(id/2-(-1)^id);
end
% odd=abs(id/2)-abs(floor(id/2));
% if imode==1
%     if odd
%         amp1(1)=(-1)^((1-id)/2)*kappa;
%         amp1(2)=(-1)^((1-id)/2)*(-id/2);
%         amp1(3)=(-1)^((1-id)/2)*(id/2-1);
%     else
%         amp1(1)=kappa;
%         amp1(2)=(-id/2);
%         amp1(3)=(id/2+1);
%     end
% elseif imode==2
%     if odd
%         amp1(1)=i*(-1)^((1-id)/2)*kappa;
%         amp1(2)=i*(-1)^((1-id)/2)*(id/2);
%         amp1(3)=-i*(-1)^((1-id)/2)*(id/2+1);
%     else
%         amp1(1)=i*kappa;
%         amp1(2)=i*(id/2);
%         amp1(3)=-i*(id/2-1);
%     end
% end
end
