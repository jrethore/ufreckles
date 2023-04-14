function U=CorrectForDistortion(sizeim,xo,yo,Uo,Up)
%%
%Nnodes=[36,1260/36];
U=Uo;
conv=1.e-6;
maxiter=5000;
for iim=1:size(U,2)
X=(xo)/(sizeim(1))-0.5;
Y=(yo)/(sizeim(2))-0.5;
[phix,phiy]=GetCordinCorrectionFied(X,Y);
    Ux=(Uo(1:size(Uo,1)/2,iim)-phix*Up(:,iim));
    Uy=(Uo(size(Uo,1)/2+(1:size(Uo,1)/2),iim)-phiy*Up(:,iim));
    res=1;
    resp=2;
    iter=1;
%    figure

    while resp>res&&res>conv&&iter<maxiter
%    while iter<maxiter
        resp=res;
X=(xo+Ux)/(sizeim(1))-0.5;
Y=(yo+Uy)/(sizeim(2))-0.5;
[phix,phiy,phixx,phixy,phiyy,phiyx]=GetCordinCorrectionFied(X,Y);
%            detJ=detJ+double(abs(detJ)<1e-6);
rx=Uo(1:size(Uo,1)/2,iim)-Ux-phix*Up(:,iim);
ry=Uo(size(Uo,1)/2+(1:size(Uo,1)/2),iim)-Uy-phiy*Up(:,iim);
detJ=((1+phixx*Up(:,iim)/sizeim(1)).*(1+phiyy*Up(:,iim)/sizeim(2))-(phixy*Up(:,iim)).*(phiyx*Up(:,iim))/prod(sizeim));
            invJ=[((1+phiyy*Up(:,iim))/sizeim(2))./detJ,-(phixy*Up(:,iim)/sizeim(2))./detJ,-(phiyx*Up(:,iim)/sizeim(1))./detJ,(1+phixx*Up(:,iim)/sizeim(1))./detJ];
            dx=(invJ(:,1).*rx+invJ(:,2).*ry);
            dy=(invJ(:,3).*rx+invJ(:,4).*ry);
res=norm([dx;dy])/norm([Ux;Uy]);
% subplot(2,3,1)
%        imagesc(reshape(Ux,Nnodes))
%        colorbar
% subplot(2,3,2)
%        imagesc(reshape(rx,Nnodes))
%        colorbar
% subplot(2,3,3)
%        imagesc(reshape(dx,Nnodes))
%        colorbar
%        subplot(2,3,4)
%        imagesc(reshape(phix*Up(:,iim),Nnodes))
%        colorbar
%        subplot(2,3,5)
%        imagesc(reshape(Uo(1:size(Uo,1)/2,iim),Nnodes))
%        colorbar
% 
%        subplot(2,3,6)
%        imagesc(reshape(1./detJ,Nnodes))
%        colorbar
% 

%pause
            Ux=Ux+dx;
             Uy=Uy+dy;

iter=iter+1;
    end
    U(:,iim)=[Ux;Uy];
end
%%
end

