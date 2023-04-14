function [k]=GSIF(model,U,xi,yi,conng,pix2m)
if nargin<6,pix2m=1;end
k=cell(2*size(model.zone,2),5);
HH=model.material_parameters.H;
AA=model.material_parameters.A;
%
ii=[1,2,3];
jj=circshift(ii,[0,-1]);
kk=circshift(ii,[0,-2]);

for iz=1:size(model.zone,2)
    k{2*(iz-1)+1,3}=0;
    k{2*(iz-1)+2,3}=0;
    zone=model.zone(:,iz);
    htip=real(zone{7});
    switch zone{4}
        case 5 %CRACK
            dzone=zone{8};
            xyc=zone{2};
            if dzone>0
                face_elts=zone{10};
                nconn=zone{12};
                dfac=max(0.25,htip/dzone);
                xoo=[-dzone;dzone;dzone;-dzone]*1.5;
                yoo=[-dzone;-dzone;dzone;dzone]*1.5;
                xo=[-dzone;dzone;dzone;-dzone];
%                yo=[-dzone;-dzone;dzone;dzone];
                yo=[0;0;dzone;dzone];
                conn=1:4;
                elt=4;
                rint=0;
                ng=0;
%                ns=[1,1]*50;
                ns=5*[2,1]*round(dzone/htip);
                ns=25*[2,1];
                Nnodes=[length(xo),1,1];
                Nelems=[length(elt),1,1];
                zo=1;
                jmesh_file=fullfile('TMP','jmesh');
                save(jmesh_file,'rint','xo','yo','zo','Nnodes','Nelems','conn','ng','ns','elt')
                [wdetJ,inde]=GetWeigthDetJ(jmesh_file,[1,1],1,'sub_cells');
                wdetJ=diag(wdetJ);
                wdetJ=[wdetJ;flipud(wdetJ)];
                phig=CreateFiniteElementBasis(jmesh_file,[1,1],1,[],'sub_cells');
                xg=phig*xo;
                xg=[xg;flipud(xg)];
                yg=phig*yo;
                yg=[yg;-flipud(yg)];
                for itip=1:1+(zone{9}>0)
                    exx=0*xg+NaN;
                    eyy=0*xg+NaN;
                    ux=0*xg+NaN;
                    uy=0*xg+NaN;
                    uxy=0*xg+NaN;
                    uyx=0*xg+NaN;
                    uxxx=0*xg+NaN;
                    uxxy=0*xg+NaN;
                    uxyy=0*xg+NaN;
                    uyxx=0*xg+NaN;
                    uyxy=0*xg+NaN;
                    uyyy=0*xg+NaN;
                    switch itip
                        case 1
                            xytip=xyc(1,:);
                            t=-diff(xyc(1:2,:),[],1);
                        case 2
                            xytip=xyc(end,:);
                            t=diff(xyc(end-1:end,:),[],1);
                    end
                    t=t/norm(t);
                    n=[-t(2),t(1)];
                    rot=[t(1),t(2);n(1),n(2)];
                    irot=inv(rot);
                    xgglo=irot(1,1)*xg+irot(1,2)*yg+xytip(1);
                    ygglo=irot(2,1)*xg+irot(2,2)*yg+xytip(2);
                    xoglo=irot(1,1)*xoo+irot(1,2)*yoo+xytip(1);
                    yoglo=irot(2,1)*xoo+irot(2,2)*yoo+xytip(2);
                    xc=mean(xi(conng),2);
                    yc=mean(yi(conng),2);
                    inbox=find(inpolygon(xc,yc,xoglo([1:4,1]),yoglo([1:4,1])));
                    cin=conng(inbox,:);
                    [~,~,ty,tx]=GetSignedDistanceToCrack(xyc,xgglo+1i*ygglo);
                    ty=-ty;
                    nt=abs(tx+1i*ty);
                    tx=tx./nt;
                    ty=ty./nt;
                    for ie=1:size(cin,1)
                        inods=cin(ie,:);
                        xn=xi(inods);
                        yn=yi(inods);
                        
                        [xig,yig,wig]=GetGaussPointsTriangle(0,1,xn,yn,xgglo,ygglo);
                        if any(wig>0)
                            xig=xig(wig);
                            yig=yig(wig);
                            Nt=[1-xig-yig,xig,yig];
                            xig=Nt*xn;
                            yig=Nt*yn;
                            
                            ng=numel(xig);
                            N=zeros(ng,18);
                            Nx=zeros(ng,18);
                            Ny=zeros(ng,18);
                            Nxx=zeros(ng,18);
                            Nyy=zeros(ng,18);
                            Nxy=zeros(ng,18);
                            
                            A=0.5*(xn(2)*yn(3)-xn(3)*yn(2)+xn(3)*yn(1)-xn(1)*yn(3)+xn(1)*yn(2)-xn(2)*yn(1));
                            ao=xn(jj).*yn(kk)-xn(kk).*yn(jj);
                            bo=yn(jj)-yn(kk);
                            co=xn(kk)-xn(jj);
                            for in=1:3
                                a=circshift(ao,-(in-1));
                                b=circshift(bo,-(in-1));
                                c=circshift(co,-(in-1));
                                
                                L=0.5*((1+0*xig)*a'+xig*b'+yig*c')/A;
                                
                                r=diag(1./(b.^2+c.^2))*(-b*b'-c*c');
                                
                                a1=[b,c]/(2*A);
                                a2=[b.^2,c.^2,b.*c;...
                                    2*b.*circshift(b,-1),2*c.*circshift(c,-1),b.*circshift(c,-1)+c.*circshift(b,-1)]/(4*A.^2); % The second derivative
                                N0=zeros(numel(xig),1);
                                N1=zeros(numel(xig),3);
                                N2=zeros(numel(xig),6);
                                for ij=1:6
                                    
                                    switch ij
                                        case 1
                                            N0= L(:,1).^5+ 5*(L(:,1).^4).*(L(:,2))+ 5*(L(:,1).^4).*(L(:,3))+10*(L(:,1).^3).*(L(:,2).^2)+ ...
                                                10*(L(:,1).^3).*(L(:,3).^2)+ 20*(L(:,1).^3).*(L(:,2).*L(:,3))+30*r(2,1)*(L(:,1).^2).*L(:,2).*(L(:,3).^2)+...
                                                30*r(3,1)*(L(:,1).^2).*L(:,3).*(L(:,2).^2);
                                            N1(:,1)=5*L(:,1).^4+20*(L(:,1).^3).*(L(:,2))+20*(L(:,1).^3).*(L(:,3))+30*(L(:,1).^2).*(L(:,2).^2)+ ...
                                                30*(L(:,1).^2).*(L(:,3).^2)+ 60*(L(:,1).^2).*(L(:,2).*L(:,3))+60*r(2,1)*(L(:,1) ).*L(:,2).*(L(:,3).^2)+...
                                                60*r(3,1)*(L(:,1)   ).*L(:,3).*(L(:,2).^2);
                                            N1(:,2)=5*(L(:,1).^4)+20*(L(:,1).^3).*(L(:,2))+ 20*((L(:,1).^3).*L(:,3))+30*r(2,1)*(L(:,1).^2).*(L(:,3).^2)+...
                                                60*r(3,1)*(L(:,1).^2).*L(:,3).*(L(:,2));
                                            N1(:,3)=5*(L(:,1).^4)+20*(L(:,1).^3).*(L(:,3))+ 20*(L(:,1).^3).*(L(:,2))+60*r(2,1)*(L(:,1).^2).*L(:,2).*(L(:,3))+...
                                                30*r(3,1)*(L(:,1).^2).*(L(:,2).^2);
                                            
                                            N2(:,1)=20*L(:,1).^3+60*(L(:,1).^2).*(L(:,2))+60*(L(:,1).^2).*(L(:,3))+60*(L(:,1)).*(L(:,2).^2)+...
                                                60*(L(:,1)).*(L(:,3).^2)+120*(L(:,1)).*(L(:,2).*L(:,3))+60*r(2,1)*L(:,2).*(L(:,3).^2)+...
                                                60*r(3,1)*L(:,3).*(L(:,2).^2);
                                            N2(:,2)=20*(L(:,1).^3)+60*r(3,1)*(L(:,1).^2).*L(:,3);
                                            N2(:,3)=20*(L(:,1).^3)+60*r(2,1)*(L(:,1).^2).*L(:,2);
                                            N2(:,4)=20*(L(:,1).^3)+60*(L(:,1).^2).*(L(:,2))+ 60*((L(:,1).^2).*L(:,3))+60*r(2,1)*(L(:,1)).*(L(:,3).^2)+...
                                                120*r(3,1)*(L(:,1)).*L(:,3).*(L(:,2));
                                            N2(:,5)=20*((L(:,1).^3))+60*r(2,1)*(L(:,1).^2).*(L(:,3))+...
                                                60*r(3,1)*(L(:,1).^2).*(L(:,2));
                                            N2(:,6)=20*(L(:,1).^3)+60*(L(:,1).^2).*(L(:,3))+ 60*(L(:,1).^2).*(L(:,2))+120*r(2,1)*(L(:,1)).*L(:,2).*(L(:,3))+...
                                                60*r(3,1)*(L(:,1)).*(L(:,2).^2);
                                            
                                            
                                        case 2
                                            
                                            N0=c(3)*(L(:,1).^4).*L(:,2)-c(2)*(L(:,1).^4).*L(:,3)+4*c(3)*(L(:,1).^3).*(L(:,2).^2)-...
                                                4*c(2).*(L(:,1).^3).*(L(:,3).^2)+4*(c(3)-c(2))*(L(:,1).^3).*L(:,2).*L(:,3)-...
                                                (3*c(1)+15*r(2,1)*c(2))*(L(:,1).^2).*L(:,2).*(L(:,3).^2)+(3*c(1)+15*r(3,1)*c(3))*(L(:,1).^2).*L(:,3).*(L(:,2).^2);
                                            N1(:,1)=4*c(3)*(L(:,1).^3).*L(:,2)-4*c(2)*(L(:,1).^3).*L(:,3)+12*c(3)*(L(:,1).^2).*(L(:,2).^2)-...
                                                12*c(2).*(L(:,1).^2).*(L(:,3).^2)+12*(c(3)-c(2))*(L(:,1).^2).*L(:,2).*L(:,3)-...
                                                2*(3*c(1)+15*r(2,1)*c(2))*(L(:,1)).*L(:,2).*(L(:,3).^2)+2*(3*c(1)+15*r(3,1)*c(3))*(L(:,1)).*L(:,3).*(L(:,2).^2);
                                            N1(:,2)=c(3)*(L(:,1).^4)+8*c(3)*(L(:,1).^3).*(L(:,2))+4*(c(3)-c(2))*(L(:,1).^3).*L(:,3)-...
                                                (3*c(1)+15*r(2,1)*c(2))*(L(:,1).^2).*(L(:,3).^2)+2*(3*c(1)+15*r(3,1)*c(3))*(L(:,1).^2).*L(:,3).*(L(:,2));
                                            N1(:,3)=-c(2)*(L(:,1).^4)-8*c(2).*(L(:,1).^3).*(L(:,3))+4*(c(3)-c(2))*(L(:,1).^3).*L(:,2)-...
                                                2*(3*c(1)+15*r(2,1)*c(2))*(L(:,1).^2).*L(:,2).*(L(:,3))+(3*c(1)+15*r(3,1)*c(3))*(L(:,1).^2).*(L(:,2).^2);
                                            
                                            
                                            N2(:,1)=12*c(3)*(L(:,1).^2).*L(:,2)-12*c(2)*(L(:,1).^2).*L(:,3)+24*c(3)*(L(:,1)).*(L(:,2).^2)-...
                                                24*c(2).*(L(:,1)).*(L(:,3).^2)+24*(c(3)-c(2))*(L(:,1)).*L(:,2).*L(:,3)-...
                                                2*(3*c(1)+15*r(2,1)*c(2)).*L(:,2).*(L(:,3).^2)+2*(3*c(1)+15*r(3,1)*c(3))*L(:,3).*(L(:,2).^2);
                                            N2(:,2)=8*c(3)*(L(:,1).^3)+2*(3*c(1)+15*r(3,1)*c(3))*(L(:,1).^2).*L(:,3);
                                            N2(:,3)=-8*c(2).*(L(:,1).^3)-2*(3*c(1)+15*r(2,1)*c(2))*(L(:,1).^2).*L(:,2);
                                            N2(:,4)=4*c(3)*(L(:,1).^3)+24*c(3)*(L(:,1).^2).*(L(:,2))+12*(c(3)-c(2))*(L(:,1).^2).*L(:,3)-...
                                                2*(3*c(1)+15*r(2,1)*c(2))*(L(:,1)).*(L(:,3).^2)+4*(3*c(1)+15*r(3,1)*c(3))*(L(:,1)).*L(:,3).*(L(:,2));
                                            N2(:,5)=4*(c(3)-c(2))*(L(:,1).^3)-2*(3*c(1)+15*r(2,1)*c(2))*(L(:,1).^2).*(L(:,3))+...
                                                2*(3*c(1)+15*r(3,1)*c(3))*(L(:,1).^2).*(L(:,2));
                                            N2(:,6)=-4*c(2)*(L(:,1).^3)-24*c(2).*(L(:,1).^2).*(L(:,3))+12*(c(3)-c(2))*(L(:,1).^2).*L(:,2)-...
                                                4*(3*c(1)+15*r(2,1)*c(2))*(L(:,1)).*L(:,2).*(L(:,3))+2*(3*c(1)+15*r(3,1)*c(3))*(L(:,1)).*(L(:,2).^2);
                                            
                                        case 3
                                            
                                            N0=-b(3)*(L(:,1).^4).*L(:,2)+b(2)*(L(:,1).^4).*L(:,3)-4*b(3)*(L(:,1).^3).*(L(:,2).^2)+4*b(2)*(L(:,1).^3).*(L(:,3).^2)+...
                                                4*(b(2)-b(3))*(L(:,1).^3).*(L(:,2)).*(L(:,3))+(3*b(1)+15*r(2,1)*b(2))*(L(:,1).^2).*(L(:,2)).*(L(:,3).^2)-...
                                                (3*b(1)+15*r(3,1)*b(3))*(L(:,1).^2).*L(:,3).*(L(:,2).^2);
                                            N1(:,1)=-4*b(3)*(L(:,1).^3).*L(:,2)+4*b(2)*(L(:,1).^3).*L(:,3)-12*b(3)*(L(:,1).^2).*(L(:,2).^2)+12*b(2)*(L(:,1).^2).*(L(:,3).^2)+...
                                                12*(b(2)-b(3))*(L(:,1).^2).*(L(:,2)).*(L(:,3))+2*(3*b(1)+15*r(2,1)*b(2))*(L(:,1)).*(L(:,2)).*(L(:,3).^2)-...
                                                2*(3*b(1)+15*r(3,1)*b(3))*(L(:,1)).*L(:,3).*(L(:,2).^2);
                                            N1(:,2)=-b(3)*(L(:,1).^4)-8*b(3)*(L(:,1).^3).*(L(:,2))+...
                                                4*(b(2)-b(3))*(L(:,1).^3).*(L(:,3))+(3*b(1)+15*r(2,1)*b(2))*(L(:,1).^2).*(L(:,3).^2)-...
                                                2*(3*b(1)+15*r(3,1)*b(3))*(L(:,1).^2).*L(:,3).*(L(:,2));
                                            N1(:,3)=b(2)*(L(:,1).^4)+8*b(2)*(L(:,1).^3).*(L(:,3))+...
                                                4*(b(2)-b(3))*(L(:,1).^3).*(L(:,2))+2*(3*b(1)+15*r(2,1)*b(2))*(L(:,1).^2).*(L(:,2)).*(L(:,3))-...
                                                (3*b(1)+15*r(3,1)*b(3))*(L(:,1).^2).*(L(:,2).^2);
                                            
                                            N2(:,1)=-12*b(3)*(L(:,1).^2).*L(:,2)+12*b(2)*(L(:,1).^2).*L(:,3)-24*b(3)*(L(:,1)).*(L(:,2).^2)+24*b(2)*(L(:,1)).*(L(:,3).^2)+...
                                                24*(b(2)-b(3))*(L(:,1)).*(L(:,2)).*(L(:,3))+2*(3*b(1)+15*r(2,1)*b(2))*(L(:,2)).*(L(:,3).^2)-...
                                                2*(3*b(1)+15*r(3,1)*b(3))*L(:,3).*(L(:,2).^2);
                                            N2(:,2)=-8*b(3)*(L(:,1).^3)-2*(3*b(1)+15*r(3,1)*b(3))*(L(:,1).^2).*L(:,3);
                                            N2(:,3)=8*b(2)*(L(:,1).^3)+2*(3*b(1)+15*r(2,1)*b(2))*(L(:,1).^2).*(L(:,2));
                                            N2(:,4)=-4*b(3)*(L(:,1).^3)-24*b(3)*(L(:,1).^2).*(L(:,2))+...
                                                12*(b(2)-b(3))*(L(:,1).^2).*(L(:,3))+2*(3*b(1)+15*r(2,1)*b(2))*(L(:,1)).*(L(:,3).^2)-...
                                                4*(3*b(1)+15*r(3,1)*b(3))*(L(:,1)).*L(:,3).*(L(:,2));
                                            N2(:,5)= 4*(b(2)-b(3))*(L(:,1).^3)+2*(3*b(1)+15*r(2,1)*b(2))*(L(:,1).^2).*(L(:,3))-...
                                                2*(3*b(1)+15*r(3,1)*b(3))*(L(:,1).^2).*(L(:,2));
                                            N2(:,6)=4*b(2)*(L(:,1).^3)+24*b(2)*(L(:,1).^2).*(L(:,3))+...
                                                12*(b(2)-b(3))*(L(:,1).^2).*(L(:,2))+4*(3*b(1)+15*r(2,1)*b(2))*(L(:,1)).*(L(:,2)).*(L(:,3))-...
                                                2*(3*b(1)+15*r(3,1)*b(3))*(L(:,1)).*(L(:,2).^2);
                                            
                                            
                                        case 4
                                            
                                            N0=c(3)^2/2*(L(:,1).^3).*(L(:,2).^2)+c(2)^2/2*(L(:,1).^3).*(L(:,3).^2)-c(2)*c(3)*(L(:,1).^3).*(L(:,2)).*(L(:,3))+...
                                                (c(1)*c(2)+5/2*r(2,1)*c(2)^2)*L(:,2).*(L(:,3).^2).*(L(:,1).^2)+...
                                                (c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,3)).*(L(:,2).^2).*(L(:,1).^2);
                                            N1(:,1)=3*c(3)^2/2*(L(:,1).^2).*(L(:,2).^2)+3*c(2)^2/2*(L(:,1).^2).*(L(:,3).^2)-3*c(2)*c(3)*(L(:,1).^2).*(L(:,2)).*(L(:,3))+...
                                                2*(c(1)*c(2)+5/2*r(2,1)*c(2)^2)*L(:,2).*(L(:,3).^2).*(L(:,1))+...
                                                2*(c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,3)).*(L(:,2).^2).*(L(:,1));
                                            N1(:,2)=2*c(3)^2/2*(L(:,1).^3).*(L(:,2))-c(2)*c(3)*(L(:,1).^3).*(L(:,3))+...
                                                (c(1)*c(2)+5/2*r(2,1)*c(2)^2)*(L(:,3).^2).*(L(:,1).^2)+...
                                                2*(c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,3)).*(L(:,2)).*(L(:,1).^2);
                                            N1(:,3)=2*c(2)^2/2*(L(:,1).^3).*(L(:,3))-c(2)*c(3)*(L(:,1).^3).*(L(:,2))+...
                                                2*(c(1)*c(2)+5/2*r(2,1)*c(2)^2)*L(:,2).*(L(:,3)).*(L(:,1).^2)+...
                                                (c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,2).^2).*(L(:,1).^2);
                                            
                                            N2(:,1)=6*c(3)^2/2*(L(:,1)).*(L(:,2).^2)+6*c(2)^2/2*(L(:,1)).*(L(:,3).^2)-6*c(2)*c(3)*(L(:,1)).*(L(:,2)).*(L(:,3))+...
                                                2*(c(1)*c(2)+5/2*r(2,1)*c(2)^2)*L(:,2).*(L(:,3).^2)+...
                                                2*(c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,3)).*(L(:,2).^2);
                                            N2(:,2)=2*c(3)^2/2*(L(:,1).^3)+2*(c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,3)).*(L(:,1).^2);
                                            N2(:,3)=2*c(2)^2/2*(L(:,1).^3)+2*(c(1)*c(2)+5/2*r(2,1)*c(2)^2)*L(:,2).*(L(:,1).^2);
                                            N2(:,4)=6*c(3)^2/2*(L(:,1).^2).*(L(:,2))-3*c(2)*c(3)*(L(:,1).^2).*(L(:,3))+...
                                                2*(c(1)*c(2)+5/2*r(2,1)*c(2)^2)*(L(:,3).^2).*(L(:,1))+...
                                                4*(c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,3)).*(L(:,2)).*(L(:,1));
                                            N2(:,5)=-c(2)*c(3)*(L(:,1).^3)+ 2*(c(1)*c(2)+5/2*r(2,1)*c(2)^2)*(L(:,3)).*(L(:,1).^2)+...
                                                2*(c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,2)).*(L(:,1).^2);
                                            N2(:,6)=6*c(2)^2/2*(L(:,1).^2).*(L(:,3))-3*c(2)*c(3)*(L(:,1).^2).*(L(:,2))+...
                                                4*(c(1)*c(2)+5/2*r(2,1)*c(2)^2)*L(:,2).*(L(:,3)).*(L(:,1))+...
                                                2*(c(1)*c(3)+5/2*r(3,1)*c(3)^2)*(L(:,2).^2).*(L(:,1));
                                            
                                            
                                            
                                        case 5
                                            
                                            N0=-b(3)*c(3)*(L(:,1).^3).*(L(:,2).^2)-b(2)*c(2)*(L(:,1).^3).*(L(:,3).^2)+...
                                                (b(2)*c(3)+b(3)*c(2))*(L(:,1).^3).*L(:,2).*L(:,3)-...
                                                (b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*L(:,2).*(L(:,3).^2).*(L(:,1).^2)...
                                                -(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*L(:,3).*(L(:,2).^2).*(L(:,1).^2);
                                            N1(:,1)=-3*b(3)*c(3)*(L(:,1).^2).*(L(:,2).^2)-3*b(2)*c(2)*(L(:,1).^2).*(L(:,3).^2)+...
                                                3*(b(2)*c(3)+b(3)*c(2))*(L(:,1).^2).*L(:,2).*L(:,3)-...
                                                2*(b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*L(:,2).*(L(:,3).^2).*(L(:,1))...
                                                -2*(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*L(:,3).*(L(:,2).^2).*(L(:,1));
                                            N1(:,2)=-2*b(3)*c(3)*(L(:,1).^3).*(L(:,2))+...
                                                (b(2)*c(3)+b(3)*c(2))*(L(:,1).^3).*L(:,3)-...
                                                (b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*(L(:,3).^2).*(L(:,1).^2)...
                                                -2*(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*L(:,3).*(L(:,2)).*(L(:,1).^2);
                                            N1(:,3)=-2*b(2)*c(2)*(L(:,1).^3).*(L(:,3))+...
                                                (b(2)*c(3)+b(3)*c(2))*(L(:,1).^3).*L(:,2)-...
                                                2*(b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*L(:,2).*(L(:,3)).*(L(:,1).^2)...
                                                -(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*(L(:,2).^2).*(L(:,1).^2);
                                            
                                            N2(:,1)=-6*b(3)*c(3)*(L(:,1)).*(L(:,2).^2)-6*b(2)*c(2)*(L(:,1)).*(L(:,3).^2)+...
                                                6*(b(2)*c(3)+b(3)*c(2))*(L(:,1)).*L(:,2).*L(:,3)-...
                                                2*(b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*L(:,2).*(L(:,3).^2)...
                                                -2*(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*L(:,3).*(L(:,2).^2);
                                            N2(:,2)=-2*b(3)*c(3)*(L(:,1).^3)-2*(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*L(:,3).*(L(:,1).^2);
                                            N2(:,3)=-2*b(2)*c(2)*(L(:,1).^3)-2*(b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*L(:,2).*(L(:,1).^2);
                                            N2(:,4)=-6*b(3)*c(3)*(L(:,1).^2).*(L(:,2))+...
                                                3*(b(2)*c(3)+b(3)*c(2))*(L(:,1).^2).*L(:,3)-...
                                                2*(b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*(L(:,3).^2).*(L(:,1))...
                                                -4*(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*L(:,3).*(L(:,2)).*(L(:,1));
                                            N2(:,5)=(b(2)*c(3)+b(3)*c(2))*(L(:,1).^3)- 2*(b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*(L(:,3)).*(L(:,1).^2)-...
                                                2*(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*(L(:,2)).*(L(:,1).^2);
                                            N2(:,6)=-6*b(2)*c(2)*(L(:,1).^2).*(L(:,3))+...
                                                3*(b(2)*c(3)+b(3)*c(2))*(L(:,1).^2).*L(:,2)-...
                                                4*(b(1)*c(2)+b(2)*c(1)+5*r(2,1)*b(2)*c(2))*L(:,2).*(L(:,3)).*(L(:,1))...
                                                -2*(b(1)*c(3)+b(3)*c(1)+5*r(3,1)*b(3)*c(3))*(L(:,2).^2).*(L(:,1));
                                            
                                            
                                        case 6
                                            
                                            N0=b(3)^2/2*(L(:,1).^3).*(L(:,2).^2)+b(2)^2/2*(L(:,1).^3).*(L(:,3).^2)-...
                                                b(2)*b(3)*(L(:,1).^3).*(L(:,2)).*(L(:,3))+...
                                                (b(1)*b(2)+5/2*r(2,1)*b(2)^2)*L(:,2).*(L(:,3).^2).*(L(:,1).^2)+(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*L(:,3).*(L(:,2).^2).*(L(:,1).^2);
                                            N1(:,1)=3*b(3)^2/2*(L(:,1).^2).*(L(:,2).^2)+3*b(2)^2/2*(L(:,1).^2).*(L(:,3).^2)-...
                                                3*b(2)*b(3)*(L(:,1).^2).*(L(:,2)).*(L(:,3))+...
                                                2*(b(1)*b(2)+5/2*r(2,1)*b(2)^2)*L(:,2).*(L(:,3).^2).*(L(:,1))+2*(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*L(:,3).*(L(:,2).^2).*(L(:,1));
                                            N1(:,2)=2*b(3)^2/2*(L(:,1).^3).*(L(:,2))-b(2)*b(3)*(L(:,1).^3).*(L(:,3))+...
                                                (b(1)*b(2)+5/2*r(2,1)*b(2)^2)*(L(:,3).^2).*(L(:,1).^2)+2*(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*L(:,3).*(L(:,2)).*(L(:,1).^2);
                                            N1(:,3)=2*b(2)^2/2*(L(:,1).^3).*(L(:,3))-b(2)*b(3)*(L(:,1).^3).*(L(:,2))+...
                                                2*(b(1)*b(2)+5/2*r(2,1)*b(2)^2)*L(:,2).*(L(:,3)).*(L(:,1).^2)+(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*(L(:,2).^2).*(L(:,1).^2);
                                            N2(:,1)=6*b(3)^2/2*(L(:,1)).*(L(:,2).^2)+6*b(2)^2/2*(L(:,1)).*(L(:,3).^2)-...
                                                6*b(2)*b(3)*(L(:,1)).*(L(:,2)).*(L(:,3))+...
                                                2*(b(1)*b(2)+5/2*r(2,1)*b(2)^2)*L(:,2).*(L(:,3).^2)+2*(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*L(:,3).*(L(:,2).^2);
                                            N2(:,2)=2*b(3)^2/2*(L(:,1).^3)+2*(b(1)*b(3)+5/2*r(3,1)*b(3)^2)*L(:,3).*(L(:,1).^2);
                                            N2(:,3)=2*b(2)^2/2*(L(:,1).^3)+ 2*(b(1)*b(2)+5/2*r(2,1)*b(2)^2)*L(:,2).*(L(:,1).^2);
                                            N2(:,4)=6*b(3)^2/2*(L(:,1).^2).*(L(:,2))-3*b(2)*b(3)*(L(:,1).^2).*(L(:,3))+...
                                                2*(b(1)*b(2)+5/2*r(2,1)*b(2)^2)*(L(:,3).^2).*(L(:,1))+4*(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*L(:,3).*(L(:,2)).*(L(:,1));
                                            N2(:,5)=-b(2)*b(3)*(L(:,1).^3)+2*(b(1)*b(2)+5/2*r(2,1)*b(2)^2)*(L(:,3)).*(L(:,1).^2)+2*(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*(L(:,2)).*(L(:,1).^2);
                                            N2(:,6)=6*b(2)^2/2*(L(:,1).^2).*(L(:,3))-3*b(2)*b(3)*(L(:,1).^2).*(L(:,2))+...
                                                4*(b(1)*b(2)+5/2*r(2,1)*b(2)^2)*L(:,2).*(L(:,3)).*(L(:,1))+2*(b(1)*b(3)+...
                                                5/2*r(3,1)*b(3)^2)*(L(:,2).^2).*(L(:,1));
                                            
                                            
                                    end
                                    N(:,6*(in-1)+ij)=N0;
                                    Nx(:,6*(in-1)+ij)=N1*a1(:,1);
                                    Ny(:,6*(in-1)+ij)=N1*a1(:,2);
                                    Nxx(:,6*(in-1)+ij)=N2*a2(:,1);
                                    Nyy(:,6*(in-1)+ij)=N2*a2(:,2);
                                    Nxy(:,6*(in-1)+ij)=N2*a2(:,3);
                                end
                            end
                            
                            uxi=0*xig;
                            uyi=0*xig;
                            exxi=0*xig;
                            eyyi=0*xig;
                            uxyi=0*xig;
                            uyxi=0*xig;
                            uxxxi=0*xig;
                            uxxyi=0*xig;
                            uxyyi=0*xig;
                            uyxxi=0*xig;
                            uyxyi=0*xig;
                            uyyyi=0*xig;
                            if any(face_elts==inbox(ie))
                                [crack,~]=GetSignedDistanceToCrack(xyc,xig+1i*yig);
                                for iii=1:2
                                    econn=nconn{iii};
                                    inods=econn(face_elts==inbox(ie),:);
                                    uxn=[U((1:6)+6*(inods(1)-1));U((1:6)+6*(inods(2)-1));U((1:6)+6*(inods(3)-1))];
                                    uyn=[U((1:6)+6*(inods(1)-1)+size(U,1)/2);U((1:6)+6*(inods(2)-1)+size(U,1)/2);U((1:6)+6*(inods(3)-1)+size(U,1)/2)];
                                    switch iii
                                        case 1
                                            mask=crack>0;
                                        case 2
                                            mask=crack<0;
                                    end
                                    mask=diag(sparse(mask));
                                    uxi=uxi+mask*(N*uxn);
                                    uyi=uyi+mask*(N*uyn);
                                    exxi=exxi+mask*(Nx*uxn);
                                    eyyi=eyyi+mask*(Ny*uyn);
                                    uxyi=uxyi+mask*(Ny*uxn);
                                    uyxi=uyxi+mask*(Nx*uyn);
                                    uxxxi=uxxxi+mask*(Nxx*uxn);
                                    uxxyi=uxxyi+mask*(Nxy*uxn);
                                    uxyyi=uxyyi+mask*(Nyy*uxn);
                                    uyxxi=uyxxi+mask*(Nxx*uyn);
                                    uyxyi=uyxyi+mask*(Nxy*uyn);
                                    uyyyi=uyyyi+mask*(Nyy*uyn);
                                end
                            else
                                uxn=[U((1:6)+6*(inods(1)-1));U((1:6)+6*(inods(2)-1));U((1:6)+6*(inods(3)-1))];
                                uyn=[U((1:6)+6*(inods(1)-1)+size(U,1)/2);U((1:6)+6*(inods(2)-1)+size(U,1)/2);U((1:6)+6*(inods(3)-1)+size(U,1)/2)];
                                uxi=N*uxn;
                                uyi=N*uyn;
                                exxi=Nx*uxn;
                                eyyi=Ny*uyn;
                                uxyi=Ny*uxn;
                                uyxi=Nx*uyn;
                                uxxxi=Nxx*uxn;
                                uxxyi=Nxy*uxn;
                                uxyyi=Nyy*uxn;
                                uyxxi=Nxx*uyn;
                                uyxyi=Nxy*uyn;
                                uyyyi=Nyy*uyn;
                            end
                            ux(wig)=uxi;
                            uy(wig)=uyi;
                            exx(wig)=exxi;
                            eyy(wig)=eyyi;
                            uxy(wig)=uxyi;
                            uyx(wig)=uyxi;
                            uxxx(wig)=uxxxi;
                            uxxy(wig)=uxxyi;
                            uxyy(wig)=uxyyi;
                            uyxx(wig)=uyxxi;
                            uyxy(wig)=uyxyi;
                            uyyy(wig)=uyyyi;
                        end
                    end
                    
                    
                    sxx=HH(1,1)*exx+HH(1,2)*eyy+0.5*HH(1,3)*(uxy+uyx)*sqrt(2);
                    syy=HH(2,1)*exx+HH(2,2)*eyy+0.5*HH(2,3)*(uxy+uyx)*sqrt(2);
                    sxy=(HH(3,1)*exx+HH(3,2)*eyy+0.5*HH(3,3)*(uxy+uyx)*sqrt(2))/sqrt(2);
                    sxxx=AA(1,1)*uxxx+AA(1,2)*uxxy+(0.5*AA(1,3))*(uxxy+uyxx)+(0.5*AA(1,4))*(uxyy+uyxy)+AA(1,5)*uyxy+AA(1,6)*uyyy;
                    sxxy=AA(2,1)*uxxx+AA(2,2)*uxxy+(0.5*AA(2,3))*(uxxy+uyxx)+(0.5*AA(2,4))*(uxyy+uyxy)+AA(2,5)*uyxy+AA(2,6)*uyyy;
                    sxyx=AA(3,1)*uxxx+AA(3,2)*uxxy+(0.5*AA(3,3))*(uxxy+uyxx)+(0.5*AA(3,4))*(uxyy+uyxy)+AA(3,5)*uyxy+AA(3,6)*uyyy;
                    sxyy=AA(4,1)*uxxx+AA(4,2)*uxxy+(0.5*AA(4,3))*(uxxy+uyxx)+(0.5*AA(4,4))*(uxyy+uyxy)+AA(4,5)*uyxy+AA(4,6)*uyyy;
                    syyx=AA(5,1)*uxxx+AA(5,2)*uxxy+(0.5*AA(5,3))*(uxxy+uyxx)+(0.5*AA(5,4))*(uxyy+uyxy)+AA(5,5)*uyxy+AA(5,6)*uyyy;
                    syyy=AA(6,1)*uxxx+AA(6,2)*uxxy+(0.5*AA(6,3))*(uxxy+uyxx)+(0.5*AA(6,4))*(uxyy+uyxy)+AA(6,5)*uyxy+AA(6,6)*uyyy;
                    
                    
                    
%                     figure
%                     scatter(xgglo,ygglo,[],uy)
%                     hold on
%                     axis equal
%                     triplot(conng,xi,yi)
%                     plot(xyc(:,1),xyc(:,2),'k')
%                     
%                     figure
%                     scatter(xgglo,ygglo,[],eyy)
%                     hold on
%                     axis equal
%                     triplot(conng,xi,yi)
%                     plot(xyc(:,1),xyc(:,2),'k')
%                     pause
                    
                    %                                         qx=-sign(xg).*((abs(yg)<=abs(xg)).*(abs(xg)>=dfac*dzone))*(1/(dzone*(1-dfac)));
                    %                                         qy=-sign(yg).*((abs(yg)>abs(xg)).*(abs(yg)>=dfac*dzone))*(1/(dzone*(1-dfac)));
                    %
                    %
                    %                                         figure
                    %                                         subplot(2,2,1)
                    %                                         scatter(xgglo,ygglo,[],qx)
                    %                                         colorbar
                    %                                         title('qx')
                    %                                         subplot(2,2,2)
                    %
                    %                                         scatter(xgglo,ygglo,[],qy)
                    %                                         colorbar
                    %                                         title('qy')
                    %
                    %
                    %                                         %TMP=IROT*GQ
                    %                                         tmpx=irot(1,1)*qx+irot(1,2)*qy;
                    %                                         tmpy=irot(2,1)*qx+irot(2,2)*qy;
                    %                                         %GQ=TMP*IROT'
                    %
                    %                                         qxx=irot(1,1)*tmpx;
                    %                                         qyx=irot(2,1)*tmpx;
                    %                                         qxy=irot(1,1)*tmpy;
                    %                                         qyy=irot(2,1)*tmpy;
                    %
                    %                                         figure
                    %                                         subplot(2,2,1)
                    %                                         scatter(xgglo,ygglo,[],qxx)
                    %                                         colorbar
                    %                                         title('qxx ref')
                    %                                         subplot(2,2,2)
                    %                                         scatter(xgglo,ygglo,[],qxy)
                    %                                         colorbar
                    %                                         title('qxy ref')
                    %                                         subplot(2,2,3)
                    %                                         scatter(xgglo,ygglo,[],qyx)
                    %                                         colorbar
                    %                                         title('qyx ref')
                    %                                         subplot(2,2,4)
                    %                                         scatter(xgglo,ygglo,[],qyy)
                    %                                         colorbar
                    %                                         title('qyy ref')
                    
                    
                    
                    
                    
                    
                    
%                    qtx=-sign(xg).*((abs(yg)<=abs(xg)).*(abs(xg)>=dfac*dzone))*(1/(dzone*(1-dfac)));
%                    qty=-sign(yg).*((abs(yg)>abs(xg)).*(abs(yg)>=dfac*dzone))*(1/(dzone*(1-dfac)));
zg=(xg+1i*yg);
                    qtx=-cos(angle(zg)).*((abs(zg)>=dfac*dzone)).*((abs(zg)<=dzone))*(1/(dzone*(1-dfac)));
                    qty=-sin(angle(zg)).*((abs(zg)>=dfac*dzone)).*((abs(zg)<=dzone))*(1/(dzone*(1-dfac)));
                    

                    tmpx=irot(1,1)*qtx+irot(1,2)*qty;
                    tmpy=irot(2,1)*qtx+irot(2,2)*qty;
%                     
%                     figure
%                     quiver(xgglo,ygglo,tx,ty)
%                     axis equal
                    
%                                                             figure
%                                                             subplot(2,2,1)
%                                                             scatter(xgglo,ygglo,[],qtx)
%                                                             colorbar
%                                                             title('qtx loc')
%                                                             subplot(2,2,2)
%                                                             scatter(xgglo,ygglo,[],qty)
%                                                             colorbar
%                                                             title('qty loc')
%                                                             subplot(2,2,3)
%                                                             scatter(xgglo,ygglo,[],tmpx)
%                                                             colorbar
%                                                             title('qtx glo')
%                                                             subplot(2,2,4)
%                                                             scatter(xgglo,ygglo,[],tmpy)
%                                                             colorbar
%                                                             title('qty glo')
                    
                    qxx=tx.*tmpx;
                    qyx=ty.*tmpx;
                    qxy=tx.*tmpy;
                    qyy=ty.*tmpy;
                    
                    
                    
%                                                             figure
%                                                             subplot(2,2,1)
%                                                             scatter(xgglo,ygglo,[],qxx)
%                                                             colorbar
%                                                             title('qxx glo2')
%                                                             subplot(2,2,2)
%                                                             scatter(xgglo,ygglo,[],qxy)
%                                                             colorbar
%                                                             title('qxy glo2')
%                                                             subplot(2,2,3)
%                                                             scatter(xgglo,ygglo,[],qyx)
%                                                             colorbar
%                                                             title('qyx glo2')
%                                                             subplot(2,2,4)
%                                                             scatter(xgglo,ygglo,[],qyy)
%                                                             colorbar
%                                                             title('qyy glo2')
                    
                    
                    
                    
                    
                    wdetJ(isnan(exx))=0;
                    welas=0.5*(sxx.*exx+syy.*eyy+(uxy+uyx).*sxy+sxxx.*uxxx+syyx.*uyxy+sxyx.*(uxxy+uyxx)+sxxy.*uxxy+syyy.*uyyy+sxyy.*(uxyy+uyxy));
                    % welas*qx,x - (sxx*ux,x+syx*uy,x)*qx,x - (sxx*ux,y+syx*uy,y)*qy,x
                    Jx=(welas-sxx.*exx-sxy.*uyx).*qxx-(sxx.*uxy+sxy.*eyy).*qyx;
                    Jgx=-(sxxx.*uxxx+sxyx.*(uxxy+uyxx)+syyx.*uyxy).*qxx-(sxxx.*uxxy+sxyx.*(uxyy+uyxy)+syyx.*uyyy).*qyx;
                    Jgx=Jgx-(sxxx.*uxxx+sxyx.*uxxy+syyx.*uyxy).*qxx-(sxxx.*uxxy+sxyx.*uxyy+syyx.*uyyy).*qyx;
                    % welas*qy,y - (sxy*ux,x+syy*uy,x)*qx,y - (sxy*ux,y+syy*uy,y)*qy,y
                    Jy=-(sxy.*exx+syy.*uyx).*qxy+(welas-sxy.*uxy-syy.*eyy).*qyy;
                    Jgy=-(sxxy.*uxxx+sxyy.*(uxxy+uyxx)+syyy.*uyxy).*qxy-(sxxy.*uxxy+sxyy.*(uxyy+uyxy)+syyy.*uyyy).*qyy;
                    Jgy=Jgy-(sxxy.*uxxx+sxyy.*uxxy+syyy.*uyxy).*qxy-(sxxy.*uxxy+sxyy.*uxyy+syyy.*uyyy).*qyy;
                    Jx(isnan(exx))=0;
                    Jy(isnan(exx))=0;
                    Jgx(isnan(exx))=0;
                    Jgy(isnan(exx))=0;
                    
                    k{2*(iz-1)+itip,3}=-((wdetJ'*(Jx+Jy))*pix2m+(wdetJ'*(Jgx+Jgy))/pix2m);
                    
                    
%                              exxo=exx;
%                             eyyo=eyy;
%                             uxyo=uxy;
%                             uyxo=uyx;
%                             uxxxo=uxxx;
%                             uxxyo=uxxy;
%                             uxyyo=uxyy;
%                             uyxxo=uyxx;
%                             uyxyo=uyxy;
%                             uyyyo=uyyy;
%                    
% % symmetrique part uy=1/2*(uy(y)-uy(-y))
% %                  ux=1/2*(ux(y)+ux(-y))
%                             exx=0.5*(exxo+flipud(exxo));
%                             eyy=0.5*(eyyo+flipud(eyyo));
%                             uxy=0.5*(uxyo+flipud(uxyo));
%                             uyx=0.5*(uyxo+flipud(uyxo));
%                             uxxx=0.5*(uxxxo+flipud(uxxxo));
%                             uxxy=0.5*(uxxyo+flipud(uxxyo));
%                             uxyy=0.5*(uxyyo+flipud(uxyyo));
%                             uyxx=0.5*(uyxxo+flipud(uyxxo));
%                             uyxy=0.5*(uyxyo+flipud(uyxyo));
%                             uyyy=0.5*(uyyyo+flipud(uyyyo));
%                   
%                     
%                     sxx=HH(1,1)*exx+HH(1,2)*eyy+0.5*HH(1,3)*(uxy+uyx)*sqrt(2);
%                     syy=HH(2,1)*exx+HH(2,2)*eyy+0.5*HH(2,3)*(uxy+uyx)*sqrt(2);
%                     sxy=(HH(3,1)*exx+HH(3,2)*eyy+0.5*HH(3,3)*(uxy+uyx)*sqrt(2))/sqrt(2);
%                     sxxx=AA(1,1)*uxxx+AA(1,2)*uxxy+(0.5*AA(1,3))*(uxxy+uyxx)+(0.5*AA(1,4))*(uxyy+uyxy)+AA(1,5)*uyxy+AA(1,6)*uyyy;
%                     sxxy=AA(2,1)*uxxx+AA(2,2)*uxxy+(0.5*AA(2,3))*(uxxy+uyxx)+(0.5*AA(2,4))*(uxyy+uyxy)+AA(2,5)*uyxy+AA(2,6)*uyyy;
%                     sxyx=AA(3,1)*uxxx+AA(3,2)*uxxy+(0.5*AA(3,3))*(uxxy+uyxx)+(0.5*AA(3,4))*(uxyy+uyxy)+AA(3,5)*uyxy+AA(3,6)*uyyy;
%                     sxyy=AA(4,1)*uxxx+AA(4,2)*uxxy+(0.5*AA(4,3))*(uxxy+uyxx)+(0.5*AA(4,4))*(uxyy+uyxy)+AA(4,5)*uyxy+AA(4,6)*uyyy;
%                     syyx=AA(5,1)*uxxx+AA(5,2)*uxxy+(0.5*AA(5,3))*(uxxy+uyxx)+(0.5*AA(5,4))*(uxyy+uyxy)+AA(5,5)*uyxy+AA(5,6)*uyyy;
%                     syyy=AA(6,1)*uxxx+AA(6,2)*uxxy+(0.5*AA(6,3))*(uxxy+uyxx)+(0.5*AA(6,4))*(uxyy+uyxy)+AA(6,5)*uyxy+AA(6,6)*uyyy;
 
                    
                    
                    anglon=(-pi/2:pi/100:pi/2)';
                    xon=htip*cos(anglon);
                    yon=htip*sin(anglon);
                    sxxon=griddata(xg,yg,sxx,xon,yon,'nearest');
                    syyon=griddata(xg,yg,syy,xon,yon,'nearest');
                    sxyon=griddata(xg,yg,sxy,xon,yon,'nearest');
                    sxn=sxxon.*cos(anglon)+sxyon.*sin(anglon);
                    syn=sxyon.*cos(anglon)+syyon.*sin(anglon);
                    snn=sxn.*cos(anglon)+syn.*sin(anglon);
                    stn=-sxn.*sin(anglon)+syn.*cos(anglon);
                    figure
                    plot(anglon,snn,'b')
                    hold on
                    plot(anglon,abs(stn),'g')
                    plot(anglon,sxxon,'r')
                    plot(anglon,syyon,'k')
                    plot(anglon,sxyon,'m')
                    legend('snn','snt','sxx','syy','sxy')

                    w=exp(-(2*abs(xg+1i*yg)/(htip)).^2);%.*(xg>=0);
                    nw=wdetJ'*w;
                    w=diag(sparse(w));
                    Et=(wdetJ'*(w*[exx,uxy,uyx,eyy]))'/nw
                    St=(wdetJ'*(w*[sxx,sxy,sxy,syy]))'/nw
                    k{2*(iz-1)+itip,1}=[Et,St];
                    
                    EEt=(wdetJ'*(w*[uxxx,uxxy,uxyy,uyxx,uyxy,uyyy]))'/nw
                    HSt=(wdetJ'*(w*[sxxx,sxxy,sxyy,sxyx,syyx,syyy]))'/nw
                    k{2*(iz-1)+itip,2}=[EEt,HSt];
                    
%                     w=exp(-(2*abs(xg+1i*yg)/(htip)).^2);%.*(xg>=0);
%                     nw=wdetJ'*w;
%                     w=diag(sparse(w));
%                     Et=(wdetJ'*(w*[exx,uxy,uyx,eyy]))'/nw
%                     St=(wdetJ'*(w*[sxx,sxy,sxy,syy]))'/nw
%                     k{2*(iz-1)+itip,1}=[Et,St];
%                     
%                     EEt=(wdetJ'*(w*[uxxx,uxxy,uxyy,uyxx,uyxy,uyyy]))'/nw
%                     HSt=(wdetJ'*(w*[sxxx,sxxy,sxyy,sxyx,syyx,syyy]))'/nw
%                     k{2*(iz-1)+itip,2}=[EEt,HSt];
                    
                end
            end
            
    end
end


