function [k]=SIF(nmod,model,U,xn,yn,conng)
k=zeros(2*size(model.zone,2),3);
matmod=model.material_parameters;
mu=0.5*matmod.young/(1+matmod.nu);
lambda=matmod.nu*matmod.young/((1.+matmod.nu)*(1.-2.*matmod.nu));
dfac=0.25;
kappa=3-4*matmod.nu;
Es=matmod.young/(1-(matmod.nu)^2);
scal=2*mu*sqrt(2*pi);

%
for iz=1:size(model.zone,2)
    zone=model.zone(:,iz);
    switch zone{4}
        case 5 %CRACK
            dzone=zone{8};
            xyc=zone{2};
            if dzone>0
                face_elts=zone{10};
                nconn=zone{12};
                
                xoo=[-dzone;dzone;dzone;-dzone]*1.25;
                yoo=[-dzone;-dzone;dzone;dzone]*1.25;
                xo=[-dzone;dzone;dzone;-dzone];
                yo=[-dzone;-dzone;dzone;dzone];
                conn=1:4;
                elt=4;
                rint=0;
                ng=0;
                ns=[1,1]*50;
                Nnodes=[length(xo),1,1];
                Nelems=[length(elt),1,1];
                zo=1;
                jmesh_file=fullfile('TMP','jmesh');
                save(jmesh_file,'rint','xo','yo','zo','Nnodes','Nelems','conn','ng','ns','elt')
                [wdetJ,inde]=GetWeigthDetJ(jmesh_file,[1,1],1,'sub_cells');
                wdetJ=diag(wdetJ);
                phig=CreateFiniteElementBasis(jmesh_file,[1,1],1,[],'sub_cells');
                xg=phig*xo;
                yg=phig*yo;
                
                for itip=1:1+(zone{9}>0)
                    exx=0*xg+NaN;
                    eyy=0*xg+NaN;
                    uxy=0*xg+NaN;
                    uyx=0*xg+NaN;
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
                    xglo=xn-xytip(1);
                    yglo=yn-xytip(2);
                    xgglo=mean(xglo(conng),2);
                    ygglo=mean(yglo(conng),2);
                    xgloc=xgglo*t(1)+ygglo*t(2);
                    ygloc=xgglo*n(1)+ygglo*n(2);
                    inbox=find(inpolygon(xgloc,ygloc,xoo([1:4,1]),yoo([1:4,1])));
                    cin=conng(inbox,:);
                    for ie=1:size(cin,1)
                        inods=cin(ie,:);
                        xi=xglo(inods)*t(1)+yglo(inods)*t(2);
                        yi=xglo(inods)*n(1)+yglo(inods)*n(2);
                        
                        [xig,yig,wig]=GetGaussPointsTriangle(0,1,xi,yi,xg,yg);
                        if any(wig>0)
                            xig=xig(wig);
                            yig=yig(wig);
                            N_r=[-1+0*xig,1+0*xig,0*yig];
                            N_s=[-1+0*yig,0*xig,1+0*yig];
                            dxdr=N_r*xi;
                            dydr=N_r*yi;
                            dxds=N_s*xi;
                            dyds=N_s*yi;
                            detJ=(dxdr.*dyds-dydr.*dxds);
                            invJ=[dyds./detJ,-dxds./detJ,-dydr./detJ,dxdr./detJ];
                            exxi=0*xig;
                            eyyi=0*xig;
                            uxyi=0*xig;
                            uyxi=0*xig;
                            if any(face_elts==inbox(ie))
                                N=[1-xig-yig,xig,yig];
                                xigglo=N*xn(inods);
                                yigglo=N*yn(inods);
                                [crack,~]=GetSignedDistanceToCrack(xyc,xigglo+1i*yigglo);
                                %                                 figure
                                %                                 plot(xg,yg,'kx')
                                %                                 hold on
                                %                                 plot(xi,yi,'b+')
                                %                                 plot(N*xi,N*yi,'ro')
                                for ii=1:2
                                    econn=nconn{ii};
                                    inods=econn(face_elts==inbox(ie),:);
                                    uxi=U(inods)*t(1)+U(inods+size(U,1)/2)*t(2);
                                    uyi=U(inods)*n(1)+U(inods+size(U,1)/2)*n(2);
                                    switch ii
                                        case 1
                                            mask=crack>0;
                                        case 2
                                            mask=crack<0;
                                    end
                                    mask=diag(sparse(mask));
                                    for in=1:3
                                        N_x=N_r(:,in).*invJ(:,1)+N_s(:,in).*invJ(:,3);
                                        N_y=N_r(:,in).*invJ(:,2)+N_s(:,in).*invJ(:,4);
                                        exxi=exxi+mask*(N_x*uxi(in));
                                        eyyi=eyyi+mask*(N_y*uyi(in));
                                        uxyi=uxyi+mask*(N_y*uxi(in));
                                        uyxi=uyxi+mask*(N_x*uyi(in));
                                    end
                                end
                            else
                                uxi=U(inods)*t(1)+U(inods+size(U,1)/2)*t(2);
                                uyi=U(inods)*n(1)+U(inods+size(U,1)/2)*n(2);
                                for in=1:3
                                    N_x=N_r(:,in).*invJ(:,1)+N_s(:,in).*invJ(:,3);
                                    N_y=N_r(:,in).*invJ(:,2)+N_s(:,in).*invJ(:,4);
                                    exxi=exxi+N_x*uxi(in);
                                    eyyi=eyyi+N_y*uyi(in);
                                    uxyi=uxyi+N_y*uxi(in);
                                    uyxi=uyxi+N_x*uyi(in);
                                end
                            end
                            exx(wig)=exxi;
                            eyy(wig)=eyyi;
                            uxy(wig)=uxyi;
                            uyx(wig)=uyxi;
                        end
                    end
                    
                    ltr=lambda*(exx+eyy);
                    sxx=2*mu*exx+ltr;
                    syy=2*mu*eyy+ltr;
                    sxy=mu*(uxy+uyx);
                    
                    
                    qx=-sign(xg).*((abs(yg)<=abs(xg)).*(abs(xg)>=dfac*dzone))*(1/(dzone*(1-dfac)));
                    qy=-sign(yg).*((abs(yg)>abs(xg)).*(abs(yg)>=dfac*dzone))*(1/(dzone*(1-dfac)));
                    
                    wdetJ(isnan(exx))=0;
                    welas=0.5*(sxx.*exx+syy.*eyy+(uxy+uyx).*sxy);
                    Jx=(welas-sxx.*exx-sxy.*uyx).*qx;
                    Jy=-(sxy.*exx+syy.*uyx).*qy;
                    Jx(isnan(exx))=0;
                    Jy(isnan(exx))=0;
                    k(2*(iz-1)+itip,3)=-(wdetJ'*(Jx+Jy));
                    heq=1;
                    dist=abs(xg+1i*yg);
                    angl=angle(xg+1i*yg);
                    for im=1:2
                        
                        [harm1,amp1]=KMPotFourrier(im,heq,kappa);
                        
                        ft=0*dist(:);
                        fr=0*dist(:);
                        
                        for kk=1:3
                            fr=fr+amp1(kk)*exp(harm1(kk)*angl(:)*1i/2);
                            ft=ft+harm1(kk)*1i/2*amp1(kk)*exp(harm1(kk)*angl(:)*1i/2);
                        end
                        ft=ft.*(dist(:).^(heq/2-1));
                        fr=heq/2*fr.*(dist(:).^(heq/2-1));
                        
                        dfndth=imag(ft);
                        dfndr=imag(fr);
                        dftdth=real(ft);
                        dftdr=real(fr);
                        
                        drdt=cos(angl(:));
                        drdn=sin(angl(:));
                        dthdt=-sin(angl(:));
                        dthdn=cos(angl(:));
                        
                        dfndt=dfndr.*drdt+dfndth.*dthdt;
                        dfndn=dfndr.*drdn+dfndth.*dthdn;
                        dftdt=dftdr.*drdt+dftdth.*dthdt;
                        dftdn=dftdr.*drdn+dftdth.*dthdn;
                        
                        altr=lambda*(dftdt+dfndn);
                        asxx=2*mu*dftdt+altr;
                        asyy=2*mu*dfndn+altr;
                        asxy=mu*(dfndt+dftdn);
                        
                        welas=0.5*(asxx.*exx+asyy.*eyy+(uxy+uyx).*asxy);
                        Jx=(welas-asxx.*exx-asxy.*uyx-sxx.*dftdt-sxy.*dfndt).*qx;
                        Jy=-(sxy.*dftdt+syy.*dfndt+asxy.*exx+asyy.*uyx).*qy;
                        Jx(isnan(exx))=0;
                        Jy(isnan(exx))=0;
                        k(2*(iz-1)+itip,im)=-Es/(2*scal)*(wdetJ'*(Jx+Jy));
                    end
                    
                end
            end
            
    end
end

