function [lvl,xyfix,indc]=GetMeshDensity(fem_model,xy,y)
if isfield(fem_model,'phantom_nodes')
    cmesh=fem_model.phantom_nodes==0;
else
    cmesh=1;
end
xyfix=[];indc=[];
dh=0.2;
lvl=ones(size(xy));
if ~isempty(fem_model.zone)
    h=mean(fem_model.mesh_size);
    xy=[xy,y]*[1;1i];
    for iz=1:size(fem_model.zone,2)
        zone=fem_model.zone(:,iz);
        dmin=real(zone{6})/h;
        if zone{4}==1
            xyp=zone{2};
            xyp=xyp*[1;1i];
            xyfix=[xyfix;[real(xyp(1:4)),imag(xyp(1:4))]];
        end
        if dmin>0
            switch zone{4}
                case {1,2,4}
                    xyp=zone{2};
                    xyp=xyp*[1;1i];
                    lvli=Inf;
                    for ip=1:length(xyp)
                        lvli=min(lvli,abs(xy-xyp(ip)));
                    end
                    for ip=1:length(xyp)-1
                        seg=diff(xyp(ip+(0:1)));
                        n=(-1i*seg)/abs(seg);
                        t=(seg)/abs(seg);
                        d=abs(real((xy-xyp(ip))'*n))';
                        along=abs(real((xy-0.5*(xyp(ip)+xyp(ip+1)))'*t))'<0.5*abs(seg);
                        lvli(along)=min(lvli(along),d(along));
                    end
                    if zone{1}==-2
                        li=~inpolygon(real(xy),imag(xy),real(xyp),imag(xyp));
                        lvli=lvli.*li;
                    end
                case 3
                    xyp=zone{5};
                    switch zone{1}
                        case {0,1,-1}
                            lvli=abs(abs(xy-(xyp(1)+1i*xyp(2)))-xyp(3));
                        case -2
                            lvli=max(0,(abs(xy-(xyp(1)+1i*xyp(2)))-xyp(3)));
                    end
                case 5
                    xyp=zone{2};
                    xyp=xyp*[1;1i];
                    rc=zone{8};
                    rc=rc(1)*cmesh;
                    dtip=real(zone{7})/h;
                    lvli=max(0,(abs(xy-xyp(1))-rc));
                    lvli=dtip+(1-dtip)*min(lvli*dh/h,1);
                    lvl=max(lvli,lvl).*(lvli>1)+min(lvli,lvl).*(lvli<=1);
                    if zone{9}
                        lvli=max(0,(abs(xy-xyp(end))-rc));
                        lvli=dtip+(1-dtip)*min(lvli*dh/h,1);
                        lvl=max(lvli,lvl).*(lvli>1)+min(lvli,lvl).*(lvli<=1);
                    end
                    
                    lvli=Inf;
                    for ip=1:length(xyp)
                        lvli=min(lvli,abs(xy-xyp(ip)));
                    end
                    for ip=1:length(xyp)-1
                        seg=diff(xyp(ip+(0:1)));
                        n=(-1i*seg)/abs(seg);
                        t=(seg)/abs(seg);
                        d=abs(real((xy-xyp(ip))'*n))';
                        along=abs(real((xy-0.5*(xyp(ip)+xyp(ip+1)))'*t))'<0.5*abs(seg);
                        lvli(along)=min(lvli(along),d(along));
                    end
                otherwise
                    lvli=0;
            end
            if any(lvli>0)
                lvli=dmin+(1-dmin)*min(lvli*dh/h,1);
                lvl=max(lvli,lvl).*(lvli>1)+min(lvli,lvl).*(lvli<=1);
            end
        end
        if imag(zone{6})
            switch zone{4}
                case {1,2,4}
                    xyp=zone{2};
                    xyp=xyp*[1;1i];
                    xyfixi=xyp;
                    for ip=1:length(xyp)-1
                        seg=diff(xyp(ip+(0:1)));
                        if dmin>0
                            nn=floor(abs(seg)/(dmin*h));
                        else
                            nn=floor(abs(seg)/h);
                        end
                        xyfixi=[xyfixi;xyp(ip)+seg/nn*((1:nn-1)')];
                    end
                    
                case 3
                    xyp=zone{5};
                    if dmin>0
                        dtheta=2*pi/round(2*pi*xyp(3)/(dmin*h));
                    else
                        dtheta=2*pi/round(2*pi*xyp(3)/(h));
                    end
                    theta=(0:dtheta:(2*pi-dtheta))';
                    xyfixi=xyp(1)+xyp(3)*exp(1i*theta)+1i*xyp(2);
                case 5
                    xyp=zone{2};
                    xyp=xyp*[1;1i];
                    xyfixi=[];
                    rc=zone{8};
                    rc=rc(1);
                    dtip=real(zone{7})/h;
                    xyo=xyp(1);
                    for ip=1:length(xyp)-1
                        
                        seg=xyp(ip+1)-xyo(end);
                        lvli=max(0,(abs(xyo(end)-xyp(1))-rc));
                        dx=(dtip+(dmin-dtip)*min(lvli*dh/h,1));
                        if zone{9}
                            lvli=max(0,(abs(xyo(end)-xyp(end))-rc));
                            dx2=(dtip+(dmin-dtip)*min(lvli*dh/h,1));
                            dx=max(dx2,dx).*(dx2>1)+min(dx2,dx).*(dx2<=1);
                        end
                        dx=dx*h;
                        iter=0;
                        while abs(seg)>0.5*dx&&iter<100
                            xyo=[xyo;xyo(end)+dx*seg/abs(seg)];
                            seg=xyp(ip+1)-xyo(end);
                            lvli=max(0,(abs(xyo(end)-xyp(1))-rc));
                            dx=(dtip+(dmin-dtip)*min(lvli*dh/h,1));
                            if zone{9}
                                lvli=max(0,(abs(xyo(end)-xyp(end))-rc));
                                dx2=(dtip+(dmin-dtip)*min(lvli*dh/h,1));
                                dx=max(dx2,dx).*(dx2>1)+min(dx2,dx).*(dx2<=1);
                            end
                            dx=dx*h;
                            iter=iter+1;
                        end
                        xyo(end)=xyp(ip+1);
                        
                        
                    end
                    xyfixi=[xyfixi;xyo];
                    indc{1,iz}=length(xyfix)+(1:length(xyo));
                    indc{2,iz}=[];
                    indc{3,iz}=[];
                    if rc>h*dtip
                        rco=rc;
                        [val,id]=min(abs(abs(xyo-xyp(1))-rco));
                        rc=abs(xyo(id)-xyp(1));
                        if dtip>0
                            dtheta=2*pi/round(2*pi*rc/(dtip*h));
                        else
                            dtheta=2*pi/round(2*pi*rc/(h));
                        end
                        theta=(dtheta:dtheta:(2*pi-dtheta))';
                        seg=diff(xyp(1:2));
                        indc{2,iz}=length(xyfix)+length(xyfixi)+(1:length(theta));
                        xyfixi=[xyfixi;xyp(1)+rc*exp(1i*theta)*seg/abs(seg)];
                        
                        if zone{9}
                            [val,id]=min(abs(abs(xyo-xyp(end))-rco));
                            rc=abs(xyo(id)-xyp(end));
                            if dtip>0
                                dtheta=2*pi/round(2*pi*rc/(dtip*h));
                            else
                                dtheta=2*pi/round(2*pi*rc/(h));
                            end
                            theta=(dtheta:dtheta:(2*pi-dtheta))';
                            seg=-diff(xyp((end-1):end));
                            indc{3,iz}=length(xyfix)+length(xyfixi)+(1:length(theta));
                            xyfixi=[xyfixi;xyp(end)+rc*exp(1i*theta)*seg/abs(seg)];
                        end
                        rc=zone{8};
                        if numel(rc)>1
                            indc{4,iz}=[];
                            indc{5,iz}=[];
                            
                            rco=rc(2);
                            [val,id]=min(abs(abs(xyo-xyp(1))-rco));
                            rc=abs(xyo(id)-xyp(1));
                            if dtip>0
                                dtheta=2*pi/round(2*pi*rc/(dtip*h));
                            else
                                dtheta=2*pi/round(2*pi*rc/(h));
                            end
                            theta=(dtheta:dtheta:(2*pi-dtheta))';
                            seg=diff(xyp(1:2));
                            indc{4,iz}=length(xyfix)+length(xyfixi)+(1:length(theta));
                            xyfixi=[xyfixi;xyp(1)+rc*exp(1i*theta)*seg/abs(seg)];
                            
                            if zone{9}
                                [val,id]=min(abs(abs(xyo-xyp(end))-rco));
                                rc=abs(xyo(id)-xyp(end));
                                if dtip>0
                                    dtheta=2*pi/round(2*pi*rc/(dtip*h));
                                else
                                    dtheta=2*pi/round(2*pi*rc/(h));
                                end
                                theta=(dtheta:dtheta:(2*pi-dtheta))';
                                seg=-diff(xyp((end-1):end));
                                indc{5,iz}=length(xyfix)+length(xyfixi)+(1:length(theta));
                                xyfixi=[xyfixi;xyp(end)+rc*exp(1i*theta)*seg/abs(seg)];
                            end
                            
                            
                        end
                    end
            end
            xyfix=[xyfix;[real(xyfixi),imag(xyfixi)]];
            
        end
    end
    if ~(cmesh==1)
        nodes=[];
        for iz=1:size(indc,2)
            for iy=1:size(indc,1)
                nodes=[nodes,cell2mat(indc(iy,iz))];
            end
        end
        xyfix(nodes(:),:)=[];
        indc=[];
    end
end