function lvl=GetSignedDistanceToZone(fem_model,roi,xy,y)
if ~isstruct(fem_model)
    xyp=fem_model;
    clear fem_model
    fem_model.zone{1,1}=1;
fem_model.zone{2,1}=xyp;
fem_model.zone{3,1}=0;
fem_model.zone{4,1}=2;
fem_model.zone{5,1}=0;
fem_model.zone{6,1}=0;
end
xy=[xy,y]*[1;1i];
np=length(xy);
if np==1
    xy=repmat(xy,2,1);
end
roi=roi+[-1,1,-1,1];
xyp=roi([1,2,2,1,1])+1i*roi([3,3,4,4,3]);
lio=inpolygon(real(xy),imag(xy),real(xyp),imag(xyp));
li=0;
lvl=Inf;
for ip=1:length(xyp)
    lvl=min(lvl,abs(xy-xyp(ip)));
end

for ip=1:length(xyp)-1
    seg=diff(xyp(ip+(0:1)));
    n=(-1i*seg)/abs(seg);
    t=(seg)/abs(seg);
    d=abs(real((xy-xyp(ip))'*n))';
    along=abs(real((xy-0.5*(xyp(ip)+xyp(ip+1)))'*t))'<0.5*abs(seg);
    lvl(along)=min(lvl(along),d(along));
end

if ~isempty(fem_model.zone)
    %    xy=[xy,y]*[1;1i];
    tzone=fem_model.zone(1,:);
    for iz=1:length(tzone)
        if isempty(tzone{iz}),tzone{iz}=0;end
    end
    
    insideout=cell2mat(tzone);
    [insideout,id]=sort(insideout,'descend');
    for iz=1:length(id)
        zone=fem_model.zone(:,id(iz));
        if ~(zone{1}<0)
            if numel(li)==1
                if any(insideout==1)
                    li=0;
                else
                    li=1;
                end
            end
            switch zone{4}
                case {1,2}
                    lvli=Inf;
                    xyp=zone{2};
                    xyp=xyp*[1;1i];
                    if zone{1}
                        li=li|inpolygon(real(xy),imag(xy),real(xyp),imag(xyp));
                    else
                        [in,on]=inpolygon(real(xy),imag(xy),real(xyp),imag(xyp));
                        li=li&(~in);
                    end
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
                    lvl=min(lvl,lvli);
                    
                case 3
                    xyp=zone{5};
                    lvli=abs(xy-(xyp(1)+1i*xyp(2)))-xyp(3);
                    if zone{1}
                        li=li|(lvli<0);
                    else
                        li=li&(lvli>0);
                    end
                    lvl=min(lvl,abs(lvli));
            end
        end
    end
end
if numel(li)==1
    li=1;
end
li=lio&li;
li=-2*(li>0)+1;
lvl=lvl.*li;
if np==1
    lvl=lvl(1);
end