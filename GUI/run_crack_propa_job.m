function run_crack_propa_job(param0,model0)
nmod=0;
xos=[];
yos=[];
conns=[];
Us=[];
Kss=[];
das=[];
roi=param0.roi;
param=param0;
model=model0;
KIc=model0.material_parameters.KIc;
filreso=param.result_file;
filres='upropa.res';
param.result_file=filres;
model.mesh_file='upropa.vtk';
copyfile(model0.mesh_file,model.mesh_file);
pix2m=param.pixel_size;
xfem=model.phantom_nodes;
dao=param.da;
modes=[1,2];
if xfem
    ind=1;
else
    ind=0:7;
end
mu=0.5*model.material_parameters.young/(1+model.material_parameters.nu);
Es=model.material_parameters.young/(1-(model.material_parameters.nu)^2);
scal=2*mu*sqrt(2*pi);
scalamp=scal*(pix2m.^(1-(ind)*.5));
found=find(ind==1);
tipso=zeros(size(model.zone,2),2);
nw=length(modes)*length(ind);
ntip=0;
zones=model.zone;
for iz=1:size(zones,2)
    zone=zones(:,iz);
    switch zone{4}
        case 5
            if (zone{8}>0)
                tipso(iz,1)=ntip+1;
                ntip=ntip+1;
                if  (zone{9}>0)
                    tipso(iz,2)=ntip+1;
                    ntip=ntip+1;
                end
            end
        case 6
            model0.zone{7,iz}=[];
    end
end
cracks=find(cell2mat(zones(4,:))==5);
tipso(~(cell2mat(zones(4,:))==5),:)=[];
step=1;
filreso=strrep(filreso,'.res','');
filres=strrep(filres,'.res','');
delete([fullfile('VTK',filreso),'-0*.vtk']);
copyfile([filreso,'.dat'],[filres,'.dat']);
remesh=0*tipso;
dacum=0*tipso;
while 1
    if xfem
        run_fdfea_job(param,model)
    else
        run_fea_job(param,model)
    end
    modeli=model;
    load([filres,'.res'],'-mat','model','U','Ks','tips','xo','yo','conn')
    if step==1
        movefile(fullfile('VTK',[filres , sprintf('-%06d.vtk',0)]),fullfile('VTK',[filreso , sprintf('-%06d.vtk',0)]));
    end
    movefile(fullfile('VTK',[filres , sprintf('-%06d.vtk',1)]),fullfile('VTK',[filreso , sprintf('-%06d.vtk',step)]));
    if any(tips(:)>0)
        xos{step}=xo;
        yos{step}=yo;
        conns{step}=conn;
        lf=0*tips+Inf;
        dasi=0*tips;
        Ksi=zeros(sum(tipso(:)>0)*nw,1);
        for ic=1:size(tipso,1)
            for it=1:2
                if (tipso(ic,it)>=1)&&(tips(ic,it)>=1)
                    if xfem
                        K1=sqrt(pix2m)*Ks(found+2*(cracks(ic)-1)+(it-1),1);
                        K2=sqrt(pix2m)*Ks(found+2*(cracks(ic)-1)+(it-1),2);
                    else
                        K1=scalamp(found)*Ks(found+(tips(ic,it)-1)*nw);
                        K2=scalamp(found)*Ks(found+length(ind)+(tips(ic,it)-1)*nw,:);
                    end
                    [K1eq,Oc]=MaxHoopStressCriterion(K1,K2);
                    lf(ic,it)=KIc/K1eq;
                end
            end
        end
        Us{step}=min(lf)*U;
        Ks=min(lf)*Ks;
        for iz=1:size(zones,2)
            zone=model0.zone(:,iz);
            switch zone{4}
                case 6
                    fu=model0.zone{7,iz};
                    for izz=1:size(model.zone,2)
                        if strcmp(get(model.zone{3,izz},'Tag'),get(model0.zone{3,iz},'Tag'))
                            
                            model0.zone{7,iz}=[fu,min(lf)*model.zone{7,izz}];
                            break
                        end
                    end
            end
        end
        model=modeli;
        for ic=1:size(tipso,1)
            for it=1:2
                if (tipso(ic,it)>=1)&&(tips(ic,it)>=1)
                    if xfem
                        Ksi((1:nw)+(tipso(ic,it)-1)*nw)=Ks(found+2*(cracks(ic)-1)+(it-1),1:2)';
                        K1=sqrt(pix2m)*Ks(found+2*(cracks(ic)-1)+(it-1),1);
                        K2=sqrt(pix2m)*Ks(found+2*(cracks(ic)-1)+(it-1),2);
                    else
                        Ksi((1:nw)+(tipso(ic,it)-1)*nw)=Ks((1:nw)+(tips(ic,it)-1)*nw);
                        K1=scalamp(found)*Ks(found+(tips(ic,it)-1)*nw);
                        K2=scalamp(found)*Ks(found+length(ind)+(tips(ic,it)-1)*nw,:);
                    end
                    [K1eq,Oc]=MaxHoopStressCriterion(K1,K2);
%                    Oc=0
                    if abs(Oc)>pi/6
                        remesh(ic,it)=1;
                    end
                    if min(lf)/lf(ic,it)>0.1
                        dasi(ic,it)=min(lf)/lf(ic,it)*dao;
                        
                        path=zones{2,cracks(ic)}*[1;1i];
                        t=-gradient(path);
                        t=t./abs(t);
                        switch it
                            case 1
                                t=t(1);
                                n=t*exp(1i*pi/2);
                                tip=path(1)+dasi(ic,it)*(cos(Oc)*t+sin(Oc)*n);
                                if GetSignedDistanceToZone(model,roi,real(tip),imag(tip))<-(2-1.*xfem)*zones{8,cracks(ic)}
                                    path=[tip;path];
                                else
                                    dasi(ic,it)=0;
                                end
                                
                            case 2
                                t=t(end);
                                n=t*exp(1i*pi/2);
                                tip=path(1)+dasi(ic,it)*(cos(Oc)*t+sin(Oc)*n);
                                if GetSignedDistanceToZone(model,roi,real(tip),imag(tip))<-(2-1.*xfem)*zones{8,cracks(ic)}
                                    path=[path;tip];
                                else
                                    dasi(ic,it)=0;
                                end
                        end
                        dacum(ic,it)=dacum(ic,it)+dasi(ic,it);
                        if dacum(ic,it)>min(5*dao,min(2*zones{8,cracks(ic)},10*zones{7,cracks(ic)}))
                            remesh(ic,it)=1;
                        end
                        zones{2,cracks(ic)}=[real(path),imag(path)];
                    end
                end
            end
        end
        Kss=[Kss,Ksi];
        das=[das,dasi];
        if any(abs(dasi)>0)
            model.zone=zones;
            if ~xfem
                model=Remeshing(param,model,~any(remesh==1));
                if any(remesh==1)
                    remesh=0*remesh;
                    dacum=0*dacum;
                end
            end
        else
            break;
        end
    else
        model=modeli;
        break;
    end
    step=step+1;
end
da=das;
Ks=Kss;
tips=tipso;
xo=xos;yo=yos;
conn=conns;
U=Us;
param=param0;
for iz=1:size(model.zone,2)
    zone=model.zone(:,iz);
    switch zone{4}
        case 6
            for izz=1:size(model0.zone,2)
                if strcmp(get(model0.zone{3,izz},'Tag'),get(model.zone{3,iz},'Tag'))
                    model.zone{7,iz}=model0.zone{7,izz};
                    break
                end
            end
    end
end

model.zone(:,~((cell2mat(model.zone(4,:))==5)|(cell2mat(model.zone(4,:))==6)))=[];
cracks=find(cell2mat(model.zone(4,:))==5);
model0.zone=model.zone;
model=model0;

load(fullfile('TMP','0_mesh_0'),'ng','rflag','rint');
save(param.result_file,'U','Ks','tips','cracks','xo','yo','param','model','nmod','conn','rint','ng','rflag','-v7.3');
end