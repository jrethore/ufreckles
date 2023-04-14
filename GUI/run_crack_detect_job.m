function run_crack_detect_job(param0,model0)
nmod=0;
param=param0;
model=model0;

morph=1;

filreso=param.result_file;
filres='upropa.res';
param.result_file=filres;
model.mesh_file='upropa.vtk';
copyfile(model0.mesh_file,model.mesh_file);

pix2m=param.pixel_size;
modes=[1,2];ind=-3:7;
mu=0.5*model.material_parameters.young/(1+model.material_parameters.nu);
scal=2*mu*sqrt(2*pi);
scalamp=scal*(pix2m.^(1-(ind)*.5));
found=find(ind==1);
founds= find(ind==-1);
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
zoneso=zones;
fildef=param0.deformed_image;
if iscell(fildef)
    nim=size(fildef,2);
else
    nim=1;
end

filreso=strrep(filreso,'.res','');
filres=strrep(filres,'.res','');
delete([fullfile('VTK',filreso),'-0*.vtk']);

ztips=zeros(sum(tipso(:)>0),nim);
das=zeros(sum(tipso(:)>0),nim);
Kss=zeros(sum(tipso(:)>0)*nw,nim);
xos=cell(1,nim);
yos=cell(1,nim);
conns=cell(1,nim);
Us=cell(1,nim);

facx=0.5;

for iim=nim:-1:1
    if iscell(fildef)
        param.deformed_image=fildef{iim};
    else
        param.deformed_image=fildef;
    end
    iter=0;
    U=[];
    while iter<20
        modeli=model;
        Ksi=zeros(sum(tipso(:)>0)*nw,1);
        dxi=zeros(sum(tipso(:)>0),1);
        if morph
            run_fem_job(param,model,U)
        else
            run_fem_job(param,model,[])
        end
        load([filres,'.res'],'-mat','U','Ks','tips','xo','yo','conn')
        if iim==nim
            movefile(fullfile('VTK',[filres , sprintf('-%05d.vtk',0)]),fullfile('VTK',[filreso , sprintf('-%05d.vtk',0)]));
        end
        movefile(fullfile('VTK',[filres , sprintf('-%05d.vtk',1)]),fullfile('VTK',[filreso , sprintf('-%05d.vtk',iim)]));
        if any(tips(:)>0)
            xos{iim}=xo;
            yos{iim}=yo;
            conns{iim}=conn;
            Us{iim}=U;
            model=modeli;
            for ic=1:size(tipso,1)
                for it=1:2
                    if (tipso(ic,it)==1)&&(tips(ic,it)==1)
                        Kss((1:nw)+(tipso(ic,it)-1)*nw,iim)=Ks((1:nw)+(tips(ic,it)-1)*nw);
                        K1=scalamp(found)*Ks(found+(tips(ic,it)-1)*nw);
                        SK1=scalamp(founds)*Ks(founds+(tips(ic,it)-1)*nw);
                        dxi(tipso(ic,it))=-2*SK1/K1/pix2m
                    end
                end
            end
            das(:,iim)=das(:,iim)+facx*dxi
            for ic=1:size(tipso,1)
                path=zoneso{2,cracks(ic)}*[1;1i];
                s=[0;cumsum(abs(diff(path)))];
                for it=1:2
                    if (tipso(ic,it)==1)&&(tips(ic,it)==1)
                        switch it
                            case 1
                                stip=s(1)-das(tipso(ic,it),iim);
                                si=sort([stip;max(s,stip)]);
                            case 2
                                stip=s(end)+das(tipso(ic,it),iim);
                                si=sort([stip;min(stip,s)]);
                        end
                        [si]=unique(si);
                        path=interp1(s,path,si,'linear','extrap');
                        s=si;
                        ztips(tipso(ic,it),iim)=stip;
                        
                    end
                end
                zones{2,cracks(ic)}=[real(path),imag(path)];
            end
            model.zone=zones;
            model=Remeshing(param,model,morph*facx*dxi*(abs(dxi)>1&&iter<20));
            zones=model.zone;
            if ~any(abs(dxi)>1)
                das(:,max(1,iim-1))=das(:,iim);
                break;
            end
        else
            model=modeli;
            break;
        end
        iter=iter+1;
    end
    da=das;
    for ic=1:size(da,1)
        da(ic,:)=da(ic,:)-da(ic,1);
    end
    modelp=model;
    paramp=param;
    cracksp=cracks;
    Ks=Kss;
    tips=tipso;
    xo=xos;yo=yos;
    conn=conns;
    U=Us;
    param=param0;
    model=model0;
    model.zone(:,~((cell2mat(model.zone(4,:))==5)))=[];
    cracks=find(cell2mat(model.zone(4,:))==5);
    
    load(fullfile('TMP','0_mesh_0'),'ng','rflag','rint');
    save(param.result_file,'da','U','Ks','ztips','tips','cracks','xo','yo','param','model','nmod','conn','rint','ng','rflag','-v7.3');
    model=modelp;
    param=paramp;
    cracks=cracksp;
end


da=das;
for ic=1:size(da,1)
    da(ic,:)=da(ic,:)-da(ic,1);
end
Ks=Kss;
tips=tipso;
xo=xos;yo=yos;
conn=conns;
U=Us;
param=param0;


model=model0;
model.zone(:,~((cell2mat(model.zone(4,:))==5)))=[];
cracks=find(cell2mat(model.zone(4,:))==5);

load(fullfile('TMP','0_mesh_0'),'ng','rflag','rint');
save(param.result_file,'da','U','Ks','ztips','tips','cracks','xo','yo','param','model','nmod','conn','rint','ng','rflag','-v7.3');
end