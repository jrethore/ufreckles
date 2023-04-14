function [xo,yo,zo,conn,elt,selected]=ReadVTK3D(meshfile)
%%
warn=0;
Nnodes=ones(1,3);
Nelems=ones(1,3);
fid=fopen(meshfile);

% toto=fgets(fid);
% toto=fgets(fid);
% toto=fgets(fid);
% toto=fgets(fid);
% Nnodes(1)=fscanf(fid,'POINTS %d double');
% toto=fscanf(fid,'%f',[prod(Nnodes)*3]);
% toto=reshape(toto,3,prod(Nnodes));
% xo=toto(1,:)';
% yo=toto(2,:)';
% zo=toto(3,:)';
% toto=fgets(fid);
% toto=fscanf(fid,'CELLS %d %d');
% Nelems(1)=toto(1);

parse = [];
while isempty(parse)
    parse=sscanf(strtrim(fgetl(fid)),'POINTS %d double');
end
Nnodes(1)=parse;

toto=fscanf(fid,'%f',[prod(Nnodes)*3]);
toto=reshape(toto,3,prod(Nnodes));
xo=toto(1,:)';
yo=toto(2,:)';
zo=toto(3,:)';

parse = [];
while isempty(parse)
    parse=sscanf(strtrim(fgetl(fid)),'CELLS %d %d');
end
Nelems(1)=parse(1);


selected=ones(prod(Nnodes),1);
conno=zeros(prod(Nelems),8);
conn=zeros(prod(Nelems),8);
elto=zeros(prod(Nelems),1);
elt=zeros(prod(Nelems),1);
[xgt,ygt,zgt,wgt]=GetGaussPointsTetrahedron(1);
Nt_r=[-1+0*xgt,1+0*xgt,0*ygt,0*zgt];
Nt_s=[-1+0*xgt,0*xgt,1+0*ygt,0*zgt];
Nt_t=[-1+0*xgt,0*xgt,0*ygt,1+0*zgt];

[xgp,ygp,zgp,wgp]=GetGaussPointsWedge(1);
Np_r=[-0.5*(1-zgp),0.5*(1-zgp),(0*ygp),...
    -0.5*(1+zgp),0.5*(1+zgp),(0*ygp)];
Np_s=[-0.5*(1-zgp),(0*xgp),0.5*(1-zgp),...
    -0.5*(1-zgp),(0*xgp),0.5*(1-zgp)];
Np_t=[-0.5*(1-xgp-ygp),-0.5*xgp,-0.5*ygp,...
    0.5*(1-xgp-ygp),0.5*xgp,0.5*ygp];
[xgq,ygq,zgq,wgq]=GetGaussPointsHexaedron(1);
Nq_r=[-0.125*(1-ygq).*(1-zgq),0.125*(1-ygq).*(1-zgq),0.125*(1+ygq).*(1-zgq),-0.125*(1+ygq).*(1-zgq),...
    -0.125*(1-ygq).*(1+zgq),0.125*(1-ygq).*(1+zgq),0.125*(1+ygq).*(1+zgq),-0.125*(1+ygq).*(1+zgq)];
Nq_s=[-0.125*(1-xgq).*(1-zgq),-0.125*(1+xgq).*(1-zgq),0.125*(1+xgq).*(1-zgq),0.125*(1-xgq).*(1-zgq),...
    -0.125*(1-xgq).*(1+zgq),-0.125*(1+xgq).*(1+zgq),0.125*(1+xgq).*(1+zgq),0.125*(1-xgq).*(1+zgq)];
Nq_t=[-0.125*(1-xgq).*(1-ygq),-0.125*(1+xgq).*(1-ygq),-0.125*(1+xgq).*(1+ygq),-0.125*(1-xgq).*(1+ygq),...
    0.125*(1-xgq).*(1-ygq),0.125*(1+xgq).*(1-ygq),0.125*(1+xgq).*(1+ygq),0.125*(1-xgq).*(1+ygq)];

for ii=1:prod(Nelems)
    toto=fscanf(fid,'%d',1);
    titi=fscanf(fid,'%d',toto);
    conno(ii,1:toto)=titi+1;
    elto(ii)=toto;
end

% toto=fgets(fid);
% toto=fgets(fid);
% toto=fscanf(fid,'CELL_TYPES %d');
parse = [];
while isempty(parse)
    parse=sscanf(strtrim(fgetl(fid)),'CELL_TYPES %d');
end
if Nelems(1)~=parse
    error('wrong VTK file')
end
for ii=1:prod(Nelems)
    toto=fscanf(fid,'%d',1);
    switch toto
        case {5,9,2}
            selected(conno(ii,1:elto(ii)))=0;
            elt(ii)=1;
        case 10
            titi=conno(ii,1:elto(ii))-1;
            dxdr=Nt_r*xo(titi+1);
            dydr=Nt_r*yo(titi+1);
            dzdr=Nt_r*zo(titi+1);
            dxds=Nt_s*xo(titi+1);
            dyds=Nt_s*yo(titi+1);
            dzds=Nt_s*zo(titi+1);
            dxdt=Nt_t*xo(titi+1);
            dydt=Nt_t*yo(titi+1);
            dzdt=Nt_t*zo(titi+1);
            J=[dxdr,dxds,dxdt;...
                dydr,dyds,dydt;...
                dzdr,dzds,dzdt];
            detJ=det(J);
            if detJ>0
                conn(ii,1:4)=titi+1;
            else
                conn(ii,1:4)=titi([1,4,3,2])+1;
                dxdr=Nt_r*xo(conn(ii,1:4));
                dydr=Nt_r*yo(conn(ii,1:4));
                dzdr=Nt_r*zo(conn(ii,1:4));
                dxds=Nt_s*xo(conn(ii,1:4));
                dyds=Nt_s*yo(conn(ii,1:4));
                dzds=Nt_s*zo(conn(ii,1:4));
                dxdt=Nt_t*xo(conn(ii,1:4));
                dydt=Nt_t*yo(conn(ii,1:4));
                dzdt=Nt_t*zo(conn(ii,1:4));
                J=[dxdr,dxds,dxdt;...
                    dydr,dyds,dydt;...
                    dzdr,dzds,dzdt];
                detJ=det(J);
                if detJ<0
                    conn(ii,1:4)=titi([1,3,2,4])+1;
                    dxdr=Nt_r*xo(conn(ii,1:4));
                    dydr=Nt_r*yo(conn(ii,1:4));
                    dzdr=Nt_r*zo(conn(ii,1:4));
                    dxds=Nt_s*xo(conn(ii,1:4));
                    dyds=Nt_s*yo(conn(ii,1:4));
                    dzds=Nt_s*zo(conn(ii,1:4));
                    dxdt=Nt_t*xo(conn(ii,1:4));
                    dydt=Nt_t*yo(conn(ii,1:4));
                    dzdt=Nt_t*zo(conn(ii,1:4));
                    J=[dxdr,dxds,dxdt;...
                        dydr,dyds,dydt;...
                        dzdr,dzds,dzdt];
                    detJ=det(J);
                    if detJ<0
                        conn(ii,1:4)=titi([1,2,4,3])+1;
                        dxdr=Nt_r*xo(conn(ii,1:4));
                        dydr=Nt_r*yo(conn(ii,1:4));
                        dzdr=Nt_r*zo(conn(ii,1:4));
                        dxds=Nt_s*xo(conn(ii,1:4));
                        dyds=Nt_s*yo(conn(ii,1:4));
                        dzds=Nt_s*zo(conn(ii,1:4));
                        dxdt=Nt_t*xo(conn(ii,1:4));
                        dydt=Nt_t*yo(conn(ii,1:4));
                        dzdt=Nt_t*zo(conn(ii,1:4));
                        J=[dxdr,dxds,dxdt;...
                            dydr,dyds,dydt;...
                            dzdr,dzds,dzdt];
                        detJ=det(J);
                        if detJ<0
                            warn=1;
                        end
                    end
                end
                
            end
            elt(ii)=4;
        case 13
            titi=conno(ii,1:elto(ii))-1;
            dxdr=Np_r*xo(titi+1);
            dydr=Np_r*yo(titi+1);
            dzdr=Np_r*zo(titi+1);
            dxds=Np_s*xo(titi+1);
            dyds=Np_s*yo(titi+1);
            dzds=Np_s*zo(titi+1);
            dxdt=Np_t*xo(titi+1);
            dydt=Np_t*yo(titi+1);
            dzdt=Np_t*zo(titi+1);
            J=[dxdr,dxds,dxdt;...
                dydr,dyds,dydt;...
                dzdr,dzds,dzdt];
            detJ=det(J);
            if detJ>0
                conn(ii,1:6)=titi+1;
            else
                warn=1;
                conn(ii,6:-1:1)=titi+1;
            end
            
            elt(ii)=6;
        case 12
            titi=conno(ii,1:elto(ii))-1;
            dxdr=Nq_r*xo(titi+1);
            dydr=Nq_r*yo(titi+1);
            dzdr=Nq_r*zo(titi+1);
            dxds=Nq_s*xo(titi+1);
            dyds=Nq_s*yo(titi+1);
            dzds=Nq_s*zo(titi+1);
            dxdt=Nq_t*xo(titi+1);
            dydt=Nq_t*yo(titi+1);
            dzdt=Nq_t*zo(titi+1);
            J=[dxdr,dxds,dxdt;...
                dydr,dyds,dydt;...
                dzdr,dzds,dzdt];
            detJ=det(J);
            if detJ>0
                conn(ii,1:8)=titi+1;
            else
                warn=1;
                conn(ii,8:-1:1)=titi+1;
            end
            
            elt(ii)=8;
            
        otherwise
            elt(ii)=1;
    end
    if warn
        keyboard
    end
end
found=find(elt>2);
Nelems(1)=length(found);
conn=conn(found,:);
elt=elt(found);
if warn
    display('CHECK FOR NEGATIVE DETJ');
end
fclose(fid);
%%
end