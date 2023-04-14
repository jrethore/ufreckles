function [xo,yo,zo,conn,elt,selected]=ReadVTK(meshfile)
%%
Nnodes=ones(1,3);
Nelems=ones(1,3);
fid=fopen(meshfile);

toto=fgets(fid);
toto=fgets(fid);
toto=fgets(fid);
toto=fgets(fid);
Nnodes(1)=fscanf(fid,'POINTS %d double');
toto=fscanf(fid,'%f',[prod(Nnodes)*3]);
toto=reshape(toto,3,prod(Nnodes));
xo=toto(1,:)';
yo=toto(2,:)';
zo=toto(3,:)';
toto=fgets(fid);
toto=fgets(fid);
toto=fscanf(fid,'CELLS %d %d');
Nelems(1)=toto(1);
selected=ones(prod(Nnodes),1);
conn=zeros(prod(Nelems),4);
elt=zeros(prod(Nelems),1);
[xgt,ygt,wgt]=GetGaussPointsTriangle(1,1);
Nt_r=[-1+0*xgt,1+0*xgt,0*ygt];
Nt_s=[-1+0*ygt,0*xgt,1+0*ygt];
[xgq,ygq,wgq]=GetGaussPointsQuadrangle(1,1);
Nq_r=[-0.25*(1-ygq),0.25*(1-ygq),0.25*(1+ygq),-0.25*(1+ygq)];
Nq_s=[-0.25*(1-xgq),-0.25*(1+xgq),0.25*(1+xgq),0.25*(1-xgq)];

for ii=1:prod(Nelems)
    toto=fscanf(fid,'%d',1);
    if toto==2
        toto=fscanf(fid,'%d',2);
        selected(toto+1)=0;
        elt(ii)=2;
    elseif toto==3
        toto=fscanf(fid,'%d',3);
        dxdr=Nt_r*xo(toto+1);
        dydr=Nt_r*yo(toto+1);
        dxds=Nt_s*xo(toto+1);
        dyds=Nt_s*yo(toto+1);
        detJ=(dxdr.*dyds-dydr.*dxds);
        if detJ>0
            conn(ii,1:3)=toto+1;
        else
            conn(ii,3:-1:1)=toto+1;
        end

        elt(ii)=3;
    elseif toto==4
        toto=fscanf(fid,'%d',4);

        dxdr=Nq_r*xo(toto+1);
        dydr=Nq_r*yo(toto+1);
        dxds=Nq_s*xo(toto+1);
        dyds=Nq_s*yo(toto+1);
        detJ=(dxdr.*dyds-dydr.*dxds);

        if detJ>0
            conn(ii,1:4)=toto+1;
        else
            conn(ii,4:-1:1)=toto+1;
        end

        elt(ii)=4;
    elseif toto==1
        toto=fscanf(fid,'%d',1);
        elt(ii)=1;

    else
        error('SCANNING ERROR');
    end
end
found=find(elt>2);
Nelems(1)=length(found);
conn=conn(found,:);
elt=elt(found);

fclose(fid);
%%
end