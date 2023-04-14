function [xo,yo,zo,conn,elt,eset,selected]=ReadVINP(meshfile)
%%
fid=fopen(meshfile);
toto=1;
Nnodes=ones(1,3);
Nelems=ones(1,3);
selected=[];
eset={};
while isempty(toto==-1)|~(toto==-1)
toto=fgets(fid);
if strcmp(toto(1),'*')
    cas=toto(2:5);
    switch cas
    case 'Node'
 titi=fscanf(fid,'%f,%f,%f,%f');
           Nnodes(1)=titi(end-3);
           titi=reshape(titi,4,prod(Nnodes))';
           xo=titi(:,2);
           yo=titi(:,3);
           zo=titi(:,4);
           clear titi
           selected=ones(Nnodes);
    case 'Elem'
 titi=fscanf(fid,'%f,%f,%f,%f');
           Nelems(1)=titi(end-3);
           titi=reshape(titi,4,prod(Nelems))';
        conn=titi(:,2:4);
        clear titi
    case 'Else'
         titi=fscanf(fid,'%d,');
         eset{numel(eset)+1}=titi;
    case 'Nset'
         titi=fscanf(fid,'%d,');
         selected(titi)=0;
    end
end
end
elt=repmat(3,prod(Nelems),1);
fclose(fid);
end