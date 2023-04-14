function writeINP(filreso)
nmod=1;

tic
[pp,filres,ext]=fileparts(filreso);
load(filreso,'-mat')
if (nargin<2)
    if ~exist('selected')
        selected=ones(length(xo),1);
    end
else
    selected=selected2;
end

roi=param.roi;
dflag=exist('zo');
rint=true;
if isfield(model,'reduced_integration')
    rint=model.reduced_integration;
end
nlflag=false;
if isfield(model,'nlgeom')
    nlflag=model.nlgeom;
end
cpflag=false;
if isfield(model,'plane_stress')
    cpflag=model.plane_stress;
end
shell=false;
if isfield(model,'shell_element')
    shell=model.shell_element;
end

filres=[filres,''];
filres0=[filres,'.inp'];

%%
fwid = fopen(filres0,'w');
count = fprintf(fwid,'*HEADING\n');
count = fprintf(fwid,'*Preprint,echo=NO,history=NO,model=NO, contact=NO\n');
count = fprintf(fwid,'UFRECKLES - Laboratoire de Mecanique des Contacts et des Structures - INSA de Lyon\n');
count = fprintf(fwid,'UFRECKLES - GeM - Ecole Centrale de Nantes\n');


count = fprintf(fwid,'*Node, System=R,Nset=allnodes\n');
if rint
    ng=ones(length(elt),1);
else
    ng=4*(elt==4)+1*(elt==3)+2*(elt==6)+8*(elt==8);
end
ngc=[0;cumsum(ng(1:length(ng)-1))];
if dflag
    foundt3=find(elt==6);
    foundq4=find(elt==8);
    data=[(1:prod(Nnodes))',xo,yo,zo]';
    count = fprintf(fwid, '%d, %f, %f , %f \n', data);
    if any(foundt3)
        if shell
            count = fprintf(fwid,'*Element, type=SC6R,Elset=MAIN\n');
        else
            if rint
                count = fprintf(fwid,'*Element, type=C3D6R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=C3D6,Elset=MAIN\n');
            end
        end

        data=[(1:length(foundt3))',conn(foundt3,1:6)]';
        count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d\n',data);
    end
    if any(foundq4)
        if shell
            count = fprintf(fwid,'*Element, type=SC8R,Elset=MAIN\n');
        else
            if rint
                count = fprintf(fwid,'*Element, type=C3D8R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=C3D8,Elset=MAIN\n');
            end
        end

        data=[(1:length(foundq4))',conn(foundq4,1:8)]';
        count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',data);
    end
else
    data=[(1:prod(Nnodes))',xo,yo]';
    foundt3=find(elt==3);
    foundq4=find(elt==4);
    count = fprintf(fwid, '%d, %f, %f \n', data);
    if any(foundt3)
        if cpflag
                count = fprintf(fwid,'*Element, type=CPS3,Elset=MAIN\n');
        else
                count = fprintf(fwid,'*Element, type=CPE3,Elset=MAIN\n');
        end
        data=[(1:length(foundt3))',conn(foundt3,1:3)]';
        count = fprintf(fwid,'%d, %d, %d, %d\n',data);
    end
    if any(foundq4)
        if cpflag
            if rint
                count = fprintf(fwid,'*Element, type=CPS4R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=CPS4,Elset=MAIN\n');
            end
        else
            if rint
                count = fprintf(fwid,'*Element, type=CPE4R,Elset=MAIN\n');
            else
                count = fprintf(fwid,'*Element, type=CPE4,Elset=MAIN\n');
            end
        end
        data=[(1:length(foundq4))',conn(foundq4,1:4)]';
        count = fprintf(fwid,'%d, %d, %d, %d, %d\n',data);
    end
end
%%
fclose(fwid);
fprintf(1,'inp export done in %5.3f s\n',toc);

end
