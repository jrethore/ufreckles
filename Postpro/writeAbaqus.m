function writeAbaqus(U,nmod,step,restart,action)
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
dflag=isfield(param,'extrusion_parameters')||(length(param0.roi)==6);
rint=true;
if isfield(param,'reduced_integration')
    rint=param.reduced_integration;
end
nlflag=false;
if isfield(param,'nlgeom')
    nlflag=param.nlgeom;
end
cpflag=false;
if isfield(param,'plane_stress')
    cpflag=param.plane_stress;
end
shell=false;
if isfield(param,'shell_element')
    shell=param.shell_element;
end
filres0=sprintf('%d_%s0.inp',nmod,param0.result_file);
if nargin<5, action='init';end
if nargin<4 ,restart=0;end
if nargin<3
    filres=sprintf('%d_%s.inp',nmod,param0.result_file);
else
    if isempty(step)
        filres=sprintf('%d_%s.inp',nmod,param0.result_file);
    else
        filres=sprintf('%d_%s_%d.inp',nmod,param0.result_file,step);
    end
end
if exist(filres0,'file')&&strcmp(action,'init')
    unix(['rm ',filres0]);
end
if exist(filres,'file')
    unix(['rm ',strrep(filres,'.inp','.*')]);
end
if dflag
    load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,0)),'xo','yo','zo','elt','conn','Nnodes');
else
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,0)),'xo','yo','elt','conn','Nnodes');
end
switch action
    case 'init'
        fwid = fopen(filres0,'w');
        count = fprintf(fwid,'*HEADING\n');
        count = fprintf(fwid,'*Preprint,echo=NO,history=NO,model=NO, contact=NO\n');
        count = fprintf(fwid,'MIC - Laboratoire de Mecanique des Contacts et des Structures - INSA de Lyon\n');


        count = fprintf(fwid,'*Node, System=R,Nset=allnodes\n');


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
            foundt3=find(elt==3);
            foundq4=find(elt==4);
            data=[(1:prod(Nnodes))',xo,yo]';
            count = fprintf(fwid, '%d, %f, %f \n', data);
            if any(foundt3)
                if cpflag
                    if rint
                        count = fprintf(fwid,'*Element, type=CPS3R,Elset=MAIN\n');
                    else
                        count = fprintf(fwid,'*Element, type=CPS3,Elset=MAIN\n');
                    end
                else
                    if rint
                        count = fprintf(fwid,'*Element, type=CPE3R,Elset=MAIN\n');
                    else
                        count = fprintf(fwid,'*Element, type=CPE3,Elset=MAIN\n');
                    end
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
        fclose(fwid);
    case {'Fint','K','S'}
        if restart
            fwid = fopen(filres,'w');
            count = fprintf(fwid,'*HEADING\n');
            count = fprintf(fwid,'*RESTART,READ\n');
        else
            unix(['cp ',filres0,' ',filres]);
            fwid = fopen(filres,'a');
            load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
            matmod=param.material_model;
            buffer=writeAbaqusMat(matmod,param.material_parameters);
            if shell
                count = fprintf(fwid,'*SHELL SECTION,ELSET=MAIN,NODAL THICKNESS,MATERIAL=%s\n',matmod);
            else
                count = fprintf(fwid,'*SOLID SECTION,ELSET=MAIN,MATERIAL=%s\n',matmod);
            end
            count = fprintf(fwid,'*MATERIAL,NAME=%s\n',matmod);
            count = fprintf(fwid,'%s',buffer);
        end
        switch action
            case {'Fint','S'}
                if nlflag
                    count = fprintf(fwid,'*STEP,NLGEOM=YES\n');
                else
                    count = fprintf(fwid,'*STEP,NLGEOM=NO\n');
                end
                count = fprintf(fwid,'*STATIC\n');
                count = fprintf(fwid,'1,1,0.0000001,1\n');
                count = fprintf(fwid,'*BOUNDARY, OP=NEW\n');

                data=[(1:prod(Nnodes))',repmat(1,prod(Nnodes),2),U((1:prod(Nnodes)))]';
                count = fprintf(fwid, '%d, %d, %d, %12.5e\n', data);

                data=[(1:prod(Nnodes))',repmat(2,prod(Nnodes),2),U(prod(Nnodes)+(1:prod(Nnodes)))]';
                count = fprintf(fwid, '%d, %d, %d, %12.5e\n', data);

                if dflag

                    data=[(1:prod(Nnodes))',repmat(3,prod(Nnodes),2),U(2*prod(Nnodes)+(1:prod(Nnodes)))]';
                    count = fprintf(fwid, '%d, %d, %d, %12.5e\n', data);

                end

                switch action
                    case 'Fint'
                        count = fprintf(fwid,'*NODE PRINT, FREQUENCY=1000000000,NSET=allnodes\nRF\n');
%                        count = fprintf(fwid,'*NODE PRINT, FREQUENCY=1000000000,NSET=allnodes\nU\n');
                    case 'S'
                        error('DO NOT DO THAT');
                        count = fprintf(fwid,'*EL PRINT, FREQUENCY=1000000000,ELSET=MAIN\nS\n');
%                        count = fprintf(fwid,'*EL PRINT, FREQUENCY=1000000000,ELSET=MAIN\nE\n');
                end
                count = fprintf(fwid,'*RESTART,WRITE OVERLAY\n');


                count = fprintf(fwid,'*END STEP\n');

            case 'K'
                count = fprintf(fwid,'*STEP\n');
                count = fprintf(fwid,'*MATRIX GENERATE, STIFFNESS\n');
                count = fprintf(fwid,'*BOUNDARY, OP=NEW\n');
                count = fprintf(fwid,'*END STEP\n');

        end
        fclose(fwid);
    otherwise
        error('UNKNOWN ACTION')

end

end
