function genINPFile(nmod,Fos,Cs,Ups,Po)
if nargin<5,Po=0;end
load(fullfile('TMP','params'),'param');
param0=param;
pix2m=param0.pixel_size;
if isfield(param0,'time_step')
    dt=param0.time_step;
else
    dt=1;
end

load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
tscale=1;
if isfield(param,'time_scale')
    tscale=param.time_scale;
end


        C=Cs;Fo=Fos;Up=Ups;
        nstep=size(Up,2);
        if 0
            tn=1:size(Up,2);
        else
            
            if size(Up,2)>1
                tn=0:size(Up,2);
                ti=dt:dt:size(Up,2);
                Up=[zeros(size(Up,1),1),Up];
                if size(Fo,2)>1
                    Fo=[zeros(size(Fo,1),1),Fo];
                end
                if size(Po,2)>1
                    Po=[0,Po];
                end
            else
                tn=0:1;
                ti=dt:dt:1;
                Up=[zeros(size(Up)),Up];
            end
        end
        dflag=isfield(param,'extrusion_parameters')||(length(param0.roi)==6);
        dyn=false;
        if isfield(param,'dynamic')
            dyn=param.dynamic;
        end
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
        membrane=false;
        if isfield(param,'membrane_element')
            membrane=param.membrane_element;
        end
        if membrane,dflag=true;end
        filres=strrep(param0.result_file,'.res','.inp');
        if dflag
            load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,0)),'xo','yo','zo','elt','conn','Nnodes','selected');
            Nddl=3*prod(Nnodes);
        else
            load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,0)),'xo','yo','elt','conn','Nnodes','selected');
            Nddl=2*prod(Nnodes);
        end
        
        fwid = fopen(filres,'w');
        count = fprintf(fwid,'*HEADING\n');
        count = fprintf(fwid,'*Preprint,echo=NO,history=NO,model=NO, contact=NO\n');
        count = fprintf(fwid,'UFRECKLES - Laboratoire de Mecanique des Contacts et des Solides - INSA de Lyon\n');
        
        if isfield(param0,'abaqus_controls')
            count = fprintf(fwid,'*SECTION CONTROLS, NAME=ctrl,%s\n',param0.abaqus_controls);
        end
        
        count = fprintf(fwid,'*Node, System=R,Nset=allnodes\n');
        
        if dflag
            foundt3=find((elt==6)|(elt==3));
            foundq4=find((elt==8)|(elt==4));
            data=[(1:prod(Nnodes))',xo,yo,zo]';
            count = fprintf(fwid, '%d, %f, %f , %f \n', data);
            if any(foundt3)
                if shell
                    count = fprintf(fwid,'*Element, type=SC6R,Elset=MAIN\n');
                    data=[(1:length(foundt3))',conn(foundt3,1:6)]';
                    count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d\n',data);
                elseif membrane
                    count = fprintf(fwid,'*Element, type=M3D3,Elset=MAIN\n');
                    data=[(1:length(foundt3))',conn(foundt3,1:3)]';
                    count = fprintf(fwid,'%d, %d, %d, %d\n',data);
                    
                else
                    if rint
                        count = fprintf(fwid,'*Element, type=C3D6R,Elset=MAIN\n');
                    else
                        count = fprintf(fwid,'*Element, type=C3D6,Elset=MAIN\n');
                    end
                    data=[(1:length(foundt3))',conn(foundt3,1:6)]';
                    count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d\n',data);
                end
            end
            if any(foundq4)
                if shell
                    count = fprintf(fwid,'*Element, type=SC8R,Elset=MAIN\n');
                    data=[length(foundt3)+(1:length(foundq4))',conn(foundq4,1:8)]';
                    count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',data);
                elseif membrane
                    count = fprintf(fwid,'*Element, type=M3D4,Elset=MAIN\n');
                    data=[length(foundt3)+(1:length(foundq4))',conn(foundq4,1:4)]';
                    count = fprintf(fwid,'%d, %d, %d, %d, %d\n',data);
                else
                    if rint
                        count = fprintf(fwid,'*Element, type=C3D8R,Elset=MAIN\n');
                    else
                        count = fprintf(fwid,'*Element, type=C3D8,Elset=MAIN\n');
                    end
                    data=[length(foundt3)+(1:length(foundq4))',conn(foundq4,1:8)]';
                    count = fprintf(fwid,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',data);
                end
            end
            
            
        else
            foundt3=find(elt==3);
            foundq4=find(elt==4);
            data=[(1:prod(Nnodes))',xo,yo]';
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
                data=[length(foundt3)+(1:length(foundq4))',conn(foundq4,1:4)]';
                count = fprintf(fwid,'%d, %d, %d, %d, %d\n',data);
            end
        end
if ~isfield(param,'material_model')
           param.material_model='elastic_homogeneous_isotropic';         
    param.material_parameters.nu=0.3;
    param.material_parameters.young=200e9;
end
        matmod=param.material_model;
        buffer=writeAbaqusMat(matmod,param.material_parameters);
        if shell
            count = fprintf(fwid,'*SHELL SECTION,ELSET=MAIN,NODAL THICKNESS,MATERIAL=%s\n',matmod);
        elseif membrane
            count = fprintf(fwid,'*SURFACE, TYPE=ELEMENT, name=TOUT\nMAIN,SPOS\n');
            count = fprintf(fwid,'*MEMBRANE SECTION,ELSET=MAIN,MATERIAL=%s\n%d,\n',matmod,param.membrane_thickness);
        else
            count = fprintf(fwid,'*SOLID SECTION,ELSET=MAIN,MATERIAL=%s',matmod);
            if isfield(param0,'abaqus_controls')
                count = fprintf(fwid,',CONTROLS=ctrl');
            end
            count = fprintf(fwid,'\n');
        end
        count = fprintf(fwid,'*MATERIAL,NAME=%s\n',matmod);
        count = fprintf(fwid,'%s',buffer);
        
        [iddl,icond,val]=find(C);
        if 0
            comp=ceil(iddl/prod(Nnodes));
            inod=iddl-(comp-1)*prod(Nnodes);
            [inods,index]=sort(inod);
            new=[1;diff(inods)];
            inew=find(new);
            for in=1:length(inew)
                count = fprintf(fwid,'*NSET,NSET=N%d\n%d\n',inods(inew(in)),inods(inew(in)));
            end
            
            for in=1:length(iddl)
                comp=ceil(iddl(in)/prod(Nnodes));
                inod=iddl(in)-(comp-1)*prod(Nnodes);
                count = fprintf(fwid,'*AMPLITUDE, NAME=AMP_N%d_%d, DEFINITION=TABULAR\n',inod,comp);
                count = fprintf(fwid,'%d, %d\n',[tn*tscale;Up(in,:)]);
                
            end
            
            if any(abs(Po))
                count = fprintf(fwid,'*AMPLITUDE, NAME=AMP_PRESSION, DEFINITION=TABULAR\n');
                count = fprintf(fwid,'%d, %d\n',[tn*tscale;Po]);
                
            end
            
            if dyn
                count = fprintf(fwid,'*STEP, NAME=ALLSTEPS\n');
                count = fprintf(fwid,'*DYNAMIC, EXPLICIT\n, %12.5e\n',max(tn)*tscale);
                count = fprintf(fwid,'*BULK VISCOSITY\n0.06, 1.2\n');
            else
                if nlflag
                    count = fprintf(fwid,'*STEP,NLGEOM=YES, NAME=ALLSTEPS\n');
                else
                    count = fprintf(fwid,'*STEP,NLGEOM=NO, NAME=ALLSTEPS\n');
                end
                count = fprintf(fwid,'*STATIC\n');
                count = fprintf(fwid,'%12.5e,%12.5e,%12.5e,%12.5e\n',dt*tscale,max(tn)*tscale,dt*1.e-6*tscale,dt*tscale);
                if isfield(param0,'abaqus_controls')
                    count = fprintf(fwid,'%s\n',param0.abaqus_controls);
                end
            end
            
            
            for in=1:length(iddl)
                comp=ceil(iddl(in)/prod(Nnodes));
                inod=iddl(in)-(comp-1)*prod(Nnodes);
                count = fprintf(fwid,'*BOUNDARY, AMPLITUDE=AMP_N%d_%d\n',inod,comp);
                count = fprintf(fwid,'N%d, %d, %d, 1.\n',[inod,comp,comp]);
                
            end
            
            
            if any(abs(Po))
                count = fprintf(fwid,'*DSLOAD, AMPLITUDE=AMP_PRESSION\n TOUT,P, 1.\n');
            end
            
            
            count = fprintf(fwid,'*OUTPUT, FIELD, NUMBER INTERVAL=%d, TIME MARKS=YES\n',length(tn));
            count = fprintf(fwid,'*NODE OUTPUT\nU\nRF\n');
            count = fprintf(fwid,'*ELEMENT OUTPUT\nS\nLE\nPEEQ\n');
            count = fprintf(fwid,'*END STEP\n');
            
        else
            for n=1:length(ti)
                Ui=(interp1(tn,Up',ti(n)))';
                if size(Fo,2)>1
                    Foi=(interp1(tn,Fo',ti(n)))';
                else
                    Foi=Fo;
                end
                if size(Po,2)>1
                    Poi=(interp1(tn,Po,ti(n)));
                else
                    Poi=Po;
                end
                if nlflag
                    count = fprintf(fwid,'*STEP,NLGEOM=YES,AMPLITUDE=RAMP\n');
                else
                    count = fprintf(fwid,'*STEP,NLGEOM=NO,AMPLITUDE=RAMP\n');
                end
                count = fprintf(fwid,'*STATIC\n');
                count = fprintf(fwid,'0.1,1,0.0000001,1\n');
                %                 count = fprintf(fwid,'*CONTROLS, PARAMETERS=LINE SEARCH\n');
                %                 count = fprintf(fwid,'10,1,0.0001,0.25,0.01\n');
                count = fprintf(fwid,'*BOUNDARY\n');
                for in=1:length(icond)
                    comp=ceil(iddl(in)/prod(Nnodes));
                    inod=iddl(in)-(comp-1)*prod(Nnodes);
                    count = fprintf(fwid, '%d, %d, %d, %12.5e\n', inod,comp,comp,Ui(in));
                end
                
                
                if any(abs(Foi)>0)
                    iconf=find(abs(Foi)>0);
                    count = fprintf(fwid,'*CLOAD, OP=NEW\n');
                    for in=1:length(iconf)
                        comp=ceil(iconf(in)/prod(Nnodes));
                        inod=iconf(in)-(comp-1)*prod(Nnodes);
                        count = fprintf(fwid, '%d, %d, %12.5e\n', inod,comp,Foi(iconf(in)));
                    end
                end
                if abs(Poi)>0
                    count = fprintf(fwid,'*DSLOAD\n TOUT,P, %d\n',Poi');
                end
                
                
                count = fprintf(fwid,'*END STEP\n');
                
            end
        end
        
        %%
        fclose(fwid);



end


