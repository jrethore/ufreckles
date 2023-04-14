function writeAnsys(U,nmod,step,restart,action)
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
dflag=isfield(param,'extrusion_parameters')||(length(param0.roi)==6);
nlflag=false;
if isfield(param,'nlgeom')
    nlflag=param.nlgeom;
end
cpflag=false;
if isfield(param,'plane_stress')
    cpflag=param.plane_stress;
end

if nargin<5, action='init';end
if nargin<4 ,restart=0;end

if dflag
    load(fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,0)),'xo','yo','zo','elt','conn','Nnodes');
else
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,0)),'xo','yo','elt','conn','Nnodes');
end
switch action
    case 'init'
        jobname=sprintf('%s_k',param0.result_file);
        jobnameu=sprintf('%s',param0.result_file);
        unix(['rm dummy']);
        unix(['rm ',jobname,'.*']);
        unix(['rm ',jobnameu,'.*']);
        unix(['rm dump_nodal_forces.mac']);
        fwid = fopen('dump_nodal_forces.mac','w');
        count = fprintf(fwid,'/NOPR\n');
        count = fprintf(fwid,'/POST1\n');
        count = fprintf(fwid,'ALLSEL\n');
        count = fprintf(fwid,'/OUTPUT, dummy\n');
        count = fprintf(fwid,'*GET,nodeCount,NODE,0,COUNT\n');
        if dflag
            count = fprintf(fwid,'*DIM,arNodalForce,ARRAY,3*nodeCount,1\n');
        else
            count = fprintf(fwid,'*DIM,arNodalForce,ARRAY,2*nodeCount,1\n');
        end
        count = fprintf(fwid,'currentNode = 0\n');
        count = fprintf(fwid,'*DO, ar30,0,nodeCount - 1\n');
        count = fprintf(fwid,'	ESEL,ALL\n');
        count = fprintf(fwid,'	NSEL,ALL\n');
        count = fprintf(fwid,'	currentNode = NDNEXT(currentNode)\n');
        count = fprintf(fwid,'	NSEL,S,NODE,,currentNode\n');
        count = fprintf(fwid,'	FSUM\n');
        count = fprintf(fwid,'	*GET,arNodalForce( ar30 + 1+0*nodeCount, 1),FSUM,0,ITEM,FX\n');
        count = fprintf(fwid,'	*GET,arNodalForce( ar30 + 1+1*nodeCount, 1),FSUM,0,ITEM,FY\n');
        if dflag
            count = fprintf(fwid,'	*GET,arNodalForce( ar30 + 1+2*nodeCount, 1),FSUM,0,ITEM,FZ\n');
        end
        count = fprintf(fwid,'*ENDDO\n');
        count = fprintf(fwid,'*CFOPEN, nodalForces,txt,,\n');
        count = fprintf(fwid,'*VWRITE, arNodalForce(1,1)\n');
        count = fprintf(fwid,'(E)\n');
        count = fprintf(fwid,'*CFCLOS\n');
        count = fprintf(fwid,'/OUTPUT\n');
        count = fprintf(fwid,'/SYS, rm dummy\n');

        fclose(fwid);

    case {'Fint'}
        unix(['rm nodalForces.txt']);
        filres=sprintf('%d_%s_%d.inp',nmod,param0.result_file,step);
        jobname=sprintf('%s',param0.result_file);

        if exist(filres,'file')
            unix(['rm ',strrep(filres,'.inp','.*')]);
        end

        fwid = fopen(filres,'w');
        count = fprintf(fwid,'fini\n/clear\n');
        count = fprintf(fwid,'/NOPR\n');
        if ~restart
            count = fprintf(fwid,'/PNUM,NODE,1\n');
            count = fprintf(fwid,'/PREP7\n');


            if dflag
                foundt3=find(elt==6);
                foundq4=find(elt==8);
                data=[(1:prod(Nnodes))',xo,yo,zo]';
                count = fprintf(fwid, 'N,%d, %f, %f , %f \n', data);
                if any(foundt3)
                    error('not yet')
                end
                if any(foundq4)
                    %                        count = fprintf(fwid,'ET,1,SOLID45,,,,,1\n');
                    count = fprintf(fwid,'ET,1,SOLID185\n');
                    count = fprintf(fwid,'KEYOPT,1,1,0\n');

                    data=[conn(foundq4,1:8)]';
                    count = fprintf(fwid,'E, %d, %d, %d, %d, %d, %d, %d, %d\n',data);
                end


            else
                foundt3=find(elt==3);
                foundq4=find(elt==4);
                data=[(1:prod(Nnodes))',xo,yo]';
                count = fprintf(fwid, 'N,%d, %f, %f \n', data);
                if any(foundt3)
                    error('not yet')
                end
                if any(foundq4)
                    count = fprintf(fwid,'ET,1,PLANE182\n');
                    count = fprintf(fwid,'KEYOPT,1,1,0\n');
                    count = fprintf(fwid,'KEYOPT,1,6,0\n');
                    if cpflag
                        count = fprintf(fwid,'KEYOPT,1,3,0\n');
                    else
                        count = fprintf(fwid,'KEYOPT,1,3,2\n');
                    end
                    data=[conn(foundq4,1:4)]';
                    count = fprintf(fwid,'E, %d, %d, %d, %d\n',data);
                end
            end
            matmod=param.material_model;
            buffer=writeAnsysMat(matmod,param.material_parameters);
            count = fprintf(fwid,'\n');
            count = fprintf(fwid,'%s',buffer);
            count = fprintf(fwid,'\n');
        end
        if restart
            count = fprintf(fwid,'RESUME\n');
        end
        count = fprintf(fwid,'/SOLU\n');
        if restart
            count = fprintf(fwid,'ANTYPE,static,RESTART,,,CONTINUE\n');
        else
            count = fprintf(fwid,'ANTYPE,static,new\n');
        end
        count = fprintf(fwid,'SOLCON,ON\n');
        %        count = fprintf(fwid,'RESCONT,DEFINE,LAST,LAST\n');
        if nlflag
            count = fprintf(fwid,'NLGEOM,ON\n');
        else
            count = fprintf(fwid,'NLGEOM,OFF\n');
        end
        count = fprintf(fwid,'TIME,%d\n',step);
        count = fprintf(fwid,'NSUBST,1,100,1\n');
        data=[(1:prod(Nnodes))',U((1:prod(Nnodes)))]';
        count = fprintf(fwid, 'D, %d, UX, %12.5e\n', data);

        data=[(1:prod(Nnodes))',U(prod(Nnodes)+(1:prod(Nnodes)))]';
        count = fprintf(fwid, 'D, %d, UY, %12.5e\n', data);

        if dflag

            data=[(1:prod(Nnodes))',U(2*prod(Nnodes)+(1:prod(Nnodes)))]';
            count = fprintf(fwid, 'D, %d, UZ, %12.5e\n', data);

        end

        %                count = fprintf(fwid,'OUTPR,BASIC,ALL\n');
        %              count = fprintf(fwid,'/OUT ,%s,out\n',filres);
        count = fprintf(fwid,'/OUTPUT ,%s,out\n',strrep(filres,'.inp',''));

        count = fprintf(fwid,'EQSLV,sparse\n');


        count = fprintf(fwid,'SOLVE\n');
        count = fprintf(fwid,'SAVE\n');
        count = fprintf(fwid,'dump_nodal_forces\n');
        count = fprintf(fwid,'/OUTPUT\n');
        fclose(fwid);
    case 'K'

        filres=sprintf('%d_%s_k_%d.inp',nmod,param0.result_file,step);
        jobname=sprintf('%s_k',param0.result_file);
        jobnameu=sprintf('%s',param0.result_file);

        if exist(filres,'file')
            unix(['rm ',strrep(filres,'.inp','.*')]);
        end

        fwid = fopen(filres,'w');
        count = fprintf(fwid,'fini\n/clear\n');
        count = fprintf(fwid,'/NOPR\n');
        if ~restart
            count = fprintf(fwid,'/PNUM,NODE,1\n');
            count = fprintf(fwid,'/PREP7\n');


            if dflag
                foundt3=find(elt==6);
                foundq4=find(elt==8);
                data=[(1:prod(Nnodes))',xo,yo,zo]';
                count = fprintf(fwid, 'N,%d, %f, %f , %f \n', data);
                if any(foundt3)
                    error('not yet')
                end
                if any(foundq4)
                    %                        count = fprintf(fwid,'ET,1,SOLID45,,,,,1\n');
                    count = fprintf(fwid,'ET,1,SOLID185\n');
                    count = fprintf(fwid,'KEYOPT,1,1,0\n');

                    data=[conn(foundq4,1:8)]';
                    count = fprintf(fwid,'E, %d, %d, %d, %d, %d, %d, %d, %d\n',data);
                end


            else
                foundt3=find(elt==3);
                foundq4=find(elt==4);
                data=[(1:prod(Nnodes))',xo,yo]';
                count = fprintf(fwid, 'N,%d, %f, %f \n', data);
                if any(foundt3)
                    error('not yet')
                end
                if any(foundq4)
                    count = fprintf(fwid,'ET,1,PLANE182\n');
                    count = fprintf(fwid,'KEYOPT,1,1,0\n');
                    count = fprintf(fwid,'KEYOPT,1,6,0\n');
                    if cpflag
                        count = fprintf(fwid,'KEYOPT,1,3,0\n');
                    else
                        count = fprintf(fwid,'KEYOPT,1,3,2\n');
                    end
                    data=[conn(foundq4,1:4)]';
                    count = fprintf(fwid,'E, %d, %d, %d, %d\n',data);
                end
            end
            matmod=param.material_model;
            buffer=writeAnsysMat(matmod,param.material_parameters);
            count = fprintf(fwid,'\n');
            count = fprintf(fwid,'%s',buffer);
            count = fprintf(fwid,'\n');
        end
        if restart
            count = fprintf(fwid,sprintf('RESUME,%s,%s\n',jobnameu,'db'));
        end


cumiter=step+1;
        count = fprintf(fwid,'/SOLU\n');
        count = fprintf(fwid,'TIME,%d\n',step);
        count = fprintf(fwid,'DDELE,ALL,ALL\n');
        count = fprintf(fwid,'EQSLV,sparse\n');
        count = fprintf(fwid,'NCNV,2,,%d\n',cumiter);
        count = fprintf(fwid,'SOLVE\n');
        count = fprintf(fwid,'FINISH\n');
        count = fprintf(fwid,'/AUX2\n');
        count = fprintf(fwid,'FILE,%s,full\n',jobname);
        count = fprintf(fwid,'HBMAT,%s,,,ascii,stiff,no,yes\n',strrep(filres,'.inp',''));
        count = fprintf(fwid,'FINISH\n');
        count = fprintf(fwid,'*SMAT,matK,D,IMPORT,FULL,%s,STIFF\n',[jobname,'.full']);
        count = fprintf(fwid,'*EXPORT,matK,MATLAB,%s\n',strrep(filres,'.inp','.txt'));
        count = fprintf(fwid,'/OUTPUT\n');
        fclose(fwid);
    otherwise
        error('UNKNOWN ACTION')


end

end
