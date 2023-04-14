function Fint=ComputeInternalForceVector(nmod,U,step)
persistent nmodo maskg epsxx epsyy epsxy epszz epsxz epsyz wdetJ doinit Fold Fprev stepold
iscale=1;
load(fullfile('TMP','params'),'param');
param0=param;
pix2m=param0.pixel_size;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
if isfield(param0,'helper')
    soft=param0.helper;
else
    soft='none';
end
dflag=isfield(param,'extrusion_parameters')||(length(param0.roi)==6);
nlflag=false;
if isfield(param,'nlgeom')
    nlflag=param.nlgeom;
end
reload=1;
if ~isempty(nmodo)
    if (nmodo==nmod)
        reload=0;
    else
        reload=1;
    end
end

reload=reload||nlflag;
switch soft
    case 'none'

        if ~isfield(param,'enrichment')
            if reload


                if dflag


                    load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,iscale-1)),'epsyy');
                    load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,iscale-1)),'epsxy');
                    load(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,iscale-1)),'epszz');
                    load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,iscale-1)),'epsxz');
                    load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,iscale-1)),'epsyz');
maskg=1;
                else
                    load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,iscale-1)),'epsyy');
                    load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,iscale-1)),'epsxy');
                    epsxz=sparse(size(epsxx,1),size(epsxx,2));epsyz=epsxz;epszz=epsxz;
                load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'maskg');
                end
                nmodo=nmod;
            end
            [sxx,syy,sxy,sxz,syz,szz]=ComputeStress(epsxx*U,epsyy*U,epsxy*U,nmod,epsxz*U,epsyz*U,epszz*U);
            Fint=(epsxx'*maskg*wdetJ*sxx+epsyy'*maskg*wdetJ*syy+epsxy'*maskg*wdetJ*sxy...
                +epsxz'*maskg*wdetJ*sxz+epsyz'*maskg*wdetJ*syz+epszz'*maskg*wdetJ*szz)*pix2m;
            if nlflag
                load(fullfile('TMP',sprintf('%d_matmod',nmod)),'DU');
                if dflag
                    meshfile=fullfile('TMP',sprintf('%d_3d_mesh_%d',nmod,0));
                else
                    meshfile=fullfile('TMP',sprintf('%d_mesh_%d',nmod,0));
                end
                load(meshfile,'xo','yo','zo','Nnodes');
nn=prod(Nnodes);
                Ux=0.5*DU((1:nn));
                Uy=0.5*DU(nn+(1:nn));
                if dflag
                    Uz=0.5*DU(2*nn+(1:nn));
                else
                    Uz=zeros(size(zi));
                end
                xo=xo+Ux;
                yo=yo+Uy;
                zo=zo+Uz;
                save(meshfile,'xo','yo','zo','-append');
                if dflag
                    CreateGradBasisFunction3D(1,nmod);
                else
                    CreateGradBasisFunction(1,nmod);
                end

            end
        else
            load(fullfile('TMP',sprintf('%d_k_operator_%d',nmod,iscale-1)),'K');
            Fint=K*U;
            if isfield(param0,'cz_length')&&isfield(param,'interface_model')
                load(fullfile('TMP',sprintf('%d_dphin_%d',nmod,iscale-1)),'dphincz');
                load(fullfile('TMP',sprintf('%d_dphit_%d',nmod,iscale-1)),'dphitcz');
                load(sprintf('TMP/%d_iphi_%d',nmod,iscale-1),'dz');
                [tn,tt]=ComputeCohesiveStress(dphincz*U,dphitcz*U,nmod);
                dz=diag(sparse(dz));
                Fint=Fint+(dphincz'*(dz*tn)+dphitcz'*(dz*tt))*pix2m;
            end

        end
    case 'abaqus'
        abqcmd='abaqus';
        if isfield(param0,'abaqus_command_line')
            abqcmd=param0.abaqus_command_line;
        end
        if nargin<3
            step=[];
        end

        if isempty(step)
            restart=0;
            filres=[sprintf('%d_',nmod),param0.result_file];
        else
            filres=sprintf('%d_%s_%d',nmod,param0.result_file,step);
            if step==1
                restart=0;
            else
                restart=1;
            end
        end
        if isempty(doinit)
            writeAbaqus(U,nmod,step,restart,'init');
            doinit=false;
        end

        writeAbaqus(U,nmod,step,restart,'Fint');

        if restart
            filresold=sprintf('%d_%s_%d',nmod,param0.result_file,step-1);
            unix([abqcmd,' job=',filres,' oldjob=',filresold,' double interactive']);
        else
            unix([abqcmd,' job=',filres,' double  interactive']);
        end
        Fint=readAbaqus(nmod,step,'Fint')*pix2m;
        %         writeAbaqus(U,nmod,step,restart,'S');
        %
        %         if restart
        %             filresold=sprintf('%d_%s_%d',nmod,param0.result_file,step-1);
        %             unix(['abaqus68 job=',filres,' oldjob=',filresold,' double interactive']);
        %         else
        %             unix(['abaqus68 job=',filres,' double  interactive']);
        %         end
        %         S=readAbaqus(nmod,step,'S');
        %         if size(S,2)==3
        %             if reload
        %                 load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,iscale-1)),'epsxx','wdetJ');
        %                 load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,iscale-1)),'epsyy');
        %                 load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,iscale-1)),'epsxy');
        %                 load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'maskg');
        %                 nmodo=nmod;
        %             end
        %         Fint=(epsxx'*maskg*wdetJ*S(:,1)+epsyy'*maskg*wdetJ*S(:,2)+epsxy'*maskg*wdetJ*S(:,3))*pix2m^2;
        %         elseif size(S,2)==6
        %             if reload
        %                 load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,iscale-1))),'epsxx','wdetJ');
        %                 load(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,iscale-1)),'epsyy');
        %                 load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,iscale-1)),'epsxy');
        %                 load(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,iscale-1)),'epszz');
        %                 load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,iscale-1)),'epsxz');
        %                 load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,iscale-1)),'epsyz');
        %                 load(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'maskg');
        %                 nmodo=nmod;
        %             end
        % %             U
        % %             [epsxx*U,epsyy*U,epszz*U,epsxy*U,epsxz*U,epsyz*U]
        %         Fint=(epsxx'*maskg*wdetJ*S(:,1)+epsyy'*maskg*wdetJ*S(:,2)+epszz'*maskg*wdetJ*S(:,3)...
        %             +epsxy'*maskg*wdetJ*S(:,4)+epsxz'*maskg*wdetJ*S(:,5)+epsyz'*maskg*wdetJ*S(:,6))*pix2m^2;
        %
        %         else
        %             error
        %         end

        %         Fprev=Fint;
        %         if ~isempty(step)
        %             if (step==1)
        %                 Fold=0;
        %                 stepold=0;
        %                 'step 1'
        %             else
        %                 'step>1'
        %             display(sprintf('step %d stepold %d',step,stepold))
        %                 if (step>stepold+1)
        %                     'step>stepold+1'
        %                     Fold=Fprev+Fold;
        %                     stepold=step-1;
        %
        %                 end
        %                 Fint=Fint+Fold;
        %             end
        %         end
    case 'ansys'
        if nargin<3
            step=1;
        end
        jobname=sprintf('%s',param0.result_file);
        filres=sprintf('%d_%s_%d',nmod,param0.result_file,step);
        if step==1
            restart=0;
            writeAnsys([],nmod,step,restart,'init');
        else
            restart=1;
        end
        %         if isempty(doinit)
        %             writeAnsys(U,nmod,step,restart,'init');
        %             doinit=false;
        %         end
        writeAnsys(U,nmod,step,restart,'Fint');


        unix(['/opt/ansys_inc/v121/ansys/bin/ansys121  -p AA_T_I  -j ',jobname,' -s read -l en-us -b < ',filres,'.inp | tee -i ',filres,'.out' ]);
        %            ['/opt/ansys_inc/v121/ansys/bin/ansys121  -p AA_T_I  -j ',jobname,' -s read -l en-us -b < ',filres,'.inp | tee -i ',filres,'.out' ]
        Fint=readAnsys(nmod,step,'Fint')*pix2m;

    otherwise
        error('Unknown helper');
end
end