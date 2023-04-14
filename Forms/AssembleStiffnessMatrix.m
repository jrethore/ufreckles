function Kp=AssembleStiffnessMatrix(iscale,nmod,step)
persistent nmodo iscaleo maskg epsxx epsyy epsxy epsxz epsyz epszz wdetJ doinit
persistent Kef
%persistent phig phip
eflag=0;
load(fullfile('TMP','params'));
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

switch soft
    case {'none','ansys'}
        %    case {'none'}
        if strcmp(soft,'ansys')
            eflag=1;
        end
        restart=1;
        if ~isempty(nmodo)
            if (nmodo==nmod)
                restart=0;
            else
                restart=1;
            end
        end
        if ~isempty(iscaleo)
            if (iscaleo==iscale)
                restart=0;
            else
                restart=1;
            end
        end
        restart=restart||nlflag;
        if isfield(param,'enrichment')&&(iscale==1)
            if ~isfield(param,'interface_model')
                restart=1;
            end
            if restart
                if dflag
                    load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
                    load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,10*(iscale-1))),'epsxy');
                    load(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,10*(iscale-1))),'epszz');
                    load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,10*(iscale-1))),'epsxz');
                    load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,10*(iscale-1))),'epsyz');
                    maskg=1;

                else
                    load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
                    load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'epsxy');
                    epsxz=sparse(size(epsxx,1),size(epsxx,2));epsyz=epsxz;epszz=epsxz;
                    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,(iscale-1))),'maskg');
                end
%                maskg=1
                nmodo=nmod;
                iscaleo=iscale;
                [sxx,syy,sxy,sxz,syz,szz]=ComputeTangentStress(epsxx,epsyy,epsxy,nmod,epsxz,epsyz,epszz,eflag);
                Kp= epsxx'*maskg*wdetJ*sxx...
                    +epsyy'*maskg*wdetJ*syy...
                    +epszz'*maskg*wdetJ*szz...
                    +epsxy'*maskg*wdetJ*sxy...
                    +epsxz'*maskg*wdetJ*sxz...
                    +epsyz'*maskg*wdetJ*syz;
                if dflag
                    load(fullfile('TMP',sprintf('%d_3d_xepsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_3d_xepsyy_%d',nmod,10*(iscale-1))),'epsyy');
                    load(fullfile('TMP',sprintf('%d_3d_xepsxy_%d',nmod,10*(iscale-1))),'epsxy');
                    load(fullfile('TMP',sprintf('%d_3d_xepszz_%d',nmod,10*(iscale-1))),'epszz');
                    load(fullfile('TMP',sprintf('%d_3d_xepsxz_%d',nmod,10*(iscale-1))),'epsxz');
                    load(fullfile('TMP',sprintf('%d_3d_xepsyz_%d',nmod,10*(iscale-1))),'epsyz');
                    masks=1;
                else
                load(fullfile('TMP',sprintf('%d_xepsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                load(fullfile('TMP',sprintf('%d_xepsyy_%d',nmod,10*(iscale-1))),'epsyy');
                load(fullfile('TMP',sprintf('%d_xepsxy_%d',nmod,10*(iscale-1))),'epsxy');
                    epsxz=sparse(size(epsxx,1),size(epsxx,2));epsyz=epsxz;epszz=epsxz;
                load(fullfile('TMP',sprintf('%d_mask_%d',nmod,(iscale-1))),'masks');
                end
                Nddl_tot=size(epsxx,2);
                [sxx,syy,sxy,sxz,syz,szz]=ComputeTangentStress(epsxx,epsyy,epsxy,nmod,epsxz,epsyz,epszz,eflag);
                Kx= epsxx'*masks*wdetJ*sxx...
                    +epsyy'*masks*wdetJ*syy...
                    +epszz'*masks*wdetJ*szz...
                    +epsxy'*masks*wdetJ*sxy...
                    +epsxz'*masks*wdetJ*sxz...
                    +epsyz'*masks*wdetJ*syz;

                [indi,indj,val]=find(Kx);
                Nddls=size(Kp,1);
                keep=~((indi<=Nddls)&(indj<=Nddls));
                Kx=sparse(indi,indj,val.*keep,Nddl_tot,Nddl_tot);
                if Nddl_tot>size(Kp,2)
                    Ko=sparse(Nddl_tot-size(Kp,1),Nddl_tot-size(Kp,2));
                    Kp=blkdiag(Kp,Ko);
                end
                Kp=Kp+Kx;
                Kef=Kp;
            else
                Kp=Kef;
            end

            if isfield(param0,'cz_length')&&isfield(param,'interface_model')
                load(fullfile('TMP',sprintf('%d_dphin_%d',nmod,10*(iscale-1))),'dphincz');
                load(fullfile('TMP',sprintf('%d_dphit_%d',nmod,10*(iscale-1))),'dphitcz');
                load(sprintf('TMP/%d_iphi_%d',nmod,10*(iscale-1)),'dz');
                [tn,tt]=ComputeTangentCohesiveStress(dphincz,dphitcz,nmod);
                dz=diag(sparse(dz));
                Kp=Kp+dphincz'*dz*tn+dphitcz'*dz*tt;
            end
        else
            if restart
                if dflag
                    load(fullfile('TMP',sprintf('%d_3d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_3d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
                    load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,10*(iscale-1))),'epsxy');
                    load(fullfile('TMP',sprintf('%d_3d_epszz_%d',nmod,10*(iscale-1))),'epszz');
                    load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,10*(iscale-1))),'epsxz');
                    load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,10*(iscale-1))),'epsyz');
                    maskg=1;
                else

                    load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
                    load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
                    load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'epsxy');
                    epsxz=sparse(size(epsxx,1),size(epsxx,2));epsyz=epsxz;epszz=epsxz;
                    load(fullfile('TMP',sprintf('%d_mask_%d',nmod,(iscale-1))),'maskg');
                end
                nmodo=nmod;
                iscaleo=iscale;
            end
            [sxx,syy,sxy,sxz,syz,szz]=ComputeTangentStress(epsxx,epsyy,epsxy,nmod,epsxz,epsyz,epszz,eflag);

            
            Kp= epsxx'*maskg*wdetJ*sxx...
                +epsyy'*maskg*wdetJ*syy...
                +epszz'*maskg*wdetJ*szz...
                +epsxy'*maskg*wdetJ*sxy...
                +epsxz'*maskg*wdetJ*sxz...
                +epsyz'*maskg*wdetJ*syz;
            if nlflag
                load(fullfile('TMP',sprintf('%d_matmod',nmod)),'S');
                if numel(S)>1
                    load(fullfile('TMP',sprintf('%d_3d_epsxz_%d',nmod,10*(iscale-1))),'Uxz','Uzx');
                    load(fullfile('TMP',sprintf('%d_3d_epsxy_%d',nmod,10*(iscale-1))),'Uxy','Uyx');
                    load(fullfile('TMP',sprintf('%d_3d_epsyz_%d',nmod,10*(iscale-1))),'Uzy','Uyz');
                    btsx=diag(S(:,1))*epsxx+diag(S(:,4))*Uxy+diag(S(:,6))*Uxz;
                    btsy=diag(S(:,4))*epsxx+diag(S(:,2))*Uxy+diag(S(:,5))*Uxz;
                    btsz=diag(S(:,6))*epsxx+diag(S(:,5))*Uxy+diag(S(:,3))*Uxz;

                    Kgx=epsxx'*maskg*wdetJ*btsx+Uxy'*maskg*wdetJ*btsy+Uxz'*maskg*wdetJ*btsz;

                    btsx=diag(S(:,1))*Uyx+diag(S(:,4))*epsyy+diag(S(:,6))*Uyz;
                    btsy=diag(S(:,4))*Uyx+diag(S(:,2))*epsyy+diag(S(:,5))*Uyz;
                    btsz=diag(S(:,6))*Uyx+diag(S(:,5))*epsyy+diag(S(:,3))*Uyz;

                    Kgy=Uyx'*maskg*wdetJ*btsx+epsyy'*maskg*wdetJ*btsy+Uyz'*maskg*wdetJ*btsz;

                    btsx=diag(S(:,1))*Uzx+diag(S(:,4))*Uzy+diag(S(:,6))*epszz;
                    btsy=diag(S(:,4))*Uzx+diag(S(:,2))*Uzy+diag(S(:,5))*epszz;
                    btsz=diag(S(:,6))*Uzx+diag(S(:,5))*Uzy+diag(S(:,3))*epszz;

                    Kgz=Uzx'*maskg*wdetJ*btsx+Uzy'*maskg*wdetJ*btsy+epszz'*maskg*wdetJ*btsz;
                    Kp=Kp+Kgx+Kgy+Kgz;
                end
            else
                if isfield(param.material_parameters,'lc')
             load(fullfile('TMP',sprintf('%d_epsijk_%d',nmod,10*(iscale-1))),'epsxxx','epsxxy','epsxyx','epsxyy','epsyyx','epsyyy');
                   load(fullfile('TMP',sprintf('%d_matmod',nmod)),'mu','lambda','lc');
     ltrx=lc^2*lambda*(epsxxx+epsyyx);
    txxx=2*lc^2*mu*epsxxx+ltrx;
    tyyx=2*lc^2*mu*epsyyx+ltrx;
    txyx=lc^2*mu*epsxyx;
    ltry=lc^2*lambda*(epsxxy+epsyyy);
    txxy=2*lc^2*mu*epsxxy+ltry;
    tyyy=2*lc^2*mu*epsyyy+ltry;
    txyy=lc^2*mu*epsxyy;
                  Kg=epsxxx'*maskg*wdetJ*txxx+epsyyx'*maskg*wdetJ*tyyx+epsxyx'*maskg*wdetJ*txyx+...
        epsxxy'*maskg*wdetJ*txxy+epsyyy'*maskg*wdetJ*tyyy+epsxyy'*maskg*wdetJ*txyy;
     Kp=Kp+Kg;
                end
                
                
                
            end

        end
    case 'abaqus'
        abqcmd='abaqus';
        if isfield(param0,'abaqus_command_line')
            abqcmd=param0.abaqus_command_line;
        end
        if nargin<3,step=[];end
        if isempty(step)
            restart=0;
            filres=sprintf('%d_%s',nmod,param0.result_file);
        else
            filres=sprintf('%d_%s_%d',nmod,param0.result_file,step);
            if step==1
                restart=0;
            else
                restart=1;
            end
        end
        if isempty(doinit)
            writeAbaqus([],nmod,step,restart,'init');
            doinit=false;
        end
        writeAbaqus([],nmod,step,restart,'K');

        if restart
            filresold=sprintf('%d_%s_%d',nmod,param0.result_file,step-1);
            unix([abqcmd,' job=',filres,' oldjob=',filresold,' double interactive']);
        else
            unix([abqcmd,' job=',filres,' double  interactive']);
        end
        Kp=readAbaqus(nmod,step,'K');
        %             load(fullfile('TMP',sprintf('%d_epsxx_%d',nmod,10*(iscale-1))),'epsxx','wdetJ');
        %             load(fullfile('TMP',sprintf('%d_epsyy_%d',nmod,10*(iscale-1))),'epsyy');
        %             load(fullfile('TMP',sprintf('%d_epsxy_%d',nmod,10*(iscale-1))),'epsxy');
        %             load(fullfile('TMP',sprintf('%d_mask_%d',nmod,10*(iscale-1))),'maskg');
        %             nmodo=nmod;
        %         [sxx,syy,sxy]=ComputeTangentStress(epsxx,epsyy,epsxy,nmod);
        %         maskg
        %         wdetJ
        %         Kpo=sxx'*maskg*wdetJ*epsxx+syy'*maskg*wdetJ*epsyy+sxy'*maskg*wdetJ*epsxy;
        %         'matlab'
        %         full(Kpo)
        %         'abaqus'
        %         full(Kp)
        %                 load(sprintf('TMP/%d_3d_epsxx_%d',nmod,10*(iscale-1)),'epsxx','wdetJ');
        %                 load(sprintf('TMP/%d_3d_epsyy_%d',nmod,10*(iscale-1)),'epsyy');
        %                 load(sprintf('TMP/%d_3d_epsxy_%d',nmod,10*(iscale-1)),'epsxy');
        %                 load(sprintf('TMP/%d_3d_epszz_%d',nmod,10*(iscale-1)),'epszz');
        %                 load(sprintf('TMP/%d_3d_epsxz_%d',nmod,10*(iscale-1)),'epsxz');
        %                 load(sprintf('TMP/%d_3d_epsyz_%d',nmod,10*(iscale-1)),'epsyz');
        %                 load(sprintf('TMP/%d_mask_%d',nmod,10*(iscale-1)),'maskg');
        %
        %  matmod=param.material_parameters;
        %
        %             mu=matmod.mu;
        %             lambda=matmod.lambda;
        %             ltr=(lambda)*(epsxx+epsyy+epszz);
        %                     sxx=2*mu*epsxx+ltr;
        %                     syy=2*mu*epsyy+ltr;
        %                     szz=2*mu*epszz+ltr;
        % sxy=mu*epsxy;
        % sxz=mu*epsxz;
        % syz=mu*epsyz;
        %             Kpo= epsxx'*maskg*wdetJ*sxx...
        %                 +epsyy'*maskg*wdetJ*syy...
        %                 +epszz'*maskg*wdetJ*szz...
        %                 +epsxy'*maskg*wdetJ*sxy...
        %                 +epsxz'*maskg*wdetJ*sxz...
        %                 +epsyz'*maskg*wdetJ*syz;
        %                         'matlab'
        %         full(Kpo(1:8,1:8))
        %         'abaqus'
        %         full(Kp(1:8,1:8))
        %         figure
        %         imagesc(full(abs(Kp-Kpo)))
        % error
        %             case 'ansys'
        %                  if nargin<3,step=1;end
        %             jobname=sprintf('%s_k',param0.result_file);
        %             filres=sprintf('%d_%s_k_%d',nmod,param0.result_file,step);
        %           if step==1
        %             restart=0;
        %          writeAnsys([],nmod,step,restart,'init');
        %         else
        %                 restart=1;
        %         end
        %                   writeAnsys([],nmod,step,restart,'K');
        %             unix(['/opt/ansys_inc/v121/ansys/bin/ansys121  -p AA_T_I  -j ',jobname,' -s read -l en-us -b < ',filres,'.inp | tee -i ',filres,'.out' ]);
        %  %           ['/opt/ansys_inc/v121/ansys/bin/ansys121  -p AA_T_I  -j ',jobname,' -s read -l en-us -b < ',filres,'.inp | tee -i ',filres,'.out' ]
        %  Kp=readAnsys(nmod,step,'K');
        % %                          load(sprintf('TMP/%d_3d_epsxx_%d',nmod,10*(iscale-1)),'epsxx','wdetJ');
        % %                         load(sprintf('TMP/%d_3d_epsyy_%d',nmod,10*(iscale-1)),'epsyy');
        % %                         load(sprintf('TMP/%d_3d_epsxy_%d',nmod,10*(iscale-1)),'epsxy');
        % %                         load(sprintf('TMP/%d_3d_epszz_%d',nmod,10*(iscale-1)),'epszz');
        % %                         load(sprintf('TMP/%d_3d_epsxz_%d',nmod,10*(iscale-1)),'epsxz');
        % %                         load(sprintf('TMP/%d_3d_epsyz_%d',nmod,10*(iscale-1)),'epsyz');
        % %                         load(sprintf('TMP/%d_mask_%d',nmod,10*(iscale-1)),'maskg');
        % %
        % %          matmod=param.material_parameters;
        % %
        % %                     mu=matmod.mu;
        % %                     lambda=matmod.lambda;
        % %                     ltr=(lambda)*(epsxx+epsyy+epszz);
        % %                             sxx=2*mu*epsxx+ltr;
        % %                             syy=2*mu*epsyy+ltr;
        % %                             szz=2*mu*epszz+ltr;
        % %         sxy=mu*epsxy;
        % %         sxz=mu*epsxz;
        % %         syz=mu*epsyz;
        % %                     Kpo= epsxx'*maskg*wdetJ*sxx...
        % %                         +epsyy'*maskg*wdetJ*syy...
        % %                         +epszz'*maskg*wdetJ*szz...
        % %                         +epsxy'*maskg*wdetJ*sxy...
        % %                         +epsxz'*maskg*wdetJ*sxz...
        % %                         +epsyz'*maskg*wdetJ*syz;
        % %                                 'matlab'
        % %                 full(Kpo(1:8,1:8))
        % %                                 'ansys'
        % % full(Kp(1:8,1:8))
        % %
        % %           error
        %
        %
        %         %         %             filres=sprintf('%d_%s',nmod,param0.result_file);
        %         %         % %            filres=sprintf('%d_%s_%d',nmod,param0.result_file,step);
        %         %         %
        %         %         %                 if step==1
        %         %         %             restart=0;
        %         %         %         else
        %         %         %                 restart=1;
        %         %         %         end
        %         %         %         if isempty(doinit)
        %         %         %             writeAnsys([],nmod,step,restart,'init');
        %         %         %             doinit=false;
        %         %         %         end
        %         %         %         writeAnsys([],nmod,step,restart,'K');
        %         %         %
        %         %         %         if restart
        %         %         %  error('todo')
        %         %         %         else
        %         %         %             toto=['/opt/ansys_inc/v121/ansys/bin/ansys121  -p AA_T_I  -j ',filres,' -s read -l en-us -b < ',filres,'.inp | tee -i ',filres,'.out' ];
        %         %         % %           /opt/ansys_inc/v110/ansys/bin/ansys110"  -p AA_T_I -dir "/home/jrethore/ansys-start" -j "toto" -s read -l en-us -b < "/home/jrethore/ansys-start/file.dat" | tee -i '/home/jrethore/ansys-start/file.out'
        %         %         % %unix(['''/opt/ansys_inc/v110/ansys/bin/ansys110''  -p AA_T_I  -j ''',path1,filres,''' -s read -l en-us -b < ''',path1,filres,'.inp'' | tee -i ''',path1,filres,'.out''' ]);
        %         %         % system(toto);
        %         %         %         end
        %         %         %         Kp=readAnsys(nmod,step,'K');
        %         %         end
    otherwise
        error('Unknown helper');
end
Kp=Kp*pix2m;
end


