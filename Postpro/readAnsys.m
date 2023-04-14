function [data]=readAnsys(nmod,step,action)
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
dflag=isfield(param,'extrusion_parameters');
cpflag=false;
if isfield(param,'plane_stress')
    cpflag=param.plane_stress;
end
switch action
    case 'S'
        error('NOT CODED YET');
    case 'Fint'
        fwid = fopen('nodalForces.txt','r');
        data=fscanf(fwid,'%e');
        data=-data(:);

        fclose(fwid);

    case 'K'
filres=sprintf('%d_%s_k_%d.out',nmod,param0.result_file,step);
        filresh=strrep(filres,'.out',sprintf('.matrix'));
        filres1=strrep(filres,'.out',sprintf('.txt'));
        filresm=strrep(filres,'.out',sprintf('.mapping'));
        %%
        fwid = fopen(filresm,'r');
        toto=textscan(fwid,'%s',4);
        if dflag
            toto=fscanf(fwid,'%d%d UX\n %d%d UY\n %d%d UZ\n');
        else
            toto=fscanf(fwid,'%d%d UX\n %d%d UY\n');
        end
        fclose(fwid);
        toto=reshape(toto,2,length(toto)/2)';
        ind=toto(:,1);
        node=toto(:,2);
        if dflag
        comp=repmat((1:3)',length(ind)/3,1);
        else
        comp=repmat((1:2)',length(ind)/2,1);
        end
        indm=node+max(node)*(comp-1);
        P=sparse(ind(:),indm(:),1,max(ind),max(ind));
        %%
        fwid = fopen(filres1,'r');

        toto=textscan(fwid,'%1c%*[^\n]',1);
        toto=textscan(fwid,'%1c%*[^\n]',1);
        sized=fscanf(fwid,'%d%d');
        toto=textscan(fwid,'%1c%*[^\n]',1);
        toto=textscan(fwid,'%1c%*[^\n]',1);
        sizec=fscanf(fwid,'%d',1);
        indc=fscanf(fwid,'%d');

        toto=textscan(fwid,'%1c%*[^\n]',1);
        toto=textscan(fwid,'%1c%*[^\n]',1);
        sizel=fscanf(fwid,'%d',1);
        indl=fscanf(fwid,'%d');
        indc=[1;1+cumsum(diff(indl)<=0)];
        toto=textscan(fwid,'%1c%*[^\n]',1);
        toto=textscan(fwid,'%1c%*[^\n]',1);
        sizev=fscanf(fwid,'%d',1);
        val=fscanf(fwid,'%e');

        fclose(fwid);
        %%
        %         datao = hb_to_msm ( filresh );
        data=sparse(indl,indc,val,sized(1),sized(2));
        data=data+data'-diag(diag(data));
        data=P'*data*P;

    otherwise
        error('UNKNOWN ACTION')

end
end
