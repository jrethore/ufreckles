function [data]=readAbaqus(nmod,step,action)
load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');
dflag=isfield(param,'extrusion_parameters');
cpflag=false;
if isfield(param,'plane_stress')
    cpflag=param.plane_stress;
end

if nargin<2
    filres=sprintf('%d_%s.dat',nmod,param0.result_file);
else
    if isempty(step)
        filres=sprintf('%d_%s.dat',nmod,param0.result_file);
    else
        filres=sprintf('%d_%s_%d.dat',nmod,param0.result_file,step);
    end
end
switch action
    case 'S'
        fwid = fopen(filres,'r');
        titi=[];
        while isempty(findstr(titi,'MAIN'))%ELEMENT  PT FOOT'))
            toto=textscan(fwid,'%4c%*[^\n]',1);
            titi=cell2mat(toto);
        end
        toto=textscan(fwid,'%1c%*[^\n]',1);
        toto=textscan(fwid,'%1c%*[^\n]',1);

        if dflag
            data=fscanf(fwid,'%d%d%f%f%f%f%f%f');
        data=reshape(data',8,length(data)/8)';
        data=data(:,3:8);
        else
            if cpflag
            data=fscanf(fwid,'%d%d%f%f%f');
        data=reshape(data',5,length(data)/5)';
        data=data(:,3:5);
            else
            data=fscanf(fwid,'%d%d%f%f%f%f');
        data=reshape(data',6,length(data)/6)';
        data=data(:,[3,4,6]);
            end
        end
        fclose(fwid);
    case 'Fint'
        fwid = fopen(filres,'r');
        titi=[];
        while isempty(findstr(titi,'NODE FOOT'))
            toto=textscan(fwid,'%10c%*[^\n]',1);
            titi=cell2mat(toto);
        end
        toto=textscan(fwid,'%1c%*[^\n]',1);

        if dflag
            data=fscanf(fwid,'%d%f%f%f');
        data=reshape(data',4,length(data)/4)';
        data=data(:,2:4);
        else
            data=fscanf(fwid,'%d%f%f');
        data=reshape(data',3,length(data)/3)';
        data=data(:,2:3);
        end
        data=data(:);
        fclose(fwid);
    case 'K'
        filres=strrep(filres,'.dat',sprintf('_STIF%d.mtx',step));
        toto=dlmread(filres,',');
        Nnods=max(toto(:,1));
        indi=toto(:,1)+Nnods*(toto(:,2)-1);
        indj=toto(:,3)+Nnods*(toto(:,4)-1);
        val=toto(:,5);
        clear toto
        data=sparse(indi,indj,val,(2+dflag)*Nnods,(2+dflag)*Nnods);
        data=data+data'-diag(diag(data));
    otherwise
        error('UNKNOWN ACTION')

end
end
