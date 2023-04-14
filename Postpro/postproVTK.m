function postproVTK(filreso,submean,do_error)
if nargin<2,submean=0;end
if nargin<3,do_error=1;end
nmod=1;
[pp,filres,ext]=fileparts(filreso);
if isempty(ext)
    filreso=[filreso,'.mat'];
ext='.mat';
end
tic
load(filreso,'-mat')
if ~exist('U','var')
    U=U1;
end
if ~exist('model1','var')
    model1=model;
end
if strcmp(model1.basis,'nurbs')||strcmp(model1.basis,'btri')
    postproNURBSVTK(filreso,submean);
    return
end
if ~exist('unmasked_nodes','var')
    unmasked_nodes=[];
else
    if ~isempty(unmasked_nodes)
        Nnodes=[length(unmasked_nodes),1,1];
    end
end
if isfield(param,'image_number')
    images=param.image_number;
else
    images=1:size(U,2);
end
if isfield(param,'pixel_size')
    pix2m=param.pixel_size;
else
    pix2m=1;
end
if strcmp(param.analysis,'mechanics')
    dfac=1;
    fac=1+0*pix2m;
else
    fac=1;dfac=1;
end
if isfield(param,'sampling_factor')
    psample=param.sampling_factor;
else
    psample=1;
end
thermo=0;
if isfield(param,'thermo')
    thermo=param.thermo;
end
roi=param.roi;
sizeim=[roi(2)-roi(1),roi(4)-roi(3)]+1;
%%
set.ascii=0;
set.remark=' computed by UFRECKLES';
nn=prod(Nnodes);
ne=length(elt);
neo=length(elt);
%load(fullfile('TMP',[num2str(nmod),'_phix_0']),'phix');
if isfield(model1,'nb_gauss_points')
    ng=model1.nb_gauss_points;
else
    ng=0;
end
if do_error&&~strcmp(param.analysis,'mechanics')
    if exist([filres,'-error',ext],'file')
        fide=fopen([filres,'-error',ext]);
    else
        fide=fopen(fullfile('TMP',[num2str(nmod),'_error_0.mat']));
    end
    erroronelt=fread(fide,1);
    dynamic=fread(fide,1);

    if ~erroronelt
                 phix=CreateFiniteElementBasis(filreso,sizeim,1,[]);
                 [wdetJ,inde]=GetWeigthDetJ(filreso,sizeim);
        M=phix'*wdetJ*phix;
        phix=wdetJ*phix;
        if ng>0
            load(fullfile('TMP',[num2str(nmod),'_phix_0']),'inde');
            gpc2cell=sparse(inde,1:length(inde),diag(wdetJ),ne,length(inde));
            acell=sum(gpc2cell,2);
        end
    end
end

clear wdetJ
rint=false;
if isfield(model,'reduced_integration')
    rint=model.reduced_integration;
end
ngq=4;
if rint,ngq=1;end

save(filreso,'rint','-append');
[dphidx,dphidy]=CreateGradFiniteElementBasis(filreso,[1,1],1,[],'Gauss_points');

phi0=sparse(size(dphidx,1),size(dphidx,2));
epsxx=[dphidx,phi0];
epsyy=[phi0,dphidy];
epsxy=[dphidy,dphidx];
Uxy=[dphidy,phi0];
Uyx=[phi0,dphidx];
% load(fullfile('TMP',[num2str(nmod),'_epsxx_0']),'Uxx','epsxx');
% load(fullfile('TMP',[num2str(nmod),'_epsyy_0']),'Uyy','epsyy');
% load(fullfile('TMP',[num2str(nmod),'_epsxy_0']),'Uyx','Uxy','epsxy');
indi=repmat((1:length(elt))',1,ngq);
val=repmat((elt==4)/ngq,1,ngq)+[(elt==3),repmat(0,length(elt),ngq-1)];
indj=ones(size(indi'));
found=find(val'>0);
indj(found)=1:length(found);
indj=indj';
foundt3=find(elt==3);
foundq4=find(elt==4);
foundt=[foundt3;foundq4];
val=val(foundt,:);
indi=indi(foundt,:);
indj=indj(foundt,:);
gp2cell=sparse(indi,indj,val);
selected_elt=1;
if exist('face_elts','var')||~isempty(unmasked_nodes)
    selected_elt=ones(length(elt),1);
    %   selected_elt(cz_elts)=1;
    %    selected_elt=diag(sparse(selected_elt));
    if exist('face_elts','var')
        assert(model1.shift_enrichment==1,'YOU MUST USE THE SHIFTED ENRICHEMNT FOR VTK POSTPROCESSING');
        selected_elt(face_elts)=0;
    end
    if ~isempty(unmasked_nodes)
        unmasked=zeros(length(xo),1);
        unmasked(unmasked_nodes)=1;
        unmasked=[unmasked;0];
        found=find(conn(:)==0);
        conn(found)=length(unmasked);
        masked_elts=find(sum(unmasked(conn),2)./elt<1);
        selected_elt(masked_elts)=0;
        id_unmasked=zeros(length(xo)+1,1);
        id_unmasked(unmasked_nodes)=1:length(unmasked_nodes);
        conn=id_unmasked(conn);
        xo=xo(unmasked_nodes);
        yo=yo(unmasked_nodes);
    end
    indj=find(selected_elt);
    indi=1:length(indj);
    selected_elt=sparse(indi,indj,1,length(indj),length(elt));
    elt=full(selected_elt*elt);
    conn=full(selected_elt*conn);
    foundt3=find(elt==3);
    foundq4=find(elt==4);
    ne=length(elt);
    postproVTKX(filres,submean);
end
if exist('U2','var')
    postproVTKW(filres,submean);
end
gp2cell=selected_elt*gp2cell;
if ~strcmp(param.analysis,'mechanics')&&do_error
    if ~erroronelt
if ng>0
    gpc2cell=selected_elt*gpc2cell;
    acell=sum(gpc2cell,2);
end
    end
end
epsxx=gp2cell*epsxx;
epsxy=gp2cell*epsxy;
epsyy=gp2cell*epsyy;
Uxx=epsxx;
Uxy=gp2cell*Uxy;
Uyx=gp2cell*Uyx;
Uyy=epsyy;
nt3=sum(elt==3);
nq4=sum(elt==4);

display(sprintf('Result file : %s',filres));
delete([fullfile('VTK',filres),'-0*.vtk']);
L=[1+0*xo,0*xo,yo;...
    0*yo,1+0*yo,-xo];
%L=[1+0*xo,0*xo;...
%    0*yo,1+0*yo];
images=[0,images];
U=[zeros(size(U,1),1),U];
for iim=1:size(U,2)
    set.vtkname=[filres , sprintf('-%06d.vtk',images(iim))];
    %% 
    if set.ascii
    fwid = fopen(fullfile('VTK',set.vtkname),'w');
    else
        fwid = fopen(fullfile('VTK',set.vtkname),'w','b'); % IMPORTANT: big endian
    end
    count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
    count = fprintf(fwid,[set.vtkname,set.remark,'\n']);
    if set.ascii
        count = fprintf(fwid,'ASCII\n');
    else
        count = fprintf(fwid,'BINARY\n');
    end
    count = fprintf(fwid,'DATASET UNSTRUCTURED_GRID\n');

    count = fprintf(fwid,'POINTS %u float\n',nn);
    %data=[xo+roi(1)-1;yo+roi(3)-1;0*yo]';
    data=[(xo-1)*psample+roi(1),(yo-1)*psample+roi(3),0*yo]';


    if set.ascii
        fprintf(fwid, '%f %f %f \n', dfac*data);
    else
        fwrite(fwid, data,'float');
    end
    count = fprintf(fwid,'CELLS %u %u\n',nt3+nq4,4*nt3+5*nq4);
    data=[repmat(3,1,length(foundt3));conn(foundt3,1:3)'-1];


    if set.ascii
        fprintf(fwid, '%d %d %d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end

    data=[repmat(4,1,length(foundq4));conn(foundq4,1:4)'-1];


    if set.ascii
        fprintf(fwid, '%d %d %d %d %d\n', data);
    else
        fwrite(fwid, data,'uint');
    end


    count = fprintf(fwid,'CELL_TYPES %u\n',ne);
    data=[repmat(5,1,nt3),repmat(9,1,nq4)];

    if set.ascii
        fprintf(fwid, '%d\n', data);
    else
        fwrite(fwid, data,'uint');
    end
    %
    Ui=U(1:(2*length(xo)),iim);
    if submean
    A=L\Ui;
    Ui=Ui-L*A;
%         Ui(1:nn)=Ui(1:nn)-mean(Ui(1:nn));
%         Ui(nn+(1:nn))=Ui(nn+(1:nn))-mean(Ui(nn+(1:nn)));
    end
    count = fprintf(fwid,'POINT_DATA %u\n',nn);
    count = fprintf(fwid,['VECTORS Displacement float\n']);
    UU=fac*[Ui(1:nn)';Ui(nn+(1:nn))';repmat(0,1,nn)];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
    count = fprintf(fwid,['VECTORS U float\n']);
    A=L\Ui;
    Ui=Ui-L*A;
    UU=fac*[Ui(1:nn)';Ui(nn+(1:nn))';repmat(0,1,nn)];
    if set.ascii
        fprintf(fwid, '%f %f %f \n', UU);
    else
        fwrite(fwid,UU,'float');
    end
    %  count = fprintf(fwid,['TENSORS Strain float\n']);
    %  Exx=M\(phix'*epsxx*Ui);
    %  Eyy=M\(phix'*epsyy*Ui);
    %  Exy=M\(0.5*phix'*epsxy*Ui);
    %  E0=0*Exx;
    %  E=[Exx,Exy,E0,Exy,Eyy,E0,E0,E0,E0]';
    %  reshape(E(:),3,3*ne);
    %
    % if set.ascii
    %    fprintf(fwid, '%f %f %f \n', E);
    % else
    %       fwrite(fwid, E,'float');
    % end
    if thermo
            count = fprintf(fwid,['SCALARS Temperature float, 1\n']);

            count = fprintf(fwid,'LOOKUP_TABLE default\n');
            UU=U((2*length(xo))+(1:length(xo)),iim)';
    if set.ascii
        fprintf(fwid, '%f\n', UU);
    else
        fwrite(fwid,UU,'float');
    end

    end
    if do_error&&~strcmp(param.analysis,'mechanics')
        if ~erroronelt&&ng==0
            if images(iim)==0
                disc=zeros(size(phix,1),1);
            else
            disc=fread(fide,size(phix,1));
            end
            disc=M\(phix'*(disc));
            disc=100*abs(disc)/dynamic;
            count = fprintf(fwid,['SCALARS Error float, 1\n']);
            count = fprintf(fwid,'LOOKUP_TABLE default\n');
            if set.ascii
                fprintf(fwid, '%f\n', disc);
            else
                fwrite(fwid,disc,'float');
            end
        end
    end
    count = fprintf(fwid,'CELL_DATA %u\n',ne);
    if iim==1&&~strcmp(param.analysis,'mechanics')&&do_error
        Selems=GetEntropy(filreso);
        count = fprintf(fwid,['SCALARS Entropy float,1\n']);
        count = fprintf(fwid,'LOOKUP_TABLE default\n');
        if set.ascii
            fprintf(fwid, '%f \n', Selems);
        else
            fwrite(fwid, Selems,'float');
        end
    end
    count = fprintf(fwid,['TENSORS Strain float\n']);
    Exx=epsxx*Ui;
    Eyy=epsyy*Ui;
    Exy=0.5*epsxy*Ui;
    E0=0*Exx;
    E=[Exx,Exy,E0,Exy,Eyy,E0,E0,E0,E0]';
    if set.ascii
        fprintf(fwid, '%f %f %f \n', E);
    else
        fwrite(fwid, E,'float');
    end
    count = fprintf(fwid,['TENSORS Green-Lagrange float\n']);
    Fxx=1+Uxx*Ui;Fxy=Uxy*Ui;
    Fyy=1+Uyy*Ui;Fyx=Uyx*Ui;
    FTFxx=Fxx.*Fxx+Fyx.*Fyx;
    FTFxy=Fxx.*Fxy+Fyx.*Fyy;
    FTFyy=Fxy.*Fxy+Fyy.*Fyy;
    E0=0*Exx;
    EGL=0.5*[(FTFxx-1),FTFxy,E0,FTFxy,(FTFyy-1),E0,E0,E0,E0]';

    if set.ascii
        fprintf(fwid, '%f %f %f \n', EGL);
    else
        fwrite(fwid, EGL,'float');
    end
%     count = fprintf(fwid,['TENSORS Hencky float\n']);
%     if iim==1
%         H=0*EGL;
%     else
%         H=2*EGL;
%         H(1,:)=H(1,:)+1;
%         H(3,:)=1;
%         H(5,:)=H(5,:)+1;
%         H(6,:)=1;
%         H(7,:)=1;
%         H(8,:)=1;
%         H(9,:)=H(9,:)+1;
%         H=0.5*log(H);
%     end
%     if set.ascii
%         fprintf(fwid, '%f %f %f \n', H);
%     else
%         fwrite(fwid, H,'float');
%     end

    if isfield(model1,'vtk_export')
        
        for iex=1:length(model1.vtk_export)
            go=0;
            switch model1.vtk_export{iex}
                case 'S'
                    go=1;
                    if iim==1
                        S=sparse(size(gp2cell,1),6);
                        vm=sparse(size(gp2cell,1),1);
                    else
                        load(sprintf('%s_%04d',strrep(param.result_file,'.res',''),iim-1),'S');
                        vm=gp2cell*sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2 ...
                            -S(:,1).*S(:,3)-S(:,2).*S(:,3)-S(:,1).*S(:,2) ...
                            +3*(S(:,4).^2+S(:,5).^2+S(:,6).^2));
                        S=gp2cell*S;
                    end


                    count = fprintf(fwid,['SCALARS VM double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%f\n', full(vm'));
                    else
                        fwrite(fwid, full(vm),'double');
                    end
                    Exx=S(:,1);Eyy=S(:,2);Ezz=S(:,3);
                    Exy=S(:,4);Eyz=S(:,5);Exz=S(:,6);
                    count = fprintf(fwid,['TENSORS Stress double\n']);
                case 'D'
                    go=0;
                    if iim==1
                        d=sparse(size(gp2cell,1),1);
                    else
                        load(sprintf('%s_%04d',param.result_file,iim-1),'d');
                        d=gp2cell*d;
                    end


                    count = fprintf(fwid,['SCALARS D double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%f\n', full(d'));
                    else
                        fwrite(fwid, d,'double');
                    end
                case 'K'
                    go=0;
                    if iim==1
                        eeqmax=sparse(size(gp2cell,1),1);
                    else
                        load(sprintf('%s_%04d',param.result_file,iim-1),'eeqmax');
                        eeqmax=gp2cell*eeqmax;
                    end


                    count = fprintf(fwid,['SCALARS K double, 1\n']);
                    count = fprintf(fwid,'LOOKUP_TABLE default\n');
                    if set.ascii
                        fprintf(fwid, '%f\n', full(eeqmax'));
                    else
                        fwrite(fwid, full(eeqmax),'double');
                    end
                    
            end
            if go
                E=[ Exx,Exy,Exz,...
                    Exy,Eyy,Eyz,...
                    Exz,Eyz,Ezz]';
                if set.ascii
                    fprintf(fwid, '%e %e %e \n', full(E));
                else
                    fwrite(fwid, full(E),'double');
                end
            end
        end
    end





    
    if do_error&&~strcmp(param.analysis,'mechanics')
        if erroronelt
            if images(iim)==0
                disc=zeros(neo,1);
            else
                disc=fread(fide,neo);
            end
            disc=100*(abs(disc))/dynamic;
            count = fprintf(fwid,['SCALARS Error float, 1\n']);
            count = fprintf(fwid,'LOOKUP_TABLE default\n');
            if set.ascii
                fprintf(fwid, '%f\n', disc);
            else
                fwrite(fwid,disc,'float');
            end
            
        else
            if ng>0
                if images(iim)==0
                    disc=zeros(neo,1);
                else
                    disc=fread(fide,neo);
                end
                disc=100*(gpc2cell*abs(disc))./acell/dynamic;
                count = fprintf(fwid,['SCALARS Error float, 1\n']);
                count = fprintf(fwid,'LOOKUP_TABLE default\n');
                if set.ascii
                    fprintf(fwid, '%f\n', disc);
                else
                    fwrite(fwid,disc,'float');
                end
            end
        end
    end
    fclose(fwid);
end
if ~strcmp(param.analysis,'mechanics')&&do_error
    fclose(fide);
end
fprintf(1,'vtkexport done in %5.3f s\n',toc);
