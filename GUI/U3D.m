ana=2;
while ana>1
    clear param model
    ana=menu('3D module','Quit','Create new data set','Modify data set','Run job','Post-processing');
    pause(0.1)
    if ana >1
        
        switch ana
            
            case 2 %new data set
                
                
                [filref,path0,filterindex]=uigetfile({'*.mat','MAT files';'*.tif;*.tiff','3D tiff files';'*.raw;*.bin','3D binary image files'},'Select the reference image');
                if filterindex>0
                    %[filref,path0,filterindex]=uigetfile({'*.mat','MAT files';'*.tif;*.tiff','3D tiff files'},'Select the reference image');
                    cd(path0);
                    [fildef,path1]=uigetfile({'*.mat','MAT files';'*.tif;*.tiff','3D tiff files';'*.raw;*.bin','3D binary image files'},'Select the deformed image(s)','MultiSelect', 'on');
                    %[fildef,path1]=uigetfile({'*.mat','MAT files';'*.tif;*.tiff','3D tiff files'},'Select the deformed image(s)','MultiSelect', 'on');
                    [filres,pathres]=uiputfile('*.dat','Save data set as ?');
                    param.result_file=strrep(filres,'.dat','.res');
                    [~, filrefo, extref] = fileparts(filref);
                    switch filterindex
                        case 3
                            answer=inputdlg({'x ?','y ?','z ?'},...
                                'Image size',1);
                            sizeim=[eval(answer{1}),eval(answer{2}),eval(answer{3})];
                            param.stack_size=sizeim;
                        case 2
                            %                         tiff2mat(filref)
                            %                         filref=strrep(filref,extref,'.mat');
                            %                         load(filref,'sizeim');
                            data = imfinfo(filref, 'TIFF');
                            sizeim = [data(1).Width, data(1).Height, length(data)];
                            
                        case 1
                            load(filref,'sizeim');
                    end
                    
                    if iscell(fildef)
                        fildef=sort(fildef);
                        
                        %                     if filterindex==2
                        %                         for iim=1:length(fildef)
                        %                             [~, ~, extdef] = fileparts(fildef{iim});
                        %                             tiff2mat(fildef{iim})
                        %                             fildef{iim}=strrep(fildef{iim},extdef,'.mat');
                        %                         end
                        %                     end
                    else
                        %                     if filterindex==2
                        %                         [~, ~, extdef] = fileparts(fildef);
                        %                         tiff2mat(fildef)
                        %                         fildef=strrep(fildef,extdef,'.mat');
                        %                     end
                    end
                    
                    
                    param.reference_image=filref;
                    
                    
                    param.deformed_image=fildef;
                    
                    param.analysis='correlation';
                    param.pixel_size=1;
                    param.restart=1;
                    
                    answer=inputdlg({'Element size ?','Coarse graining ?','Pixel skip'},...
                        'Parameters',1,{'32','3','0'});
                    
                    mesh_size=str2num(answer{1})*ones(3,1);
                    nscale=str2num(answer{2});
                    param.psample=str2num(answer{3})+1;
                    model.mesh_size=mesh_size;
                    model.basis='fem';
                    
                    
                    
                    
                    
                    roidef=menu('ROI definition...','from entropy','from coordinates','import mesh');
                    close;pause(0.1);
                    try model=rmfield(model,'mesh_file');catch, end
                    
                    switch roidef
                        case 3
                            [meshfile,~,filterindex]=uigetfile({'*.vtk','VTK (ASCII) files'},'Select a mesh file');
                            model.mesh_file=meshfile;
                            [xo,yo,zo,conn,elt,selected]=ReadVTK3D(meshfile);
                            roi=[max(1,floor(min(xo))-1),min(sizeim(1),ceil(max(xo))+1),...
                                max(1,floor(min(yo))-1),min(sizeim(2),ceil(max(yo))+1),...
                                max(1,floor(min(zo))-1),min(sizeim(3),ceil(max(zo))+1)];
                        case 2
                            
                            answer=inputdlg({'xmin ?','xmax ?','ymin ?','ymax ?','zmin ?','zmax ?'},...
                                'Region of interest',1,{'1',num2str(sizeim(1)),'1',num2str(sizeim(2)),'1',num2str(sizeim(3))});
                            roi=[];
                            for ii=1:6
                                roi(ii)=str2num(answer{ii});
                            end
                            close;pause(0.1);
                        case 1
                            roi=[1,sizeim(1),1,sizeim(2),1,sizeim(3)];
                            param.roi=roi;
                            filent=sprintf('%s-%02d-entropy.mat',filrefo,mesh_size(1));
                            if ~exist(filent,'file')
                                LoadParameters(param);
                                model.nscale=1;
                                LoadParameters(model,0);
                                ReferenceImage(0);
                                LoadMeshes(0);
                                LoadMask(0);
                                S=GetEntropy3D(fullfile('TMP','0_3d_mesh_0'));
                                Sm=median3(S);
                                VTKExportScalarVol(sprintf('%s-%02d',filrefo,mesh_size(1)),'entropy',1:size(S,1),1:size(S,2),1:size(S,3),Sm,mesh_size(1));
                                load(fullfile('TMP','0_3d_mesh_0'),'xo','yo','zo','conn');
                                save(filent,'xo','yo','zo','conn','Sm');
                            else
                                load(filent,'xo','yo','zo','conn','Sm');
                                
                            end
                            thres=inputdlg({'smin ?'},'Entropy threshold',1,{num2str(min(Sm(:)))});
                            thres=eval(thres{1});
                            keep=Sm>thres;
                            model.entropy_threshold=thres;
                            
                            
                            
                            conn=conn(keep,:);
                            [pind,~,j1]=unique(conn);
                            conn=reshape(j1,size(conn));
                            xo=xo(pind)+roi(1)-1;
                            yo=yo(pind)+roi(3)-1;
                            zo=zo(pind)+roi(5)-1;
                            Nnodes=[length(xo),1,1];
                            Nelems=[size(conn,1),1,1];
                            elt=8*ones(size(conn,1),1);
                            filmesh=sprintf('%s-%02d-mesh.dat',filrefo,mesh_size(1));
                            save(filmesh,'xo','yo','zo','conn','elt','Nnodes','Nelems','thres')
                            roi=[max(1,floor(min(xo)-1)),min(sizeim(1),ceil(max(xo)+1)),max(1,floor(min(yo)-1)),min(sizeim(2),ceil(max(yo)+1)),max(1,floor(min(zo)-1)),min(sizeim(3),ceil(max(zo)+1))];
                            writeVTKmesh3D(filmesh)
                            model.mesh_file=strrep(filmesh,'.dat','.vtk');
                    end
                    model.nscale=nscale;
                    param.roi=roi;
                    
                    
                    reg=menu('Smoothing ?','None','Strain','Median');
                    switch reg
                        case 1
                            param.regularization_type='none';
                        case 2
                            param.regularization_type='tiko';
                            lc=inputdlg({'lc'},'Cut-off wave length',1,{num2str(model.mesh_size(1))});
                            param.regularization_parameter=eval(lc{1});
                            
                        case 3
                            param.regularization_type='median';
                            lc=inputdlg({'nb'},'Nb of neighboor',1,{'1'});
                            param.regularization_parameter=eval(lc{1});
                    end
                    
                    
                    param.iter_max=30;
                    param.convergance_limit=1.e-3;
                    save(filres,'param','model');
                    writeINPFile(strrep(filres,'.dat','.ufr'),param,model)
                    
                    display(sprintf('Data set %s saved...',filres))
                end
            case 3 %modifydata set
                [filres,pathres,filterindex]=uigetfile({'*.dat','Ufreckles input files (*.dat)';'*.ufr','Ufreckles input ascii files (*.ufr)'},'Load data set');
                if filterindex>0
                    cd(pathres);
                    switch filterindex
                        case 1
                            load(filres,'param','model','-mat')
                        case 2
                            [param,model]=readINPFile(filres);
                    end
                    [filres,pathres]=uiputfile('*.dat','Save data set as ?');
                    param.result_file=strrep(filres,'.dat','.res');
                    if isfield(param,'stack_size')
                        sizeim=param.stack_size;
                    else
                        [~, ~, extref] = fileparts(param.reference_image);
                        switch extref
                            case '.mat'
                                load(param.reference_image,'-mat','sizeim')
                            case {'.tif','.tiff'}
                                data = imfinfo(param.reference_image, 'TIFF');
                                sizeim = [data(1).Width, data(1).Height, length(data)];
                                
                        end
                    end
                    [~, filrefo, ~] = fileparts(param.reference_image);
                    answer=inputdlg({'Element size ?','Coarse graining ?','Pixel skip'},...
                        'Parameters',1,{num2str(model.mesh_size(1)),num2str(model.nscale),num2str(param.psample-1)});
                    
                    param.psample=str2num(answer{3})+1;
                    mesh_size=str2num(answer{1})*ones(3,1);
                    nscale=str2num(answer{2});
                    model.mesh_size=mesh_size;
                    
                    
                    roidef=menu('ROI definition...','from entropy','from coordinates','import mesh');
                    close;pause(0.1);
                    try model=rmfield(model,'mesh_file');catch, end
                    switch roidef
                        case 3
                            [meshfile,~,filterindex]=uigetfile({'*.vtk','VTK (ASCII) files'},'Select a mesh file');
                            model.mesh_file=meshfile;
                            [xo,yo,zo,conn,elt,selected]=ReadVTK3D(meshfile);
                            roi=[max(1,floor(min(xo))-1),min(sizeim(1),ceil(max(xo))+1),...
                                max(1,floor(min(yo))-1),min(sizeim(2),ceil(max(yo))+1),...
                                max(1,floor(min(zo))-1),min(sizeim(3),ceil(max(zo))+1)];
                        case 2
                            
                            answer=inputdlg({'xmin ?','xmax ?','ymin ?','ymax ?','zmin ?','zmax ?'},...
                                'Region of interest',1,{num2str(param.roi(1)),num2str(param.roi(2)),num2str(param.roi(3)),num2str(param.roi(4)),num2str(param.roi(5)),num2str(param.roi(6))});
                            roi=[];
                            for ii=1:6
                                roi(ii)=str2num(answer{ii});
                            end
                            try model=rmfield(model,'mesh_file');
                            catch
                            end
                            try model=rmfield(model,'entropy_threshold');
                            catch
                            end
                            close;pause(0.1);
                        case 1
                            roi=[1,sizeim(1),1,sizeim(2),1,sizeim(3)];
                            param.roi=roi;
                            filent=sprintf('%s-%02d-entropy.mat',filrefo,mesh_size(1));
                            if ~exist(filent,'file')
                                LoadParameters(param);
                                model.nscale=1;
                                LoadParameters(model,0);
                                ReferenceImage(0);
                                LoadMeshes(0);
                                LoadMask(0);
                                S=GetEntropy3D(fullfile('TMP','0_3d_mesh_0'));
                                Sm=median3(S);
                                
                                load(fullfile('TMP','0_3d_mesh_0'),'xo','yo','zo','conn','elt');
                                ExportVTKField3D(sprintf('%s-%02d',filrefo,mesh_size(1)),1,xo,yo,zo,elt,conn,{Sm(:)},{'entropy'});
                                delete(fullfile('VTK',sprintf('%s-%02d-%06d.vtk',filrefo,mesh_size(1),0)));
                                movefile(fullfile('VTK',sprintf('%s-%02d-%06d.vtk',filrefo,mesh_size(1),1)),fullfile('VTK',sprintf('%s-%02d.vtk',filrefo,mesh_size(1))));
                                %                                VTKExportScalarVol(sprintf('%s-%02d',filrefo,mesh_size(1)),'entropy',1:size(S,1),1:size(S,2),1:size(S,3),Sm,mesh_size(1));
                                save(filent,'xo','yo','zo','conn','Sm');
                            else
                                load(filent,'xo','yo','zo','conn','Sm');
                                
                            end
                            if isfield(model,'entropy_threshold')
                                thres=inputdlg({'smin ?'},'Entropy threshold',1,{num2str(model.entropy_threshold)});
                            else
                                thres=inputdlg({'smin ?'},'Entropy threshold',1,{num2str(min(Sm(:)))});
                            end
                            thres=eval(thres{1});
                            keep=Sm>thres;
                            
                            model.entropy_threshold=thres;
                            
                            conn=conn(keep,:);
                            [pind,~,j1]=unique(conn);
                            conn=reshape(j1,size(conn));
                            xo=xo(pind)+roi(1)-1;
                            yo=yo(pind)+roi(3)-1;
                            zo=zo(pind)+roi(5)-1;
                            Nnodes=[length(xo),1,1];
                            Nelems=[size(conn,1),1,1];
                            elt=8*ones(size(conn,1),1);
                            filmesh=sprintf('%s-%02d-mesh.dat',filrefo,mesh_size(1));
                            save(filmesh,'xo','yo','zo','conn','elt','Nnodes','Nelems')
                            roi=[max(1,floor(min(xo)-1)),min(sizeim(1),ceil(max(xo)+1)),max(1,floor(min(yo)-1)),min(sizeim(2),ceil(max(yo)+1)),max(1,floor(min(zo)-1)),min(sizeim(3),ceil(max(zo)+1))];
                            writeVTKmesh3D(filmesh);
                            model.mesh_file=strrep(filmesh,'.dat','.vtk');
                    end
                    model.nscale=nscale;
                    param.roi=roi;
                    
                    
                    reg=menu('Smoothing ?','None','Strain','Median');
                    switch reg
                        case 1
                            param.regularization_type='none';
                        case 2
                            param.regularization_type='tiko';
                            if isfield(param,'regularization_parameter')
                                lc=inputdlg({'lc'},'Cut-off wave length',1,{num2str(param.regularization_parameter)});
                            else
                                lc=inputdlg({'lc'},'Cut-off wave length',1,{num2str(model.mesh_size(1))});
                            end
                            param.regularization_parameter=eval(lc{1});
                            
                        case 3
                            param.regularization_type='median';
                            if isfield(param,'regularization_parameter')
                                lc=inputdlg({'nb'},'Nb of neighboor',1,{num2str(param.regularization_parameter)});
                            else
                                lc=inputdlg({'nb'},'Nb of neighboor',1,{'1'});
                            end
                            param.regularization_parameter=eval(lc{1});
                    end
                    
                    
                    save(filres,'param','model');
                    writeINPFile(strrep(filres,'.dat','.ufr'),param,model)
                    display(sprintf('Data set %s saved...',filres))
                end
            case 4 %run job
                [filres,pathres,filterindex]=uigetfile({'*.dat','Ufreckles input files (*.dat)';'*.ufr','Ufreckles input ascii files (*.ufr)'},'Load data set');
                if filterindex>0
                    cd(pathres);
                    switch filterindex
                        case 1
                            load(filres,'param','model','-mat')
                        case 2
                            [param,model]=readINPFile(filres);
                    end
                    run_fem_job_3D(param,model)
                end
                
            case 5 %Post-processing results
                [filres,pathres,filterindex]=uigetfile('*.res','Load result file');
                if filterindex>0
                    cd(pathres);
                    load(filres,'-mat','param','model')
                    postproVTK3D(filres,0,exist(strrep(filres,'.res','-error.res'),'file'));
                    doerror=menu('Compute residuals ?','Yes','No');
                    doim=menu('Export advected deformed image ?','Yes','No');
                    
                    if (doerror==1)||(doim==1)
                        answer=inputdlg({'xmin ?','xmax ?','ymin ?','ymax ?','zmin ?','zmax ?'},...
                            'Zone for computing residuals',1,{num2str(param.roi(1)),num2str(param.roi(2)),num2str(param.roi(3)),num2str(param.roi(4)),num2str(param.roi(5)),num2str(param.roi(6))});
                        close;pause(0.1);
                        zone=[];
                        for ii=1:6
                            zone(ii)=str2num(answer{ii});
                        end
                        load(filres,'-mat','U')
                        xi=zone(1):zone(2);
                        yi=zone(3):zone(4);
                        zi=zone(5):zone(6);
                        model.nscale=1;
                        LoadParameters(param);
                        LoadParameters(model,0);
                        ReferenceImage(0);
                        LoadMeshes(0);
                        LoadMask(0);
                        
                        
                        filreso=strrep(filres,'.res','');
                        for iim=1:size(U,2)
                            if doerror==1
                                [disc,dynamic]=GetDiscrepancy3D(U,0,zone,iim,1);
                                VTKExportScalarVol([filreso,sprintf('-error-%05d',iim)],'error',xi,yi,zi,disc,1);
                                save([filreso,sprintf('-%05d',iim),'-voxel-error.res'],'disc','zone','-v7.3');
                                disc=uint8(disc/100*dynamic);
                                fid=fopen(sprintf('%s-%05d-voxel-error-%d-%d-%d.raw',filreso,iim,size(disc)),'w');
                                fwrite(fid,disc(:));
                                fclose(fid);
                                clear disc
                            end
                            if doim==1
                                [im1,~]=GetDiscrepancy3D(U,0,zone,iim,0);
                                VTKExportScalarVol([filreso,sprintf('-advdefim-%05d',iim)],'im',xi,yi,zi,im1,1);
                                save([filreso,sprintf('-%05d',iim),'-advdefim.res'],'im1','zone','-v7.3');
                                if max(im1(:))<256
                                    im1=uint8(im1);
                                    fid=fopen(sprintf('%s-%05d-advdefim-%d-%d-%d.raw',filreso,iim,size(im1)),'w');
                                    fwrite(fid,im1(:));
                                    fclose(fid);
                                end
                            end
                        end
                    end
                end
        end
    end
end




