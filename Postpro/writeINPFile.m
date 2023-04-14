function writeINPFile(filreso,param,model)
[pp,filres,ext]=fileparts(filreso);
if isempty(ext) filreso=[filres,'.inp'];end

fid=fopen(filreso,'w');
fprintf(fid,'%% UFRECKLES INPut File\n');
fprintf(fid,'%% by Julien Rethore\n');
fprintf(fid,'%% Julien.Rethore@ec-nantes.fr\n');
fprintf(fid,'%% LaMCoS, INSA LYON, CNRS, UMR 5259\n');
fprintf(fid,'%% GeM, EC Nantes, CNRS, UMR 6183\n\n\n');
fprintf(fid,'%% Main parameters\n');
fields=fieldnames(param);
buff={};
for ii=1:length(fields)
    switch fields{ii}
        case 'analysis'
            buff{1}=sprintf('%%    Analysis type\n    param.analysis=''%s'';\n',param.analysis);
        case 'reference_image'
            buff{2}=sprintf('%%    Reference image file name\n    param.reference_image=''%s'';\n',param.reference_image);
        case 'deformed_image'
            if iscell(param.deformed_image)
                
                tmp=sprintf('%%    Deformed images file names\n    param.deformed_image={...\n');
                for ijm=1:size(param.deformed_image,1)
                    for iim=1:size(param.deformed_image,2)
                        tmp=[tmp,sprintf('                           ''%s''',param.deformed_image{ijm,iim})];
                        if ~(iim==size(param.deformed_image,2))
                            tmp=[tmp,sprintf(',...\n')];
                        end
                    end
                    if ~(ijm==size(param.deformed_image,1))
                        tmp=[tmp,sprintf(';...\n')];
                    else
                        tmp=[tmp,sprintf('};\n')];
                    end
                end
%                tmp=[tmp,sprintf('                      };\n')];
                
                buff{3}=tmp;
            else
                buff{3}=sprintf('%%    Deformed image file name\n    param.deformed_image=''%s'';\n',param.deformed_image);
            end
        case 'result_file'
            buff{4}=sprintf('%%    Result file name\n    param.result_file=''%s'';\n',param.result_file);
        case 'roi'
            sroi=sprintf('%d,',param.roi);
            sroi(end)=[];
            buff{5}=sprintf('%%    Region of interest\n    param.roi=[%s];\n',sroi);
        case 'pixel_size'
            buff{6}=sprintf('%%    Physical pixel size\n    param.pixel_size=%e;\n',param.pixel_size);
        case 'restart'
            if isfield(param,'time_step'), param.restart=2;end
            buff{7}=sprintf('%%    Type of analysis for multiple deformed images\n%%    0: sequential 1: independent 2: space-time\n    param.restart=%d;\n',param.restart);
            if param.restart==2
                buff{7}=[buff{7},sprintf('    param.time_step=%f;\n',param.time_step)];
            end
        case 'convergance_limit'
            buff{8}=sprintf('%%    Convergence criterion\n    param.convergence_limit=%e;\n',param.convergance_limit);
        case 'iter_max'
            buff{9}=sprintf('%%    Maximum iteration number\n    param.iter_max=%d;\n',param.iter_max);
        case 'do_pgd_prediction'
            buff{10}=sprintf('%%    For sequential analysis\n%%    0: first order explicit 1: automatic\n    param.do_pgd_prediction=%d;\n',param.do_pgd_prediction);
        case 'regularization_type'
            buff{11}=sprintf('%%    Regularization\n%%    ''none'': no regularization ''tiko'': Tichonov ''equilibrium_gap'': Equilibrium gap ''median'': continuous median filtering\n    param.regularization_type=''%s'';\n',param.regularization_type);
        case 'regularization_parameter'
            buff{12}=sprintf('%%    Cut-off wave length of the regularization (in pixel), except for ''median'', number of neighboors\n    param.regularization_parameter=%e;\n',param.regularization_parameter);
        case 'normalize_grey_level'
            buff{13}=sprintf('%%    Normalization of grey levels\n%%    0: global 1: element-wise\n    param.normalize_grey_level=%d;\n',param.normalize_grey_level);
        case 'number_of_frames'
            buff{14}=sprintf('%%    Number of frame to analyze\n    param.number_of_frames=%d;\n',param.number_of_frames);
        case 'video_sampling'
            buff{15}=sprintf('%%    Sampling of video file\n    param.video_sampling=%d;\n',param.video_sampling);
        case 'psample'
            buff{16}=sprintf('%%    Pixel skip\n    param.psample=%d;\n',param.psample);
        case 'stack_size'
            buff{17}=sprintf('%%    Stack size in case of 3D raw images\n    param.stack_size=[%d,%d,%d];\n',param.stack_size);
            
    end
end
for ii=1:length(buff)
    fprintf(fid,'%s',buff{ii});
end
fprintf(fid,'%% Model parameters\n');
fields=fieldnames(model);
buff={};
for ii=1:length(fields)
    switch fields{ii}
        case 'basis'
            buff{1}=sprintf('%%    Basis function\n%%    ''fem'': finite element\n%%    ''uni'': strain cage\n%%    ''beam'': 4-point bending\n%%    ''beam-nurbs'': B-splines Euler-Bernoulli flexural beam elements\n    model.basis=''%s'';\n',model.basis);
        case 'nscale'
            buff{2}=sprintf('%%    Number of coarse graining scales\n    model.nscale=%d;\n',model.nscale);
        case 'element_type'
            buff{3}=sprintf('%%    Element type for finite element basis\n%%    4: quadrangles (default), 3: triangles\n    model.element_type=%d;\n',model.element_type);
        case 'mesh_size'
            smesh=sprintf('%f,',model.mesh_size);
            smesh(end)=[];
            buff{4}=sprintf('%%    Mesh size (in pixel) along x and y for finite element basis\n    model.mesh_size=[%s];\n',smesh);
        case 'mesh_file'
            buff{5}=sprintf('%%    VTK input mesh file for finite element basis\n    model.mesh_file=''%s'';\n',model.mesh_file);
        case 'gluing_parameters'
            tmp=sprintf('%%    Parameters for gluing the VTK mesh to the reference frame\n');
            for ig=1:length(model.gluing_parameters)
                prm=model.gluing_parameters{ig};
                switch prm{1}
                    case 'translate'
                        tmp=[tmp,sprintf('%%    Translation (in pixel)\n    model.gluing_parameters{%d}={''%s'',[%f,%f,%f]};\n',ig,prm{1},prm{2})];
                    case 'rotate'
                        tmp=[tmp,sprintf('%%    Rotation (axis, center, angle)\n    model.gluing_parameters{%d}={''%s'',[%f,%f,%,f],[%f,%f,%,f],%f};\n',ig,prm{1},prm{2},prm{3},prm{4})];
                    case 'scale'
                        tmp=[tmp,sprintf('%%    Scale (center, factor)\n    model.gluing_parameters{%d}={''%s'',[%f,%f,%,f],%f};\n',ig,prm{1},prm{2},prm{3})];
                end
            end
            buff{6}=tmp;
        case 'nb_element'
            buff{7}=sprintf('%%    Number of B-splines flexural elements\n    model.nb_element=%d;\n',model.nb_element);
        case 'degree'
            buff{8}=sprintf('%%    Degree of B-splines flexural elements\n    model.degree=%d;\n',model.degree);
        case 'exx'
            buff{9}=sprintf('%%    To add an axial strain field to beam kinematics\n    model.exx=%d;\n',model.ess);
            
    end
end
for ii=1:length(buff)
    fprintf(fid,'%s',buff{ii});
end
fclose(fid);
end