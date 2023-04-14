function ExportImageToVTK(filreso)

[pp,filres,ext]=fileparts(filreso);
if isempty(ext) filreso=[filreso,'.mat'];end
load(filreso,'-mat','param')

videok=~isfield(param,'deformed_image');
if videok
    reader=VideoReader(param.reference_image);
    toto=readim(reader,1);
else
    toto=readim(param.reference_image);
end
toto=mean(toto,3);
xi=1:size(toto,1);
yi=1:size(toto,2);

delete([fullfile('VTK',[filres,'-im']),'-0*.vtk']);
toto=(double(toto)-min(toto(:)))/(max(toto(:))-min(toto(:)))*255;

VTKExportIntMap(filres,'im',xi,yi,uint8(toto));
movefile(fullfile('VTK',[filres,'-im.vtk']),fullfile('VTK',sprintf('%s-im-%06d.vtk',filres,0)));
if isfield(param,'image_number')
    images=param.image_number;
else
    if videok
        nbf=reader.NumberOfFrames-1;
        if isfield(param,'number_of_frames')
            nbf=param.number_of_frames;
        end
        dim=1;
        if isfield(param,'video_sampling')
            dim=param.video_sampling;
        end
        images=(2:dim:nbf)-1;
    else
        if iscell(param.deformed_image)
            images=1:size(param.deformed_image,2);
        else
            images=1;
        end
    end
end
for iim=1:length(images)
    if videok
        toto=readim(reader,images(iim)+1);
    else
        if length(images)>1
            toto=readim(param.deformed_image{iim});
        else
            toto=readim(param.deformed_image);
        end
    end
    toto=mean(toto,3);
toto=(double(toto)-min(toto(:)))/(max(toto(:))-min(toto(:)))*255;
    VTKExportIntMap(filres,'im',xi,yi,uint8(toto));
    movefile(fullfile('VTK',[filres,'-im.vtk']),fullfile('VTK',sprintf('%s-im-%06d.vtk',filres,images(iim))));
end

end