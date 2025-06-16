function varargout=Umaster(varargin)
% UMASTER MATLAB code for Umaster.fig
%      UMASTER, by itself, creates a GUI for UFreckles

% Begin initialization code - DO NOT EDIT
gui_Singleton=1;
gui_State=struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Umaster_OpeningFcn, ...
    'gui_OutputFcn',  @Umaster_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback=str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}]=gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT
% --- Executes just before Umaster is made visible.
function Umaster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Umaster (see VARARGIN)

% Choose default command line output for Umaster
handles.output=hObject;
warning('off');
if nargin<4
    logiciel=[1,2,3,4,5,6,7,8,9];
    version='Development release';
else
    logiciel=varargin{1};
    version=varargin{2};
    if logiciel==0
        logiciel=1;
        version='Demo version';
        set(handles.run_analysis,'Enable','off');
        set(handles.cmenu_mesh_elt_type_adapt,'Enable','off');
        set(handles.cmenu_mesh_load_mesh_load,'Enable','off');
        
    end
end
version='V 2.2';
handles.logiciel=logiciel;
set(handles.version_text,'String',version);
set(handles.message_text,'String','Welcome to UFreckles !');
% Update handles structure

% UIWAIT makes Umaster wait for user response (see UIRESUME)
% uiwait(handles.figure1);
ax1=annotation('arrow',0.008+[0,1/20]-0.002,0.08+[0,0]);
ax2=annotation('arrow',0.008+[0,0]-0.002,0.08+[0,1/20*3/2]);
set(ax1,'Color','blue','LineWidth',2,'HeadStyle','cback3');
set(ax2,'Color','blue','LineWidth',2,'HeadStyle','cback3');
tx1=annotation('textbox',[1/20-0.02,0.09,0.02,0.02],'String','X');
set(tx1,'Color','blue','LineStyle','none','FitBoxToText','on','FontSize',12);
tx2=annotation('textbox',[0.008-0.002,0.05+1/20*3/2,0.02,0.02],'String','Y');
set(tx2,'Color','blue','LineStyle','none','FitBoxToText','on','FontSize',12);
handles.axis=[ax1,ax2,tx1,tx2];

nb_cols_cscale=50;
strain_on_elts=1;
median_filter_cell=0;
if isdeployed
    try
        logoim=(fullfile(ctfroot,'ufreckles','logo.png'));
        image(permute(readim(logoim),[2,1,3]))
    catch
        logoim=(fullfile(ctfroot,'UFreckles','logo.png'));
        image(permute(readim(logoim),[2,1,3]))
    end
    try
        fid=fopen('ufreckles.par','r');
        buf=textscan(fid,'%s');
        buf=buf{1};
        for ii=1:length(buf)
            eval(buf{ii});
        end
        fclose(fid);
    catch
    end
else
    logoim='logo.png';
    image(permute(readim(logoim),[2,1,3]))
    try
        fid=fopen('ufreckles.par','r');
        buf=textscan(fid,'%s');
        buf=buf{1};
        for ii=1:length(buf)
            eval(buf{ii});
        end
        fclose(fid);
    catch
    end
end
axis off
axis xy
axis auto

cmap=[
    0.317647,0.341176,0.431373;...
    0,0,1;...
    0,1,1;...
    0,1,0;...
    1,1,0;...
    1,0,0;...
    0.878431,0,1;...
    ];
cmap=[
    0.0504    0.0298    0.5280;...
    0.2579    0.0136    0.6167;...
    0.4176    0.0006    0.6584;...
    0.5655    0.0537    0.6405;...
    0.6928    0.1651    0.5645;...
    0.7964    0.2780    0.4713;...
    0.8814    0.3925    0.3832;...
    0.9481    0.5152    0.2974;...
    0.9883    0.6523    0.2114;...
    0.9891    0.8064    0.1459;...
    0.9400    0.9752    0.1313;...
    ];
%handles.cmap=interp1(cmap,1:0.05:size(cmap,1));

handles.cmap=interp1(cmap,linspace(1,size(cmap,1),nb_cols_cscale));
handles.wfac=1;
handles.zfac=1;
handles.stereo=0;
handles.gim=0;
handles.sonelt=strain_on_elts;
handles.medcell=median_filter_cell;
handles.fsolver=0;
handles.umaster=handles.figure1;
handles.gpoint=[];
handles.gmesh=[];
handles.gedge=[];
handles.fbeam=0;
handles.fgage=0;
handles.fparams=0;
handles.ganot=[];
handles.showim=1;
handles.showcb=1;
handles.rmrbm=0;
handles.fstrain=0;
handles.vref=1;
handles.gcbar=0;
handles.erzone={};
handles.view=[0,90];
handles.scale_mode=ones(10,10);
set(handles.figure1,'uicontextmenu','');
set(handles.figure1,'WindowScrollWheelFcn','');
set(handles.figure1,'WindowButtonMotionFcn', @show_position);
set(handles.figure1,'WindowButtonDownFcn',@double_clic);

set(handles.axis,'Visible','off')
handles.pframe=1;
handles.wfac=1;
handles.zfac=1;
handles.field=1;
handles.minmax={};
handles.ominmax=repmat({[Inf,-Inf]},10,10);
handles.ucomp=0;
handles.ecomp=0;
handles.ccomp=0;
handles.ercomp=0;
handles.erroronelt=0;
set(handles.rotate_button,'Enable','on');
set(handles.switch_view_button,'Enable','off');
set(handles.switch_frame_button,'Enable','off');

handles=reset_param(handles);
handles=set_param(handles);
function move_cb(hObject,~)
if strcmp(get(gcf,'selectiontype'),'extend')
    
    set(hObject,'ButtonDownFcn','selectmoveresize');
end
set(hObject,'ButtonDownFcn',@move_cb);
function show_position(hObject,eventdata)
handles=guidata(hObject);
if (handles.ana>0)
    xy=getposition(handles);
end
function xy=getposition(handles,step)
if nargin<2,step=handles.gstep;end
set(0,'CurrentFigure',handles.figure1);
[az,el] = view;
dview=((az==0)&&(el==90));
buf='';
xy=get(handles.axes1,'CurrentPoint');
xy=xy(1,1:2);
xy(1)=step*round(xy(1)/step);
xy(2)=step*round(xy(2)/step);
roi=[1,handles.sizeim(1),1,handles.sizeim(2)];
if dview
    if inpolygon(xy(1),xy(2),roi([1,2,2,1,1,]),roi([3,3,4,4,3]))
        buf=sprintf('X: %0.1f   Y: %0.1f',xy(1),xy(2));
    end
end
set(handles.version_text,'String',buf);

function double_clic(hObject,eventdata)
handles=guidata(hObject);
if (handles.ana>0)
    if strcmp(handles.param.analysis,'correlation')
        if strcmp(get(handles.figure1,'selectiontype'),'open')%||(length(eventdata)==1&&eventdata{1})
            model=[];
            switch handles.fbasis
                case 'fem'
                    model=handles.fem_model;
                case 'uni'
                    model=handles.uni_model;
                case 'beam'
                    model=handles.beam_model;
                case 'vic'
                    %                    model=handles.vic_model;
            end
            if ~isempty(model)
                model.nscale=handles.nscale;
                handles.fparams=UParams(handles.param,model);
            end
        end
    end
    guidata(handles.figure1,handles)
end

function handles=reset_param(handles)
reset(gca);
clear functions
axis off
axis xy
axis auto
set(handles.figure1,'Name','UFreckles')
param.analysis='correlation';
param.pixel_size=1;
param.normalize_grey_level=0;
param.iter_max=50;
param.convergance_limit=1.e-3;
param.result_file='ufreckles.res';
param.onflight=1;
param.restart=1;
param.da=0;
param.do_pgd_prediction=0;
param.deformed_image=[];
param.regularization_parameter=0;
param.regularization_type='none';
param.detect=0;
param.ulive=0;
param.thermo=0;
param.vic=0;
param.psample=1;
fem_model.basis='fem';
fem_model.mesh_size=[16,16];
fem_model.nscale=1;
fem_model.element_type=4;
fem_model.mesh_type=1;
fem_model.smin=0.5;
fem_model.mask=1;
fem_model.zone={};
fem_model.phantom_nodes=0;

beam_model.zone={};
beam_model.exx=1;
beam_model.nscale=1;
beam_model.degree=3;
beam_model.nb_element=10;
beam_model.basis='beam';
beam_model.beam_type='euler';


vic_model.zone={};
vic_model.basis='vic';
vic_model.nscale=3;

uni_model.nscale=1;
uni_model.basis='uni';
uni_model.zone={};
% if isfield(handles,'ganot')
%     for ig=1:length(handles.ganot)
%         delete(handles.ganot{ig})
%     end
% end
try
    for iz=1:numel(handles.ganot)
        delete(handles.ganot{iz});
    end
catch
end
try delete(handles.gim);catch, end
try delete(handles.gpoint);catch, end
try delete(handles.gmesh);catch, end
try delete(handles.gedge);catch, end
for id=1:length(handles.fgage)
    try delete(handles.fgage(id));catch, end
    try delete(100+handles.fgage(id));catch, end
    try delete(200+handles.fgage(id));catch, end
    try delete(300+handles.fgage(id));catch, end
    try delete(400+handles.fgage(id));catch, end
    try delete(4000+handles.fgage(id));catch, end
    try delete(500+handles.fgage(id));catch, end
    try delete(600+handles.fgage(id));catch, end
end
try delete(handles.fbeam);catch, end
try delete(handles.fparams);catch, end
try delete(handles.gcbar);catch, end
try delete(handles.fsolver);catch, end
try delete((handles.calculator{3,:}));catch, end
try delete((handles.fem_model.zone{3,:}));catch, end
try delete((handles.fem_model.zone{11,:}));catch, end
try delete((handles.uni_model.zone{3,:}));catch, end
try delete((handles.beam_model.zone{3,:}));catch, end
try delete((handles.vic_model.zone{3,:}));catch, end
try delete((handles.erzone{3,:}));catch, end
try handles=rmfield(handles,'mvisu');catch, end
try handles=rmfield(handles,'reader');catch, end

handles.pframe=1;
handles.wfac=1;
handles.zfac=1;
handles.stereo=0;
handles.ganot={};
handles.gmesh=[];
handles.gpoint=[];
handles.gedge=[];
handles.fgage=0;
handles.fbeam=0;
handles.fsolver=0;
handles.fparams=0;
handles.gcbar=[];
handles.showim=1;
handles.erzone={};
handles.calculator={};
handles.uvisu=[];
handles.evisu=[];
handles.scroll_adjusted='';
handles.rmrbm=0;
handles.fstrain=0;
handles.vref=1;
handles.gstep=1;
handles.cracked=0;
handles.nmod=1;
handles.nscale=3;
handles.fbasis='none';
handles.param=param;
handles.fem_model=fem_model;
handles.beam_model=beam_model;
handles.uni_model=uni_model;
handles.vic_model=vic_model;
handles.ana=0;
handles.animation.iim=0;
handles.animation.nbstep=0;
handles.animation.frames=0;
handles.animation.playing=0;
handles.preview=0;
handles.ondefimage=1;
handles.showedge=0;
handles.upreview=[];
handles.field=1;
handles.ucomp=0;
handles.ecomp=0;
handles.ccomp=0;
handles.minmax={};
handles.ominmax=repmat({[Inf,-Inf]},10,10);
handles.erroronelt=0;
handles.ercomp=0;
handles.view=[0,90];
handles.showplot=0;
handles.scale_mode=ones(10,10);
set(handles.axis,'Visible','off')
set(handles.image_info_text,'Visible','off')
set(handles.field_text,'Visible','off')
set(handles.figure1,'uicontextmenu','');
set(handles.figure1,'WindowScrollWheelFcn','');
set(handles.figure1,'WindowButtonMotionFcn', @show_position);
set(handles.figure1,'WindowButtonDownFcn',@double_clic);
view(handles.view);
guidata(handles.figure1,handles)

set(handles.preview_button,'Enable','off');
set(handles.preview_button,'State','off');
set(handles.show_data_button,'Enable','off');
%set(handles.rotate_button,'State','off');
%set(handles.rotate_button,'Enable','off');
%set(handles.switch_view_button,'Enable','off');
set(handles.switch_frame_button,'Enable','off');
set(handles.cmenu_mesh_entropy,'Enable','on');

function update_param(hObject)
handles=guidata(hObject);
handles=set_param(handles);

function handles=set_param(handles)
try
    handles.usolver=guidata(handles.fsolver);
    handles.param.pixel_size=eval(get(handles.usolver.pix2m_edit,'String'));
    handles.param.psample=eval(get(handles.usolver.psample_edit,'String'));
    handles.param.iter_max=eval(get(handles.usolver.maxiter_edit,'String'));
    handles.param.convergance_limit=eval(get(handles.usolver.convergence_edit,'String'));
    handles.nscale=eval(get(handles.usolver.nscale_edit,'String'));
    handles.param.do_pgd_prediction=get(handles.usolver.do_pgd_edit,'Value')-1;
    handles.param.restart=~get(handles.usolver.restart_edit,'Value');
    handles.param.normalize_grey_level=get(handles.usolver.normalize_grey_levels_edit,'Value');
    
    if ~isfield(handles.param,'deformed_image')
        %     set(handles.menu_video,'Visible','on')
        %     set(handles.nbf_edit,'String',num2str(handles.param.number_of_frames));
        %     set(handles.nbf_edit,'TooltipString',sprintf('indice of the last frame (max=%d)',handles.reader.NumberOfFrames));
        %     set(handles.sampling_edit,'String',num2str(handles.param.video_sampling));
        set(handles.usolver.restart_edit,'Enable','on');
        set(handles.usolver.restart_edit,'Value',~(handles.param.restart));
        set(handles.usolver.do_pgd_edit,'Enable','on');
        set(handles.usolver.do_pgd_edit,'Value',handles.param.do_pgd_prediction+1);
    else
        %    set(handles.menu_video,'Visible','off')
        if iscell(handles.param.deformed_image)
            set(handles.usolver.restart_edit,'Enable','on');
            set(handles.usolver.restart_edit,'Value',~(handles.param.restart));
            set(handles.usolver.do_pgd_edit,'Enable','on');
            set(handles.usolver.do_pgd_edit,'Value',handles.param.do_pgd_prediction+1);
        else
            set(handles.usolver.restart_edit,'Enable','off');
            set(handles.usolver.do_pgd_edit,'Enable','off');
        end
        
    end
    switch handles.fbasis
        case 'fem'
            set(handles.usolver.nscale_edit,'Enable','on');
            set(handles.usolver.nscale_edit,'String',num2str(handles.nscale));
            set(handles.usolver.psample_edit,'Enable','on');
            set(handles.usolver.psample_edit,'String',num2str(handles.param.psample));
            set(handles.usolver.normalize_grey_levels_edit,'Enable','on')
            
        case 'uni'
            set(handles.usolver.nscale_edit,'Enable','on');
            set(handles.usolver.nscale_edit,'String',num2str(handles.nscale));
            set(handles.usolver.normalize_grey_levels_edit,'Enable','off')
            set(handles.usolver.psample_edit,'Enable','off');
            set(handles.usolver.psample_edit,'String',num2str(1));
            handles.param.normalize_grey_level=0;
        case 'beam'
            handles.param.normalize_grey_level=0;
            set(handles.usolver.normalize_grey_levels_edit,'Enable','off')
            set(handles.usolver.nscale_edit,'String',num2str(1));
            set(handles.usolver.nscale_edit,'Enable','off');
            set(handles.usolver.psample_edit,'Enable','off');
            set(handles.usolver.psample_edit,'String',num2str(1));
        case 'vic'
            handles.param.normalize_grey_level=0;
            set(handles.usolver.normalize_grey_levels_edit,'Enable','off')
            set(handles.usolver.nscale_edit,'Enable','on');
            set(handles.usolver.nscale_edit,'String',num2str(handles.nscale));
            set(handles.usolver.psample_edit,'Enable','off');
            set(handles.usolver.psample_edit,'String',num2str(1));
    end
    set(handles.usolver.normalize_grey_levels_edit,'Value',handles.param.normalize_grey_level);
catch
end
set(0,'CurrentFigure',handles.figure1);
[az,el] = view;
dview=((az==0)&&(el==90));
set(handles.rotate_button,'Enable','on');
set(handles.switch_view_button,'Enable','off');
set(handles.switch_frame_button,'Enable','off');
set(handles.cmenu_field_data_disp,'Visible','on')
set(handles.cmenu_field_data_strain,'Visible','on')
set(handles.cmenu_field_vref,'Visible','off')
set(handles.cmenu_field_data_axialstrain,'Visible','off')
set(handles.cmenu_field_data_shearstrain,'Visible','off')
set(handles.cmenu_mesh_factors,'Visible','off')
set(handles.cmenu_mesh_colorscale,'Visible','off')

set(handles.cmenu_field_data_temp,'Visible','off')
set(handles.cmenu_field_data_calculator,'Visible','off')
set(handles.cmenu_field_data_error,'Visible','off')
set(handles.cmenu_field_data_erroronelt,'Visible','off')
set(handles.cmenu_field_data_error,'Label','Correlation error')
set(handles.first_frame_button,'Enable','off')
set(handles.pprevious_frame_button,'Enable','off')
set(handles.previous_frame_button,'Enable','off')
set(handles.play_button,'Enable','off')
set(handles.nnext_frame_button,'Enable','off')
set(handles.next_frame_button,'Enable','off')
set(handles.last_frame_button,'Enable','off')
set(handles.cmenu_bcs_load,'Visible','off')
set(handles.cmenu_animation_images,'Visible','off');
set(handles.cmenu_animation_export_field,'Visible','off')
set(handles.cmenu_animation_field_vtk,'Visible','off')
set(handles.cmenu_mesh_show_edge,'Visible','off')
set(handles.cmenu_animation_goto,'Visible','off');
set(handles.cmenu_animation_images_vtk,'Visible','off');

set(handles.preview_button,'Enable','off');
set(handles.show_data_button,'Enable','off');

try set(handles.gim,'uicontextmenu',handles.cmenu_animation);catch , end

if handles.ana>0,set(handles.axis,'Visible','on');end

set(handles.cmenu_crack_detection,'Visible','off')
set(handles.cmenu_crack_propagation,'Visible','off')
switch  handles.ana
    case 1
        handles.showim=1;
        set(handles.figure1,'uicontextmenu',handles.cmenu_param)
        set(handles.cmenu_crack_detection,'Visible','on')
        set(handles.cmenu_mesh_zone,'Visible','on')
        set(handles.cmenu_field_data_strain,'Visible','off')
        set(handles.save_button,'Enable','on')
        set(handles.cmenu_field_showim,'Visible','off')
        set(handles.first_frame_button,'Enable','on')
        set(handles.previous_frame_button,'Enable','on')
        set(handles.pprevious_frame_button,'Enable','on')
        set(handles.play_button,'Enable','on')
        set(handles.next_frame_button,'Enable','on')
        set(handles.nnext_frame_button,'Enable','on')
        set(handles.last_frame_button,'Enable','on')
        if (~(handles.stereo))
            if ~(strcmp(handles.fbasis,'beam'))&&~(strcmp(handles.fbasis,'vic'))
                %                set(handles.preview_button,'Enable','on');
            end
        end
        set(handles.cmenu_animation_images,'Visible','on');
    case 2
        switch handles.param.analysis
            case 'correlation'
                
                set(handles.cmenu_animation_field_vtk,'Visible','on')
                set(handles.cmenu_animation_images_vtk,'Visible','on');
                if ~(handles.param.detect>0)&&~(handles.param.ulive==1)
                    if handles.erroronelt
                        set(handles.cmenu_field_data_erroronelt,'Visible','on')
                    else
                        set(handles.cmenu_field_data_error,'Visible','on')
                    end
                else
                    if handles.cracked==2
                        set(handles.cmenu_field_data_error,'Visible','on')
                        set(handles.cmenu_field_data_error,'Label','Fit error')
                    end
                end
                set(handles.figure1,'uicontextmenu','')
                set(handles.cmenu_mesh_zone,'Visible','off')
                if strcmp(handles.fbasis,'uni')
                    set(handles.cmenu_field_data_strain,'Visible','off')
                    set(handles.save_button,'Enable','off')
                end
                if handles.param.thermo==1
                    set(handles.cmenu_field_data_temp,'Visible','on')
                end
                switch handles.fbasis
                    case 'uni'
                        set(handles.cmenu_mesh_colorscale,'Visible','on')
                    case 'vic'
                        set(handles.cmenu_mesh_colorscale,'Visible','on')
                        set(handles.cmenu_mesh_factors,'Visible','on')
                        set(handles.cmenu_field_data_strain,'Visible','off')
                        set(handles.cmenu_field_finite_strain,'Visible','off')
                        set(handles.cmenu_field_vref,'Visible','on')
                    case 'beam'
                        set(handles.cmenu_field_data_strain,'Visible','off')
                        set(handles.cmenu_field_finite_strain,'Visible','off')
                        set(handles.cmenu_field_data_axialstrain,'Visible','on')
                        if strcmp(handles.beam_model.beam_type,'timoshenko')
                            set(handles.cmenu_field_data_shearstrain,'Visible','on')
                        end
                        set(handles.cmenu_mesh_factors,'Visible','on')
                        set(handles.cmenu_mesh_colorscale,'Visible','on')
                        set(handles.cmenu_animation_export_field,'Visible','on')
                    case 'fem'
                        set(handles.cmenu_mesh_colorscale,'Visible','on')
                        set(handles.cmenu_mesh_factors,'Visible','on')
                        set(handles.rotate_button,'Enable','on');
                        set(handles.switch_view_button,'Enable','on');
                        if handles.stereo,set(handles.switch_frame_button,'Enable','on');end
                        set(handles.cmenu_field_data_calculator,'Visible','on')
                        set(handles.cmenu_animation_export_field,'Visible','on')
                end
                if handles.preview
                    set(handles.cmenu_field_showim,'Visible','on')
                    set(handles.cmenu_field_ondefimage,'Visible','on')
                else
                    set(handles.cmenu_field_showim,'Visible','off')
                    set(handles.cmenu_field_ondefimage,'Visible','off')
                end
            case 'mechanics'
                
                handles.showim=0;
                set(handles.figure1,'uicontextmenu','')
                set(handles.cmenu_mesh_zone,'Visible','off')
                set(handles.cmenu_mesh_colorscale,'Visible','on')
                set(handles.cmenu_mesh_factors,'Visible','on')
                set(handles.rotate_button,'Enable','on');
                set(handles.switch_view_button,'Enable','on');
                set(handles.cmenu_field_data_calculator,'Visible','on')
                set(handles.cmenu_field_ondefimage,'Visible','on')
                
                set(handles.save_button,'Enable','on')
                set(handles.cmenu_field_showim,'Visible','off')
                set(handles.cmenu_bcs_load,'Visible','off')
        end
        set(handles.first_frame_button,'Enable','on')
        set(handles.previous_frame_button,'Enable','on')
        set(handles.pprevious_frame_button,'Enable','on')
        set(handles.play_button,'Enable','on')
        set(handles.next_frame_button,'Enable','on')
        set(handles.nnext_frame_button,'Enable','on')
        set(handles.last_frame_button,'Enable','on')
        set(handles.preview_button,'Enable','on');
        set(handles.show_data_button,'Enable','on');
        
    case 3
        handles.showim=0;
        set(handles.cmenu_mesh_zone,'Visible','on')
        set(handles.save_button,'Enable','on')
        set(handles.cmenu_field_showim,'Visible','off')
        set(handles.cmenu_bcs_load,'Visible','on')
        
        set(handles.cmenu_crack_propagation,'Visible','on')
        
end
set(handles.cmenu_field_ondefimage,'Enable','on')
if handles.field==3
    if handles.ercomp==0&&handles.cracked<2
        handles.ondefimage=0;
        set(handles.cmenu_field_ondefimage,'Enable','off')
    end
end
if handles.showedge
    set(handles.cmenu_mesh_show_edge,'Checked','on')
else
    set(handles.cmenu_mesh_show_edge,'Checked','off')
end
if handles.showcb
    set(handles.cmenu_field_colorbar,'Checked','on')
else
    set(handles.cmenu_field_colorbar,'Checked','off')
end
if handles.ondefimage
    set(handles.cmenu_field_ondefimage,'Checked','on')
else
    set(handles.cmenu_field_ondefimage,'Checked','off')
end
if handles.rmrbm
    set(handles.cmenu_field_rmrbm,'Checked','on')
else
    set(handles.cmenu_field_rmrbm,'Checked','off')
end
if handles.fstrain
    set(handles.cmenu_field_finite_strain,'Checked','on')
else
    set(handles.cmenu_field_finite_strain,'Checked','off')
end
if handles.vref
    set(handles.cmenu_field_vref,'Checked','on')
else
    set(handles.cmenu_field_vref,'Checked','off')
end
if handles.showim
    set(handles.cmenu_field_showim,'Check','on')
    try set(handles.gim,'Visible','on'); catch, end
else
    set(handles.cmenu_field_showim,'Check','off')
    set(handles.gim,'Visible','off')
    set(handles.gim,'Visible','off')
end
try set(handles.gpoint,'Visible','off');catch, end
try set(handles.gmesh,'Visible','off');catch, end
try set(handles.gedge,'Visible','off');catch, end
for iz=1:size(handles.fem_model.zone,2)
    try
        set((handles.fem_model.zone{3,iz}),'Visible','off');
    catch
    end
end
for iz=1:size(handles.fem_model.zone,2)
    try
        set((handles.fem_model.zone{11,iz}),'Visible','off');
    catch
    end
end
for iz=1:size(handles.uni_model.zone,2)
    try
        set((handles.uni_model.zone{3,iz}),'Visible','off');
    catch
    end
end
for iz=1:size(handles.beam_model.zone,2)
    try
        set((handles.beam_model.zone{3,iz}),'Visible','off');
    catch
    end
end
for iz=1:size(handles.vic_model.zone,2)
    try
        set((handles.vic_model.zone{3,iz}),'Visible','off');
    catch
    end
end
for iz=1:size(handles.erzone,2)
    try
        set((handles.erzone{3,iz}),'Visible','off');
    catch
    end
end
try set(handles.gcbar,'uicontextmenu',handles.cmenu_colorbar);catch, end
if handles.preview&&handles.ana
    set(handles.field_text,'Visible','on');
else
    set(handles.field_text,'Visible','off');
end

switch handles.ana
    case 1
        switch handles.fbasis
            case 'fem'
                if handles.animation.iim>0&&handles.preview
                    set(handles.gmesh,'Visible','on');
                end
                if handles.animation.iim==0||(handles.preview&&handles.showedge)
                    set(handles.gedge,'Visible','on');
                end
                if handles.animation.iim==0||(handles.preview&&~(handles.ondefimage))
                    for iz=1:size(handles.fem_model.zone,2)
                        try
                            set((handles.fem_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                end
                set(handles.cmenu_mesh_show_edge,'Visible','off')
            case 'uni'
                if handles.animation.iim==0||(handles.preview&&~(handles.ondefimage))
                    for iz=1:size(handles.uni_model.zone,2)
                        try
                            set((handles.uni_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                end
                if handles.preview&&handles.animation.iim>0
                    set(handles.gmesh,'Visible','on')
                end
            case 'beam'
                if handles.animation.iim==0||(handles.preview&&~(handles.ondefimage))
                    for iz=1:size(handles.beam_model.zone,2)
                        try
                            set((handles.beam_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                end
                if handles.preview&&handles.animation.iim>0
                    set(handles.gmesh,'Visible','on')
                end
            case 'vic'
                if handles.animation.iim==0||(handles.preview&&~(handles.ondefimage))
                    for iz=1:size(handles.vic_model.zone,2)
                        try
                            set((handles.vic_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                end
                if handles.preview&&handles.animation.iim>0
                    set(handles.gmesh,'Visible','on')
                end
                
            case 'none'
                if handles.preview&&handles.animation.iim>0
                    set(handles.gmesh,'Visible','on');
                end
        end
    case 2
        switch handles.fbasis
            case 'fem'
                if handles.animation.iim>0&&handles.preview
                    try set(handles.gmesh,'Visible','on');catch, end
                end
                if handles.animation.iim==0||(handles.preview&&handles.showedge)
                    try set(handles.gedge,'Visible','on');catch, end
                end
                if ~isfield(handles.fem_model,'degree')
                    set(handles.cmenu_mesh_show_edge,'Visible','on')
                end
                if handles.preview&&dview
                    for iz=1:size(handles.fem_model.zone,2)
                        try
                            set((handles.fem_model.zone{11,iz}),'Visible','on');
                        catch
                        end
                    end
                    
                end
                if (handles.animation.iim==0||(handles.preview&&~(handles.ondefimage)))&&dview
                    for iz=1:size(handles.fem_model.zone,2)
                        try
                            set((handles.fem_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                    try set(handles.gpoint,'Visible','on');catch, end
                    if handles.field==3&&handles.ercomp==0&&(handles.erroronelt||isfield(handles.fem_model,'degree'))
                        for iz=1:size(handles.erzone,2)
                            try
                                set((handles.erzone{3,iz}),'Visible','on');
                            catch
                            end
                        end
                    end
                end
                if handles.field==4
                    if handles.animation.iim
                        set(handles.gmesh,'Visible','off');
                    else
                        set(handles.gmesh,'Visible','on');
                    end
                end
            case 'uni'
                if handles.animation.iim==0||(handles.preview&&~(handles.ondefimage))
                    for iz=1:size(handles.uni_model.zone,2)
                        try
                            set((handles.uni_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                end
                if handles.preview&&handles.animation.iim>0+handles.param.ulive
                    set(handles.gmesh,'Visible','on')
                end
            case 'beam'
                if handles.animation.iim==0||(handles.preview&&~(handles.ondefimage))
                    for iz=1:size(handles.beam_model.zone,2)
                        try
                            set((handles.beam_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                end
                if handles.preview&&handles.animation.iim>0
                    set(handles.gmesh,'Visible','on')
                end
            case 'vic'
                if handles.animation.iim==0||(handles.preview&&~(handles.ondefimage))
                    for iz=1:size(handles.vic_model.zone,2)
                        try
                            set((handles.vic_model.zone{3,iz}),'Visible','on');
                        catch
                        end
                    end
                end
                if handles.preview&&handles.animation.iim>0
                    set(handles.gmesh,'Visible','on')
                end
        end
    case 3
        set(handles.gedge,'Visible','on');
        for iz=1:size(handles.fem_model.zone,2)
            try
                set((handles.fem_model.zone{3,iz}),'Visible','on');
            catch
            end
        end
        
end
if ~isfield(handles.fem_model,'mesh_type')
    handles.fem_model.mesh_type=1;
end
set_fem_param(handles)
set_uni_param(handles)
set_beam_param(handles)
set(handles.cmenu_param_model,'Visible','off')
set(handles.cmenu_param_vic,'Visible','off')
set(handles.cmenu_param_model_uni,'Visible','off')
set(handles.cmenu_param_model_beam,'Visible','off')
set(handles.cmenu_param_model_fem,'Visible','off')
set(handles.cmenu_param_model_nurbs,'Visible','off')
set(handles.new_fem_analysis,'Enable','off')
set(handles.new_thermo_analysis,'Enable','off')
set(handles.cmenu_mesh_crack,'Enable','off')
set(handles.cmenu_mesh_crack,'Visible','on')
set(handles.new_vic_analysis,'Enable','off');
set(handles.new_3d_analysis,'Enable','off');

if any(handles.logiciel==1),set(handles.cmenu_param_model_fem,'Visible','on');end
if any(handles.logiciel==2)&&(~handles.stereo),set(handles.cmenu_param_model_uni,'Visible','on');end
if any(handles.logiciel==3)&&(~handles.stereo),set(handles.cmenu_param_model_beam,'Visible','on');end
if any(handles.logiciel==4),set(handles.new_fem_analysis,'Enable','on');end
if any(handles.logiciel==5)
    set(handles.cmenu_mesh_crack,'Enable','on');
    if (isfield(handles.fem_model,'degree')||(handles.param.ulive==1))||~strcmp(handles.fbasis,'fem')
        set(handles.cmenu_mesh_crack,'Enable','off');
    end
    if (handles.ana==2)&&(iscell(handles.uvisu)||strcmp(handles.param.analysis,'mechanics'))
        set(handles.cmenu_mesh_crack,'Visible','off');
    end
end
if any(handles.logiciel==7), set(handles.new_thermo_analysis,'Enable','on');end
if any(handles.logiciel==8)
    set(handles.new_vic_analysis,'Enable','on');
    set(handles.cmenu_param_vic_point,'Visible','off');
end
if any(handles.logiciel==9)
    set(handles.new_3d_analysis,'Enable','on');
end
if any(handles.logiciel==6),set(handles.cmenu_param_model_nurbs,'Visible','on');end
if handles.stereo
    set(handles.cmenu_disp_data_uz,'Visible','on')
    set(handles.cmenu_field_data_topo,'Visible','on')
else
    set(handles.cmenu_field_data_topo,'Visible','off')
    set(handles.cmenu_disp_data_uz,'Visible','off')
end
if handles.param.vic==1
    set(handles.cmenu_param_vic,'Visible','on')
else
    set(handles.cmenu_param_model,'Visible','on')
end
guidata(handles.figure1,handles);
if ishghandle(handles.fparams)&&~isreal(handles.fparams)
    double_clic(handles.figure1,{1})
end
handles=guidata(handles.figure1);
set(0,'CurrentFigure',handles.figure1)

function set_fem_param(handles)
set(handles.cmenu_mesh_mesh_size_text,'Label',num2str(handles.fem_model.mesh_size(1)))
set(handles.cmenu_mesh_degree,'Visible','off');
set(handles.cmenu_mesh_elt_type_tri,'Checked','off')
set(handles.cmenu_mesh_elt_type_quad,'Checked','off')
set(handles.cmenu_mesh_elt_type_adapt,'Checked','off')
set(handles.cmenu_mesh_load_mesh_load,'Checked','off')
set(handles.cmenu_zone_attractor,'Visible','off')
set(handles.cmenu_mesh_attractors,'Visible','off')
set(handles.cmenu_mesh_bcs,'Visible','off')
set(handles.cmenu_mesh_material,'Visible','off')
switch handles.fem_model.mesh_type
    case 1
        switch handles.fem_model.element_type
            case 3
                set(handles.cmenu_mesh_elt_type_tri,'Checked','on')
            case 4
                set(handles.cmenu_mesh_elt_type_quad,'Checked','on')
        end
    case 2
        set(handles.cmenu_mesh_elt_type_adapt,'Checked','on')
        set(handles.cmenu_zone_attractor,'Visible','on')
        set(handles.cmenu_mesh_attractors,'Visible','on')
        
    case 3
        set(handles.cmenu_mesh_load_mesh_load,'Checked','on')
        
end
set(handles.cmenu_mesh_smoothing_median,'Checked','off')
set(handles.cmenu_mesh_smoothing_strain,'Checked','off')
set(handles.cmenu_mesh_smoothing_stress,'Checked','off')
set(handles.cmenu_mesh_smoothing_lc,'Label',num2str(handles.param.regularization_parameter))

if handles.param.regularization_parameter
    switch handles.param.regularization_type
        case 'median'
            set(handles.cmenu_mesh_smoothing_median,'Checked','on')
        case 'tiko'
            set(handles.cmenu_mesh_smoothing_strain,'Checked','on')
        case 'equilibrium_gap'
            set(handles.cmenu_mesh_smoothing_stress,'Checked','on')
    end
end
set(handles.cmenu_mesh_load_mesh,'Visible','off')
set(handles.cmenu_mesh_save_mesh,'Visible','on')

switch handles.ana
    case 1
        set(handles.cmenu_mesh_create_gage,'Visible','off')
        set(handles.cmenu_mesh_plot_line,'Visible','off')
        set(handles.cmenu_mesh_save_inp,'Visible','off')
        set(handles.cmenu_mesh_elt_type,'Visible','on')
        set(handles.cmenu_mesh_mesh_size,'Visible','on')
        set(handles.cmenu_mesh_smoothing,'Visible','on')
        set(handles.cmenu_mesh_attractors,'Visible','on')
        set(handles.cmenu_mesh_entropy,'Visible','on')
        if isfield(handles.fem_model,'material_parameters')
            set(handles.cmenu_mesh_material,'Visible','on')
        end
        if isfield(handles.fem_model,'degree')
            set(handles.cmenu_mesh_degree,'Visible','on');
        end
    case 2
        switch handles.fbasis
            case 'fem'
                if ~iscell(handles.uvisu)&&(handles.param.ulive==0)
                    set(handles.cmenu_mesh_create_gage,'Visible','on')
                    set(handles.cmenu_mesh_plot_line,'Visible','on')
                    set(handles.cmenu_mesh_save_inp,'Visible','on')
                else
                    set(handles.cmenu_mesh_create_gage,'Visible','off')
                    set(handles.cmenu_mesh_plot_line,'Visible','off')
                    set(handles.cmenu_mesh_save_inp,'Visible','off')
                end
                set(handles.cmenu_mesh_elt_type,'Visible','off')
                set(handles.cmenu_mesh_mesh_size,'Visible','off')
                set(handles.cmenu_mesh_smoothing,'Visible','off')
                set(handles.cmenu_mesh_attractors,'Visible','off')
                set(handles.cmenu_mesh_entropy,'Visible','off')
                set(handles.cmenu_mesh_degree,'Visible','off');
            case {'vic','beam','uni'}
                set(handles.cmenu_mesh_create_gage,'Visible','off')
                set(handles.cmenu_mesh_plot_line,'Visible','off')
                set(handles.cmenu_mesh_save_inp,'Visible','off')
                set(handles.cmenu_mesh_elt_type,'Visible','off')
                set(handles.cmenu_mesh_mesh_size,'Visible','off')
                set(handles.cmenu_mesh_smoothing,'Visible','off')
                set(handles.cmenu_mesh_attractors,'Visible','off')
                set(handles.cmenu_mesh_entropy,'Visible','off')
                set(handles.cmenu_mesh_degree,'Visible','off');
                set(handles.cmenu_mesh_save_mesh,'Visible','off')
                
        end
    case 3
        set(handles.cmenu_mesh_create_gage,'Visible','off')
        set(handles.cmenu_mesh_plot_line,'Visible','off')
        set(handles.cmenu_mesh_save_inp,'Visible','off')
        set(handles.cmenu_mesh_elt_type,'Visible','off')
        set(handles.cmenu_mesh_mesh_size,'Visible','on')
        set(handles.cmenu_mesh_smoothing,'Visible','off')
        set(handles.cmenu_mesh_attractors,'Visible','on')
        set(handles.cmenu_mesh_entropy,'Visible','off')
        set(handles.cmenu_mesh_bcs,'Visible','on')
        set(handles.cmenu_mesh_material,'Visible','on')
        set(handles.cmenu_mesh_degree,'Visible','off');
end
if isfield(handles.fem_model,'degree')
    set(handles.cmenu_mesh_degree_text,'Label',num2str(handles.fem_model.degree(1)))
end
function set_uni_param(handles)
%set(handles.gage_id_edit,'String','')
%gage_id_edit_Callback(handles.gage_id_edit,[],handles)
%if handles.preview
switch handles.ana
    case 1
        set(handles.cmenu_gage_add_gage,'Visible','on')
        set(handles.cmenu_gage_remove,'Visible','on')
        if handles.preview
            set(handles.cmenu_gage_plot,'Visible','on')
        else
            set(handles.cmenu_gage_plot,'Visible','off')
        end
    case 2
        set(handles.cmenu_gage_remove,'Visible','off')
        set(handles.cmenu_gage_add_gage,'Visible','off')
end


function set_beam_param(handles)
if handles.beam_model.exx
    set(handles.cmenu_beam_exx,'Checked','on')
else
    set(handles.cmenu_beam_exx,'Checked','off')
end

set(handles.cmenu_beam_deg_2,'Checked','off')
set(handles.cmenu_beam_deg_3,'Checked','off')
set(handles.cmenu_beam_deg_4,'Checked','off')
set(handles.cmenu_beam_deg_5,'Checked','off')
set(handles.cmenu_beam_type_beam,'Checked','off')
set(handles.cmenu_beam_type_beam_nurbs,'Checked','off')
set(handles.cmenu_beam_kine_euler,'Checked','off')
set(handles.cmenu_beam_kine_timosh,'Checked','off')

switch handles.beam_model.basis
    case 'beam'
        set(handles.cmenu_beam_type_beam,'Checked','on')
        set(handles.cmenu_beam_deg_2,'Checked','on')
        set(handles.cmenu_beam_nbelt_edit,'Label','1')
        set(handles.cmenu_beam_deg,'Enable','off')
        set(handles.cmenu_beam_nbelt,'Enable','off')
        set(handles.cmenu_beam_kine,'Enable','off')
    case 'beam-nurbs'
        set(handles.cmenu_beam_nbelt_edit,'Label',num2str(handles.beam_model.nb_element))
        set(handles.cmenu_beam_type_beam_nurbs,'Checked','on')
        switch handles.beam_model.degree
            case 2
                set(handles.cmenu_beam_deg_2,'Checked','on')
            case 3
                set(handles.cmenu_beam_deg_3,'Checked','on')
            case 4
                set(handles.cmenu_beam_deg_4,'Checked','on')
            case 5
                set(handles.cmenu_beam_deg_5,'Checked','on')
        end
        switch handles.beam_model.beam_type
            case 'euler'
                set(handles.cmenu_beam_kine_euler,'Checked','on')
            case 'timoshenko'
                set(handles.cmenu_beam_kine_timosh,'Checked','on')
        end
        set(handles.cmenu_beam_deg,'Enable','on')
        set(handles.cmenu_beam_nbelt,'Enable','on')
        set(handles.cmenu_beam_kine,'Enable','on')
end

function set_vic_param(handles,iz)
zone=handles.vic_model.zone(:,iz);
if zone{1}>0
    set(handles.cmenu_vic_type_solid,'Checked','off')
    set(handles.cmenu_vic_type_wire,'Checked','on')
    set(handles.cmenu_vic_t,'Visible','on')
else
    set(handles.cmenu_vic_type_solid,'Checked','on')
    set(handles.cmenu_vic_type_wire,'Checked','off')
    set(handles.cmenu_vic_t,'Visible','off')
end

set(handles.cmenu_vic_nbelt_edit,'Label',num2str(zone{6}))
set(handles.cmenu_vic_deg_edit,'Label',num2str(zone{7}))
set(handles.cmenu_vic_ep_edit,'Label',num2str(zone{8}))
set(handles.cmenu_vic_t_edit,'Label',num2str(zone{9}))


% --- Outputs from this function are returned to the command line.
function varargout=Umaster_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=handles.output;
% --------------------------------------------------------------------
function new_thermo_analysis_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to new_thermo_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filc,pathc,FilterIndexc]=uigetfile({...
    '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2;*.raw','Image files';...
    '*.cal','Calibration files';...
    '*.m2v;*.avi;*.AVI;*.mov;*.MOV','Video files'},'Open calibration file(s)','MultiSelect', 'on');
if FilterIndexc
    cd(pathc);
    [fil0,path0,FilterIndex]=uigetfile({...
        '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2;*.raw','Image files';...
        '*.m2v;*.avi;*.AVI;*.mov;*.MOV','Video files'},'Open file(s)','MultiSelect', 'on');
    if FilterIndex
        handles=reset_param(handles);
        switch FilterIndexc
            case {1,3}
                filc=sort(filc);
                fo=0;
                for iim=(length(filc)-9):length(filc)
                    tmp=readim(filc{iim});
                    fo=fo+tmp;
                end
                fo=double(fo)/10;
                if max(fo(:))>255,fo=uint16(fo);else,fo=uint8(fo);end
                handles.param.ir_calibration_data={pathc,filc};
            case 2
                load(filc,'-mat','T');
                fo=T.fo;
                handles.param.ir_calibration_data=T;
        end
        cd(path0);
        [pp,ff,ext]=fileparts(fil0{1});
        writeim(['00',ext],fo);
        fil0=[{['00',ext]},fil0];
        handles.param.thermo=1;
        switch FilterIndex
            case {1,2}
                start_new_analysis(hObject,handles,fil0,path0,FilterIndex);
        end
        
        
        %    set(handles.animation_panel,'Visible','on')
        %    set(handles.end_zone,'Visible','on')
        %set(handles.begin_zone,'Visible','on')
        set(handles.image_info_text,'Visible','on')
        
    end
end
% --------------------------------------------------------------------
function new_vic_analysis_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to new_vic_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fil0,path0,FilterIndex]=uigetfile({...
    '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2;*.raw','Image files';...
    '*.m2v;*.avi;*.AVI;*.mov;*.MOV','Video files'},'Open file(s)','MultiSelect', 'on');
if FilterIndex
    handles=reset_param(handles);
    cd(path0);
    handles.fbasis='vic';
    handles.param.vic=1;
    handles=basis_functions(handles);
    switch FilterIndex
        case {1,2}
            start_new_analysis(hObject,handles,fil0,path0,FilterIndex);
    end
    
    set(handles.image_info_text,'Visible','on')
    
end

% --------------------------------------------------------------------
function new_dic_analysis_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to new_dic_analysis_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fil0,path0,FilterIndex]=uigetfile({...
    '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2;*.raw','Image files';...
    '*.m2v;*.avi;*.AVI;*.mov;*.MOV','Video files'%;...
    %'*.mat;*.tiff;*.raw;*.bin','3D image files'
    },'Open file(s)','MultiSelect', 'on');
if FilterIndex
    handles=reset_param(handles);
    cd(path0);
    switch FilterIndex
        case {1,2}
            start_new_analysis(hObject,handles,fil0,path0,FilterIndex);
    end
    
    
    %    set(handles.animation_panel,'Visible','on')
    %    set(handles.end_zone,'Visible','on')
    %set(handles.begin_zone,'Visible','on')
    set(handles.image_info_text,'Visible','on')
    
end

% --------------------------------------------------------------------
function open_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to open_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fil0,path0,FilterIndex]=uigetfile({...
    '*.res;','Ufreckles result files (*.res)';...
    '*.dat;','Ufreckles input files (*.dat)';...
    '*.ufr;','Ufreckles input ascii files (*.ufr)'...
    },'Open file(s)','MultiSelect', 'on');
if FilterIndex
    handles=reset_param(handles);
    
    set(handles.figure1,'Name',sprintf('UFreckles: %s',fil0))
    cd(path0);
    switch FilterIndex
        case {1,2}
            loading_mat_file(hObject,handles,fil0,path0,FilterIndex);
        case 3
            loading_ufr_file(hObject,handles,fil0,path0,FilterIndex);
            
    end
    
    
    %    set(handles.animation_panel,'Visible','on')
    %    set(handles.end_zone,'Visible','on')
    %set(handles.begin_zone,'Visible','on')
    set(handles.image_info_text,'Visible','on')
    
end
function loading_mat_file(hObject,handles,fil0,path0,FilterIndex)
couls='rbmycg';
load(fil0,'-mat','model','param','U','cn');
if ~isfield(param,'vic'), param.vic=0;end
if ~isfield(param,'ulive'), param.ulive=0;end
if ~isfield(param,'detect'), param.detect=0;end
if ~isfield(param,'thermo'), param.thermo=0;end
if ~isfield(param,'psample'), param.psample=1;end
handles.param=param;
switch handles.param.analysis
    case 'correlation'
        handles.ana=1;
        if isfield(param,'deformed_image')
            nim=1;
            filref=param.reference_image;
            try im00=(readim(filref));
            catch, im00=[];end
            if iscell(param.deformed_image)
                nim=size(param.deformed_image,2);
                handles.stereo=size(param.deformed_image,1)==2;
            end
            handles.animation.nbstep=nim-handles.stereo;
            handles.animation.frames=1:(nim-handles.stereo);
        else
            filref=param.reference_image;
            reader=VideoReader(param.reference_image);
            handles.reader=reader;
            im00=readim(reader,1);
            handles.animation.frames=2:handles.param.video_sampling:handles.param.number_of_frames;
            handles.animation.nbstep=length(handles.animation.frames);
        end
        if size(im00,3)==1%Cas d'une image NB
            im00=repmat(im00,[1,1,3]);
            rgb='MN';
        else
            rgb='RGB';
        end
        handles.sizeim=size(im00);
        img_msg=sprintf('%s %dx%d %d-bit %s',filref,handles.sizeim(1:2),8*(1+double(max(im00(:))>255)),rgb);
        set(handles.image_info_text,'String',img_msg);
    case 'mechanics'
        set(handles.image_info_text,'String','');
        
        handles.ana=3;
        sizeim=handles.param.roi([2,4])+handles.param.roi([1,3])-1;
        im00=ones([sizeim,3]);
        handles.sizeim=size(im00);
end
if any(im00(:)>255)
    im00=double(im00);
    im00=255*(im00-min(im00(:)))/(max(im00(:))-min(im00(:)));
    %     if any(im00(:)>1023)
    %         if any(im00(:)>4095)
    %             im00=im00*(255/65535);
    %         else
    %             im00=im00*(255/4095);
    %         end
    %     else
    %         im00=im00*(255/1023);
    %     end
    im00=uint8(floor(im00));
end
set(0,'CurrentFigure',handles.figure1);
gim=imagesc(permute(im00,[2,1,3]));
handles.gim=gim;
hold on;
axis equal
axis xy;
axis off
set(handles.axes1,'Clim',[0,255]);
set(handles.axes1,'ClimMode','manual');

colormap(handles.cmap);

handles.gcbar=colorbar('location','East','FontSize',20,'FontWeight','normal','LineWidth',1,'uicontextmenu',handles.cmenu_colorbar);
%set(handles.gcbar,'Visible','off','Ycolor',[0,0,1]*0.8);
set(handles.gcbar,'Visible','off','Ycolor',[0.7964    0.2780    0.4713]);
%set(handles.gcbar,'Visible','off','Ycolor',[0.9883    0.6523    0.2114]);
%set(handles.gcbar,'ButtonDownFcn','selectmoveresize');

switch model.basis
    case 'fem'
        if ~isfield(model,'phantom_nodes'), model.phantom_nodes=0;end
        if ~isfield(model,'zone'), model.zone={};end
        handles.fem_model=model;
        handles.fbasis='fem';
    case {'beam','beam-nurbs'}
        
        if ~isfield(model,'beam_type');
            model.beam_type='euler';
        end
        handles.beam_model=model;
        handles.fbasis='beam';
        gage=handles.beam_model.zone(:,1);
        xyp=gage{2};
        gz=plot(xyp(:,1),xyp(:,2),'Color',[0.75,0,0],'LineStyle','-','LineWidth',2,'uicontextmenu',handles.cmenu_export_beam);
        gz.Tag=num2str(rand(1));
        handles.beam_model.zone{3,1}=gz;
        
    case 'uni'
        handles.uni_model=model;
        handles.fbasis='uni';
        for ig=1:size(handles.uni_model.zone,2)
            gage=handles.uni_model.zone(:,ig);
            xyp=gage{2};
            gz=plot(xyp(:,1),xyp(:,2),'Color',couls(ig),'LineStyle','-','LineWidth',2,'uicontextmenu',handles.cmenu_gage);
            gz.Tag=num2str(rand(1));
            handles.uni_model.zone{3,ig}=gz;
            ax1=annotation('line',0.006+[0,1/40],0.95+[0,0]-(ig-1)*0.025);
            set(ax1,'Color',couls(ig),'LineWidth',2);
            tx1=annotation('textbox',[0.007+1/40,0.95-(ig-1)*0.025,0.02,0.02],'String',sprintf('Gage #%d',ig'));
            set(tx1,'Color',couls(ig),'LineStyle','none','FitBoxToText','on','FontSize',12);
            handles.ganot{ig}=[ax1,tx1];
        end
    case 'vic'
        handles.vic_model=model;
        handles.fbasis='vic';
        for ig=1:size(handles.vic_model.zone,2)
            gage=handles.vic_model.zone(:,ig);
            xyp=gage{2};
            gz=plot(xyp(:,1),xyp(:,2),'Color',[1,0.5,0],'LineStyle','-','LineWidth',2);
            gz.Tag=num2str(rand(1));
            if exist('U','var')
                set(gz,'uicontextmenu',handles.cmenu_export_vic);
            else
                set(gz,'uicontextmenu',handles.cmenu_vic,'ButtonDownFcn',@zone_adjust);
            end
            handles.vic_model.zone{3,ig}=gz;
        end
        
end
if FilterIndex==1%   exist('U','var')
    handles.ana=2;
    if strcmp(handles.param.analysis,'mechanics')
        handles.animation.nbstep=size(U,2);
        handles.animation.frames=1:size(U,2);
        if exist('cn','var')
            handles.fem_model.cn=cn;
        end
    end
    handles.animation.iim=handles.animation.nbstep;
    set(handles.message_text,'String','Loading result file.....')
    handles.preview=1;
    try
        if any(cell2mat(handles.fem_model.zone(4,:))==5)
            handles.cracked=1;
        end
    catch
    end
    if handles.stereo
        load(fil0,'-mat','U1','Xo','Yo','Zo');
        handles.uvisu=U1(:,2:end,1);
        handles.uxyz=U;
        handles.mvisu.Xo=Xo;
        handles.mvisu.Yo=Yo;
        handles.mvisu.Zo=Zo;
    else
        handles.uvisu=U;
    end
    load(fil0,'-mat','xo','yo','Nnodes','Nelems','elt','conn','rint');
    if iscell(xo)
        handles.mvisu.xos=xo;
        handles.mvisu.yos=yo;
        handles.mvisu.conns=conn;
        xo=xo{end};
        yo=yo{end};
        conn=conn{end};
        Nnodes=size(xo);
        Nelems=[size(conn,1),1,1];
        elt=3*ones(prod(Nelems),1);
    end
    handles.mvisu.xo=xo;
    handles.mvisu.yo=yo;
    handles.mvisu.Nnodes=Nnodes;
    handles.mvisu.Nelems=Nelems;
    handles.mvisu.conn=conn;
    handles.mvisu.elt=elt;
    handles.mvisu.smin=0;
    handles.mvisu.entropy=1;%ones(size(xo));
    val=double(conn>0);
    conn=max(conn,1);
    
    if 1
        delete(fullfile('TMP','ufreckles.res'))
        clear ComputeStrain MedianFilterCell
        %         save(fullfile('TMP','ufreckles.res'),'xo','yo','conn','elt','Nnodes','Nelems','rint');
        %         if handles.stereo
        %             save(fullfile('TMP','ufreckles.res'),'Xo','Yo','Zo','-append');
        %         end
        %         handles=ComputeStrain(handles,1);
        
    else
        if handles.stereo
            [dphidxo,dphidyo,dphidzo]=CreateGradFiniteElementBasis25D(fil0,handles.sizeim,1,[],'Gauss_points');
            [dphidx,dphidy]=CreateGradFiniteElementBasis(fil0,handles.sizeim,1,[],'Gauss_points');
            [phi]=CreateFiniteElementBasis(fil0,handles.sizeim,1,[],'Gauss_points');
            [wdetJ]=GetWeigthDetJ(fil0,handles.sizeim,1,'Gauss_points');
            M=phi'*(wdetJ*phi);
            R=dphidx'*(wdetJ*dphidx)+dphidy'*(wdetJ*dphidy);
            
            lc=(min(handles.sizeim)/10);
            lm=sqrt((max(xo)-min(xo))*(max(yo)-min(yo))/length(elt))/20;
            V=cos(2*pi*(xo)/lc).*cos(2*pi*(yo)/lc);
            a=(V'*M*V)/(V'*R*V);
            a=a*(pi*lm/lc)^2;
            M=M+a*R;
            %   handles.mvisu.dphidx=M\(phi'*(wdetJ*dphidxo));
            %   handles.mvisu.dphidy=M\(phi'*(wdetJ*dphidyo));
            %              tic
            handles.evisu.xx=M\(phi'*(wdetJ*(dphidxo*U((1:prod(Nnodes)),:))));
            handles.evisu.yy=M\(phi'*(wdetJ*(dphidyo*U(prod(Nnodes)+(1:prod(Nnodes)),:))));
            handles.evisu.xy=M\(phi'*(wdetJ*(dphidyo*U((1:prod(Nnodes)),:))));
            handles.evisu.yx=M\(phi'*(wdetJ*(dphidyo*U((1:prod(Nnodes)),:))));
        else
            [dphidx,dphidy]=CreateGradFiniteElementBasis(fil0,handles.sizeim,1,[],'Gauss_points');
            [phi]=CreateFiniteElementBasis(fil0,handles.sizeim,1,[],'Gauss_points');
            [wdetJ]=GetWeigthDetJ(fil0,handles.sizeim,1,'Gauss_points');
            M=phi'*(wdetJ*phi);
            R=dphidx'*(wdetJ*dphidx)+dphidy'*(wdetJ*dphidy);
            
            lc=(min(handles.sizeim)/10);
            lm=sqrt((max(xo)-min(xo))*(max(yo)-min(yo))/length(elt))/20;
            V=cos(2*pi*(xo)/lc).*cos(2*pi*(yo)/lc);
            a=(V'*M*V)/(V'*R*V);
            a=a*(pi*lm/lc)^2;
            M=M+a*R;
            
            handles.evisu.xx=M\(phi'*(wdetJ*(dphidx*U((1:prod(Nnodes)),:))));
            handles.evisu.yy=M\(phi'*(wdetJ*(dphidy*U(prod(Nnodes)+(1:prod(Nnodes)),:))));
            handles.evisu.xy=M\(phi'*(wdetJ*(dphidy*U((1:prod(Nnodes)),:))));
            handles.evisu.yx=M\(phi'*(wdetJ*(dphidx*U(prod(Nnodes)+(1:prod(Nnodes)),:))));
            %toc
            
        end
    end
    switch handles.fbasis
        case 'uni'
            load(fil0,'-mat','Up');
            handles.uni_model.Up=Up;
        case 'beam'
            load(fil0,'-mat','Up','Eax');
            handles.beam_model.Up=Up;
            handles.beam_model.Eax=Eax;
            if strcmp(handles.beam_model.beam_type,'timoshenko')
                load(fil0,'-mat','Esh');
                handles.beam_model.Esh=Esh;
            end
            
        case 'vic'
            handles.nscale=handles.vic_model.nscale;
        case 'fem'
            handles.nscale=handles.fem_model.nscale;
            %            handles.fem_model.zone={};
    end
    
    if (~strcmp(handles.fbasis,'uni'))&&(~strcmp(handles.fbasis,'vic'))&&strcmp(handles.param.analysis,'correlation')&&(~param.detect)
        if isfield(handles.fem_model,'degree')||(handles.param.ulive==1)
            handles.erroronelt=0;
        else
            fide=fopen([strrep(handles.param.result_file,'.res',''),'-error.res'],'r');
            try
                erroronelt=fread(fide,1);
                fclose(fide);
                handles.erroronelt=erroronelt;
                if erroronelt
                    ie=repmat((1:length(elt))',1,size(conn,2));
                    eton=sparse(conn,ie,val,prod(Nnodes),prod(Nelems));
                    eton=diag(sparse(1./sum(eton,2)))*eton;
                    handles.eton=eton;
                end
            catch
                %    keyboard
            end
        end
    end
    set(handles.message_text,'String','Loading result file.....done')
    set(handles.preview_button,'Enable','on')
    set(handles.preview_button,'State','on')
    if isfield(model,'visu')
        handles.fem_model=rmfield(handles.fem_model,'visu');
        handles.animation.iim=model.visu.iim;
        dview=(model.visu.view(1)==0)&&(model.visu.view(2)==90);
        if dview
            view(model.visu.view);
        end
        fields=fieldnames(model.visu);
        for ii=1:length(fields)
            handles=setfield(handles,fields{ii},getfield(model.visu,fields{ii}));
        end
    end
    handles=set_param(handles);
    preview_button_OnCallback(handles.preview_button,[],handles);
    handles=guidata(hObject);
else
    
    handles.animation.iim=0;
    switch handles.fbasis
        case 'fem'
            if ~strcmp([strrep(handles.param.result_file,'.res',''),'.vtk'],handles.fem_model.mesh_file)
                handles.fem_model.element_type=3;
            else
                handles.fem_model=rmfield(handles.fem_model,'mesh_file');
            end
            handles.nscale=handles.fem_model.nscale;
            handles.fem_model.nscale=1;
        case 'vic'
            handles.nscale=handles.vic_model.nscale;
            handles.vic_model.nscale=1;
    end
    handles.param.roi=[1,handles.sizeim(1),1,handles.sizeim(2)];
    
    sizeim=handles.sizeim(1:2);
    save(fullfile('TMP','sample0_0'),'sizeim','-v7.3');
    save(fullfile('TMP','sample0'),'sizeim','-v7.3');
    handles=set_param(handles);
    %     if handles.preview==1
    iscale=1;nmod=handles.nmod;
    LoadParameters(handles.param)
    LoadParameters(handles.fem_model,nmod)
    %ReferenceImage(nmod)
    sizeim=handles.sizeim(1:2);
    save(fullfile('TMP','sample0_0'),'sizeim','-v7.3');
    save(fullfile('TMP','sample0'),'sizeim','-v7.3');
    
    %         LoadMask(nmod);
    %         LoadMeshes(nmod);
    %         handles=run_preview(handles);
    %         handles=guidata(hObject);
    %         set(handles.mesh_type_edit,'Value',1)
    %         set(handles.mesh_type_edit,'Enable','on')
    %     else
    if handles.fem_model.mesh_type==3
        handles.mvisu.mesh_file= handles.fem_model.mesh_file;
        handles.mvisu.gluing_parameters=handles.fem_model.gluing_parameters;
    end
    handles=remesh(handles);
    %    end
    set(handles.message_text,'String','Loading input file.....done')
    
end
handles=set_param(handles);
handles=display_frame(handles);
switch handles.fbasis
    case 'uni'
        for ig=1:size(handles.uni_model.zone,2)
            plot_gage_data(handles,handles.uni_model.zone{2,ig},ig);
            handles=guidata(hObject);
        end
    case 'fem'
        %         if handles.ana==2
        %             for ig=1:size(handles.fem_model.zone,2)
        %                 gage=handles.fem_model.zone(:,ig);
        %                 xyp=gage{2};
        %                 gz=plot(xyp(:,1),xyp(:,2),'LineStyle','-','LineWidth',2);
        %                 gz.Tag=num2str(rand(1));
        %                 handles.fem_model.zone{3,ig}=gz;
        %                 guidata(hObject,handles);
        %                 switch handles.fem_model.zone{4,ig}
        %                     case 5
        %                         set(gz,'Color',[1,0,0],'uicontextmenu',handles.cmenu_export_crack);
        %                     otherwise
        %                         set(gz,'Color',[0,0.5,0],'uicontextmenu',handles.cmenu_export_gage);
        %                         %                    plot_gage_data(handles,handles.fem_model.zone{2,ig},ig);
        %                 end
        %
        %             end
        %
        %         end
        %        handles=guidata(hObject);
        
end
handles=set_param(handles);


function loading_ufr_file(hObject,handles,fil0,path0,FilterIndex)
error('not coded yet');
[param,model]=readINPFile(fil0);
handles.param=param;
if isfield(param,'deformed_image')
    nim=1;
    filref=param.reference_image;
    im00=(readim(filref));
    if iscell(param.deformed_image)
        nim=size(param.deformed_image,2);
    end
    handles.animation.nbstep=nim;
    handles.animation.frames=1:nim;
else
    filref=param.reference_image;
    reader=VideoReader(param.reference_image);
    handles.reader=reader;
    handles.param.video_sampling=max(1,ceil(reader.NumberOfFrames/100));
    
    im00=readim(reader,1);
    handles.animation.frames=2:handles.param.video_sampling:reader.NumberOfFrames;
    handles.animation.nbstep=length(handles.animation.frames);
end

if size(im00,3)==1%Cas d'une image NB
    im00=repmat(im00,[1,1,3]);
    rgb='MN';
else
    rgb='RGB';
end
handles.sizeim=size(im00);
set(0,'CurrentFigure',handles.figure1);
img_msg=sprintf('%s %dx%d %d-bit %s',filref,handles.sizeim(1:2),8*(1+double(max(im00(:))>255)),rgb);
set(handles.image_info_text,'String',img_msg);
if any(im00(:)>255)
    im00=double(im00);
    im00=255*(im00-min(im00(:)))/(max(im00(:))-min(im00(:)));
    %     if any(im00(:)>1023)
    %         if any(im00(:)>4095)
    %             im00=im00*(255/65535);
    %         else
    %             im00=im00*(255/4095);
    %         end
    %     else
    %         im00=im00*(255/1023);
    %     end
    im00=uint8(floor(im00));
end
gim=imagesc(permute(im00,[2,1,3]));
handles.gim=gim;
hold on;
axis equal
axis xy;
axis off
set(handles.axes1,'Clim',[0,255]);
set(handles.axes1,'ClimMode','manual');

colormap(handles.cmap);

handles.gcbar=colorbar('location','East','FontSize',20,'FontWeight','normal','LineWidth',1,'uicontextmenu',handles.cmenu_colorbar);
set(handles.gcbar,'Visible','off','Ycolor',[0,0,1]*0.8);
%set(handles.gcbar,'ButtonDownFcn','selectmoveresize');
assert(strcmp(model.basis,'fem'),'Invalid basis function for .ufr file loading !')
handles.fem_model=model;
handles.fbasis='fem';
handles.animation.iim=0;
if ~strcmp([strrep(handles.param.result_file,'.res',''),'.vtk'],handles.fem_model.mesh_file)
    handles.fem_model.element_type=3;
else
    handles.fem_model=rmfield(handles.fem_model,'mesh_file');
end
handles.nscale=handles.fem_model.nscale;
handles.ana=1;
handles.param.roi=[1,handles.sizeim(1),1,handles.sizeim(2)];

sizeim=handles.sizeim(1:2);
save(fullfile('TMP','sample0_0'),'sizeim','-v7.3');
save(fullfile('TMP','sample0'),'sizeim','-v7.3');
set(handles.figure1,'uicontextmenu',handles.cmenu_param);
handles=set_param(handles);
%     if handles.preview==1
%         iscale=1;nmod=handles.nmod;
%         LoadParameters(handles.param)
%         LoadParameters(handles.fem_model,nmod)
%         ReferenceImage(nmod)
%         LoadMask(nmod);
%         LoadMeshes(nmod);
%         handles=run_preview(handles);
%         handles=guidata(hObject);
%         set(handles.mesh_type_edit,'Value',1)
%         set(handles.mesh_type_edit,'Enable','on')
%     else
handles=remesh(handles);
%    end
set(handles.message_text,'String','Loading input file.....done')

handles=set_param(handles);
handles=display_frame(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function new_fem_analysis_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to new_fem_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer=inputdlg({'Along X:','Along Y:','Unit:'},'Stekch size',1,...
    {num2str(1000),num2str(1000),num2str(1.e-3)});
if ~isempty(answer)
    handles=reset_param(handles);
    sizeim=[eval(answer{1}),eval(answer{2})];
    handles.param.pixel_size=eval(answer{3});
    dx=max(1,round(min(sizeim)/100));
    dx=max(1,5*round(dx/5));
    handles.gstep=dx;
    handles.ana=3;
    handles.nscale=1;
    handles.sizeim=sizeim;
    handles.fbasis='fem';
    handles.param.analysis='mechanics';
    handles.fem_model.mesh_size=(round(mean(handles.sizeim(1:2)/50))*[1,1]);
    handles.fem_model.element_type=4;
    handles.fem_model.mesh_type=2;
    handles.fem_model.phantom_nodes=0;
    handles.fem_model.material_model='elastic_homogeneous_isotropic';
    matmod.young=200e9;
    matmod.nu=0.3;
    matmod.KIc=100e6;
    handles.fem_model.material_parameters=matmod;
    set(0,'CurrentFigure',handles.figure1);
    handles.param.roi=[1,handles.sizeim(1),1,handles.sizeim(2)];
    gim=imagesc(ones([sizeim,3]));
    handles.gim=gim;
    hold on
    axis equal
    axis xy
    axis off
    
    set(handles.axes1,'Clim',[0,255]);
    set(handles.axes1,'ClimMode','manual');
    set(handles.axes1,'GridLineStyle','--');
    set(handles.axes1,'XLim',[handles.param.roi(1:2)]);
    set(handles.axes1,'YLim',[handles.param.roi(3:4)]);
    set(handles.axes1,'Xgrid','on');
    set(handles.axes1,'Ygrid','on');
    LoadParameters(handles.param)
    LoadParameters(handles.fem_model,handles.nmod)
    save(fullfile('TMP','sample0_0'),'sizeim');
    save(fullfile('TMP','sample0'),'sizeim');
    handles=remesh(handles);
    handles=set_param(handles);
    handles=display_frame(handles);
end

function start_new_analysis(hObject,handles,fil0,path0,FilterIndex)
handles.ana=1;
switch FilterIndex
    case 1
        if iscell(fil0)
            %             if ispc
            %                 fil0=fliplr(fil0);
            %             end
            fils=SortImageFiles(fil0);
            fil0=fils;
            filref=fil0{1,1};
            if size(fils,1)==2
                handles.stereo=1;
            else
                fil0(:,1)=[];
            end
            handles.param.reference_image=filref;
            im00=(readim(filref));
            handles.animation.nbstep=size(fil0,2)-handles.stereo;
            handles.animation.frames=1:size(fil0,2)-handles.stereo;
            if size(fil0,2)>2
                handles.param.restart=0;
                handles.nscale=3;
            else
                handles.param.restart=1;
                handles.nscale=3;
            end
            if numel(fil0)>1
                handles.param.deformed_image=fil0;
            else
                handles.param.deformed_image=fil0{1};
            end
            set(handles.message_text,'String','Loading images.....done')
        else
            filref=fil0;
            handles.param.reference_image=filref;
            im00=(readim(filref));
            handles.animation.nbstep=1;
            handles.animation.frames=1;
            handles.param.restart=0;
            handles.nscale=3;
            handles.param.deformed_image=fil0;
            handles.param.ulive=1;
            set(handles.message_text,'String','Loading image.....done')
            
            %           error('not coded yet')
        end
    case 2
        handles.param=rmfield(handles.param,'deformed_image');
        set(handles.message_text,'String','Loading video.....')
        reader=VideoReader(fil0);
        set(handles.message_text,'String','Loading video.....done')
        filref=fil0;
        handles.param.reference_image=filref;
        im00=readim(reader,1);
        handles.param.restart=0;
        handles.nscale=1;
        handles.param.number_of_frames=reader.NumberOfFrames;
        handles.param.video_sampling=max(1,ceil(reader.NumberOfFrames/100));
        handles.reader=reader;
        handles.animation.frames=2:handles.param.video_sampling:reader.NumberOfFrames;
        handles.animation.nbstep=length(handles.animation.frames);
end
if size(im00,3)==1%Cas d'une image NB
    im00=repmat(im00,[1,1,3]);
    rgb='MN';
else
    rgb='RGB';
end
handles.sizeim=size(im00);
handles.fem_model.mesh_size=max(16,round(mean(handles.sizeim(1:2)/50))*[1,1]);
set(0,'CurrentFigure',handles.figure1);
handles.param.roi=[1,size(im00,1),1,size(im00,2)];
img_msg=sprintf('%s %dx%d %d-bit %s',filref,handles.sizeim(1:2),8*(1+double(max(im00(:))>255)),rgb);
set(handles.image_info_text,'String',img_msg);
if any(im00(:)>255)
    im00=double(im00);
    im00=255*(im00-min(im00(:)))/(max(im00(:))-min(im00(:)));
    %     if any(im00(:)>1023)
    %         if any(im00(:)>4095)
    %             im00=im00*(255/65535);
    %         else
    %             im00=im00*(255/4095);
    %         end
    %     else
    %         im00=im00*(255/1023);
    %     end
    im00=uint8(floor(im00));
end
gim=imagesc(permute(im00,[2,1,3]));
handles.gim=gim;
hold on;
axis equal
axis xy;
axis off

set(handles.axes1,'Clim',[0,255]);
set(handles.axes1,'ClimMode','manual');

colormap(handles.cmap);

handles.gcbar=colorbar('location','East','FontSize',20,'FontWeight','normal','LineWidth',1,'uicontextmenu',handles.cmenu_colorbar);
set(handles.gcbar,'Visible','off','Ycolor',[0,0,1]*0.8);
%set(handles.gcbar,'ButtonDownFcn','selectmoveresize');
handles=set_param(handles);
if handles.preview==1
    handles=run_preview(handles);
    set(handles.mesh_type_edit,'Value',1)
    set(handles.mesh_type_edit,'Enable','on')
    handles=set_param(handles);
end
LoadParameters(handles.param);
LoadParameters(handles.uni_model,handles.nmod);
%ReferenceImage();
sizeim=handles.sizeim(1:2);
save(fullfile('TMP','sample0_0'),'sizeim','-v7.3');
save(fullfile('TMP','sample0'),'sizeim','-v7.3');
set(handles.figure1,'uicontextmenu',handles.cmenu_param);

function handles=run_preview(handles)
set(handles.message_text,'String','Running preview...');
sizeim=handles.sizeim(1:2);
nscale=max(3,ceil(log(max(sizeim)/512)/log(2)));
hmax=round((min(sizeim(1:2)/25)));
hmax=hmax+(hmax==round(hmax/10)*10);
model=handles.fem_model;
model.mesh_size=ones(2,1)*hmax;
if isfield(model,'mesh_file'),model=rmfield(model,'mesh_file');end
param=handles.param;
model.nscale=nscale;
param.regularization_type='tiko';
param.regularization_parameter=16/hmax;
param.preview=1;
param.restart=1;
nmod=2;
LoadParameters(param);
LoadParameters(model,nmod);
ReferenceImage(nmod);
LoadMeshes(nmod);
LoadMask(nmod);
CreateBasisFunction(nscale,nmod);
ComputeGradFPhi(nscale,nmod);
AssembleCorrelationOperator(nscale,nmod);

CreateGradBasisFunction(nscale,nmod);

load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,nscale-1)),'xo','yo','Nnodes','Nelems','conn','elt');
xg=2^(nscale-1)*(xo-0.5)+0.5;
yg=2^(nscale-1)*(yo-0.5)+0.5;
meshg.Nnodes=Nnodes;
meshg.Nelems=Nelems;
meshg.conn=conn;
meshg.elt=elt;
meshg.xo=xg;
meshg.yo=yg;
meshg.h=hmax;
iscale=1;
load(fullfile('TMP',sprintf('sample0_%d',iscale-1)),'im0');
[S]=GetEntropy(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),im0,1);
load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','Nnodes','Nelems','conn','elt');
% entropy=zeros(Nnodes);
% inde=repmat((1:prod(Nelems))',1,size(conn,2));
% for in=1:prod(Nnodes)
%     entropy(in)=mean(S(inde(conn==in)));
% end
% entropy=(entropy)/(max(entropy(:)));
S=S/max(S);

handles.preview=1;
handles.mpreview=meshg;

if ~isfield(handles,'mvisu')
    handles.mvisu.Nelems=Nelems;
    handles.mvisu.Nnodes=Nnodes;
    handles.mvisu.conn=conn;
    handles.mvisu.elt=elt;
    handles.mvisu.xo=xo;
    handles.mvisu.yo=yo;
    handles.mvisu.entropy=S;
    handles.mvisu.smin=0.5;
    ie=repmat((1:length(elt))',1,size(conn,2));
    eton=sparse(conn,ie,1,prod(Nnodes),prod(Nelems));
    eton=diag(sparse(1./sum(eton,2)))*eton;
    handles.eton=eton;
end

Uini=[];

handles.upreview=0;%zeros(2*prod(Nnodes),size(param.deformed_image,2));
handles.animation.iim=0;
handles.animation.playing=1;
handles.fem_model.mesh_size=hmax*ones(2,1);
handles=display_frame(handles);
iim=1;
while iim<=handles.animation.nbstep
    tic;
    [U]=SolvePreview(Uini,nscale,nmod,iim);
    Uini=U;
    handles.upreview=U;
    handles.animation.iim=iim;
    handles=display_frame(handles);
    if handles.animation.playing
        iim=iim+1;
    else
        iim=size(Uini,2);
    end
    pause(max(0.01,0.05-toc))
end
[dphidx,dphidy]=CreateGradFiniteElementBasis(fullfile('TMP',sprintf('%d_mesh_%d',nmod,nscale-1)),handles.sizeim,1,[],'Gauss_points');
[phi]=CreateFiniteElementBasis(fullfile('TMP',sprintf('%d_mesh_%d',nmod,nscale-1)),handles.sizeim,1,[],'Gauss_points');
[wdetJ]=GetWeigthDetJ(fullfile('TMP',sprintf('%d_mesh_%d',nmod,nscale-1)),handles.sizeim,1,'Gauss_points');
M=phi'*(wdetJ*phi);
R=dphidx'*(wdetJ*dphidx)+dphidy'*(wdetJ*dphidy);

lc=(min(handles.sizeim)/10);
lm=sqrt((max(xo)-min(xo))*(max(yo)-min(yo))/length(elt))/10;
V=cos(2*pi*(xg)/lc).*cos(2*pi*(yg)/lc);
a=(V'*M*V)/(V'*R*V);
a=a*(pi*lm/lc)^2;
M=M+a*R;
%handles.mpreview.dphidx=M\(phi'*(wdetJ*dphidx));
%handles.mpreview.dphidy=M\(phi'*(wdetJ*dphidy));
handles.epreview.xx=M\(phi'*(wdetJ*(dphidx*U((1:prod(meshg.Nnodes)),:))));
handles.epreview.yy=M\(phi'*(wdetJ*(dphidy*U(prod(meshg.Nnodes)+(1:prod(meshg.Nnodes)),:))));
handles.epreview.xy=M\(phi'*(wdetJ*(dphidy*U((1:prod(meshg.Nnodes)),:))));
handles.epreview.yx=M\(phi'*(wdetJ*(dphidx*U(prod(meshg.Nnodes)+(1:prod(meshg.Nnodes)),:))));

handles.animation.playing=0;
guidata(handles.figure1,handles)

set(handles.message_text,'String','Preview run...');


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function handles=basis_functions(handles)
switch handles.fbasis
    case 'uni'
        name='strain gage';
    case 'beam'
        handles.preview=0;
        name='BEAM';
    case 'vic'
        handles.preview=0;
        name='VIC';
    case 'fem'
        if handles.preview
            handles.ondefimage=0;
        else
            handles.animation.iim=0;
        end
        handles=remesh(handles);
        if isfield(handles.fem_model,'degree')
            name='NURBS';
            if handles.param.regularization_parameter==0
                handles.param.regularization_parameter=round(handles.fem_model.mesh_size(1)/4);
            end
            if strcmp(handles.param.regularization_type,'none')
                handles.param.regularization_type='tiko';
            end
        else
            name='FEM';
        end
    case 'none'
        name='none';
end
set(handles.message_text,'String',sprintf('Selecting basis function: %s.....done',name))
% if handles.preview||strcmp(handles.fbasis,'fem')
%     handles=plot_mesh(handles);
% end
handles=display_frame(handles);


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function sampling_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sampling_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampling_edit as text
%        str2double(get(hObject,'String')) returns contents of sampling_edit as a double
%handles=guidata(hObject);
frq=eval(get(hObject,'String'));
handles.param.video_sampling=frq;
frames=2:frq:handles.param.number_of_frames+1;
handles.animation.frames=frames;
handles.animation.nbstep=length(frames);
handles=set_param(handles);
set(handles.message_text,'String','Setting sampling frequency.....done')

% --- Executes during object creation, after setting all properties.
function sampling_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampling_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nbf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nbf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbf_edit as text
%        str2double(get(hObject,'String')) returns contents of nbf_edit as a double
%handles=guidata(hObject);
nbf=eval(get(hObject,'String'));
nbf=min(nbf,handles.reader.NumberOfFrames);
handles.param.number_of_frames=nbf;
frames=2:handles.param.video_sampling:nbf;
handles.animation.nbstep=length(frames);
handles=set_param(handles);
set(handles.message_text,'String','Setting last of frame to analyse.....done')


% --- Executes during object creation, after setting all properties.
function nbf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mesh_type_edit.
function handles=remesh(handles)
iscale=1;nmod=handles.nmod;
ho=handles.fem_model.mesh_size;
switch handles.fem_model.mesh_type
    case 1 %built in
        
        if isfield(handles.fem_model,'mesh_file')
            handles.fem_model=rmfield(handles.fem_model,'mesh_file');
        end
        if isfield(handles.fem_model,'gluing_parameters')
            handles.fem_model=rmfield(handles.fem_model,'gluing_parameters');
        end
    case 2
        if isfield(handles.fem_model,'mesh_file')
            handles.fem_model=rmfield(handles.fem_model,'mesh_file');
        end
        if isfield(handles.fem_model,'gluing_parameters')
            handles.fem_model=rmfield(handles.fem_model,'gluing_parameters');
        end
        handles.fem_model.mesh_size=(ho/1);
    case 3 %external
        nmod=1;
        if ~isfield(handles.mvisu,'mesh_file')
            [filename,pathres,FilterIndex]=uigetfile('*.vtk','VTK Mesh File (*.vtk)');
            if ~exist(filename,'file'),copyfile(fullfile(pathres,filename),filename);end
            handles.fem_model.mesh_file=filename;
            handles.mvisu.mesh_file=filename;
        else
            handles.fem_model.mesh_file=handles.mvisu.mesh_file;
        end
        handles=gluemeshonimage(handles);
        handles.fem_model.gluing_parameters=handles.mvisu.gluing_parameters;
        handles.fem_model.element_type=4;
        ReferenceImage(nmod);
end
LoadParameters(handles.param)
LoadParameters(handles.fem_model,nmod)
LoadMeshes(nmod);
LoadMask(nmod);

load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'Nnodes','Nelems','conn','elt','xo','yo');
handles.fem_model.mask=1;
handles.mvisu.smin=handles.fem_model.smin;
handles.mvisu.xo=xo;
handles.mvisu.yo=yo;
handles.mvisu.conn=conn;
handles.mvisu.Nnodes=Nnodes;
handles.mvisu.Nelems=Nelems;
handles.mvisu.elt=elt;

if handles.fem_model.mesh_type<3
    if handles.ana==1
        if prod(handles.sizeim(1:2))<8e6
            im0=get(handles.gim,'CData');
            im0=permute(mean(im0,3),[2,1]);
            [entropy]=GetEntropy(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),im0,1);
            %    Smax=prod(handles.fem_model.meshsize)
            entropy=(entropy-0*min(entropy(:)))/(max(entropy(:))-0*min(entropy(:)));
            
            handles.mvisu.entropy=reshape(entropy,Nelems);
        else
            handles.mvisu.entropy=1;
            set(handles.cmenu_mesh_entropy,'Enable','off');
        end
    end
    switch handles.fem_model.mesh_type
        case 1
            if numel(handles.mvisu.entropy)>1
                conn=conn(:,1:handles.fem_model.element_type);
                ie=repmat((1:length(elt))',1,size(conn,2));
                eton=sparse(conn,ie,1,prod(Nnodes),prod(Nelems));
                eton=diag(sparse(1./sum(eton,2)))*eton;
                handles.eton=eton;
            end
        case 2
            handles=GenUnstructuredMesh(handles);
            handles.fem_model.mesh_size=ho;
    end
end


handles=UpdateBCSZone(handles);
handles=plot_mesh(handles);

function handles=UpdateBCSZone(handles)
toremove=0;
for iz=1:size(handles.fem_model.zone,2)
    if handles.fem_model.zone{4,iz}==6
        zone=handles.fem_model.zone{2,iz};
        in=inpolygon(handles.mvisu.xo+handles.param.roi(1)-1,handles.mvisu.yo+handles.param.roi(3)-1,zone(:,1),zone(:,2));
        in=find(in);
        if ~isempty(in)
            xon=handles.mvisu.xo(in)+handles.param.roi(1)-1;
            yon=handles.mvisu.yo(in)+handles.param.roi(3)-1;
            handles.fem_model.zone{5,iz}=[xon,yon];
            handles.fem_model.zone{6,iz}=in;
        else
            toremove(iz)=1;
        end
    end
end
if any(toremove)
    toremove=find(toremove);
    delete(handles.fem_model.zone{3,toremove});
end
function handles=gluemeshonimage(handles)
hObject=handles.figure1;
if handles.animation.iim
    handles.animation.iim=0;
    handles=display_frame(handles);
end
[xo,yo,zo,conn,elt,selected]=ReadVTK(handles.fem_model.mesh_file);
try delete(handles.gmesh); catch, end
try delete(handles.gedge); catch, end
roi=handles.param.roi;
sizeim=handles.sizeim(1:2);
if ~isfield(handles.mvisu,'gluing_parameters')
    mxo=mean(xo);
    myo=mean(yo);
    lx=max(xo)-min(xo);
    ly=max(yo)-min(yo);
    tx=roi(1)+0.5*sizeim(1);
    ty=roi(3)+0.5*sizeim(2);
    scale=min(diff(roi(1:2))/lx,diff(roi(3:4))/ly);
    angl=0;
    glues{1,:}={'translate',[-mxo;-myo;0]};
    glues{2,:}={'scale',[0;0;0],scale};
    glues{3,:}={'rotate',[0;0;1],[0;0;0],angl};
    glues{4,:}={'translate',[tx;ty;0]};
    handles.mvisu.gluing_parameters=glues;
else
    glues=handles.mvisu.gluing_parameters;
end
guidata(handles.figure1,handles);
[xo,yo,~]=GlueMesh(glues,xo,yo);
seg3=reshape(conn(elt==3,[1,2,1,3,2,3])',2,3*sum(elt==3))';
seg4=reshape(conn(elt==4,[1,2,2,3,3,4,4,1])',2,4*sum(elt==4))';

seg=unique(sort([seg3;seg4],2),'rows');

xo=xo(seg);
yo=yo(seg);
handles.gedge=plot(xo'+roi(1)-1,yo'+roi(3)-1,'r-');
set(handles.message_text,'String','Left click to rotate, scroll to scale, middle click to translate, right click to end')
handles.gadjust.xo=xo;
handles.gadjust.yo=yo;
guidata(hObject,handles);
set(hObject,'WindowScrollWheelFcn',@mesh_scroll_scale);
set(hObject,'WindowButtonDownFcn',@mesh_adjust)
set(hObject,'KeyPressFcn',@mesh_fixe)
uiwait
set(gcf,'WindowButtonUpFcn','');
set(gcf,'WindowButtonMotionFcn',@show_position);
set(gcf,'WindowButtonDownFcn',@double_clic);
set(gcf,'WindowScrollWheelFcn','');
set(hObject,'KeyPressFcn','')

set(handles.message_text,'String','Adjusting mesh.....done')
try delete(handles.gedge);catch, end
handles=guidata(hObject);

function mesh_adjust(hObject,eventdata)
if ~strcmp(get(gcf,'selectiontype'),'alt')
    handles=guidata(hObject);
    pt=getposition(handles);
    handles.initial_point=pt(1,1:2);
    guidata(hObject,handles);
    set(gcf,'WindowButtonUpFcn',@mesh_stop_follow);
else
    uiresume
end
function mesh_fixe(hObject,eventdata)
uiresume



function mesh_stop_follow(hObject,eventdata)
set(hObject,'WindowButtonMotionFcn',@show_position);
handles=guidata(hObject);
po=handles.initial_point;
pt=getposition(handles);

xp=handles.gadjust.xo;
yp=handles.gadjust.yo;
try delete(handles.gedge);catch end
switch get(handles.figure1,'selectiontype')
    case 'extend'
        xp=xp+pt(1,1)-po(1);
        yp=yp+pt(1,2)-po(2);
        Txy=handles.mvisu.gluing_parameters{4,:};
        Txy=Txy{2};
        handles.mvisu.gluing_parameters{4,:}={'translate',[Txy(1)+pt(1,1)-po(1);Txy(2)+pt(1,2)-po(2);0]};
    case 'normal'
        po=po*[1;1i];
        pt=pt*[1;1i];
        Zcp=xp+1i*yp;
        pc=mean(Zcp(:));
        rot=((pt-pc)/abs(pt-pc)*abs(po-pc)/(po-pc));
        Zcp=(Zcp-pc)*rot+pc;
        a=handles.mvisu.gluing_parameters{3,:};
        a=a{4};
        handles.mvisu.gluing_parameters{3,:}={'rotate',[0;0;1],[0;0;0],a-angle(rot)};
        xp=real(Zcp);
        yp=imag(Zcp);
        
end
handles.gedge=plot(xp',yp','r-');
%set(handles.gedge,'Xdata',xp(:),'Ydata',yp(:));
handles.gadjust.xo=xp;
handles.gadjust.yo=yp;
guidata(handles.figure1,handles);
set(gcf,'WindowButtonUpFcn','');

function mesh_scroll_scale(hObject,eventdata)
handles=guidata(hObject);
scale=handles.mvisu.gluing_parameters{2,:};
scale=scale{3};
scale=max(0,scale*(1-0.005*eventdata.VerticalScrollCount));
Txy=handles.mvisu.gluing_parameters{4,:};
Txy=Txy{2};
handles.mvisu.gluing_parameters{2,:}={'scale',[0;0;0],scale};
xp=handles.gadjust.xo-Txy(1);
yp=handles.gadjust.yo-Txy(2);
xp=xp*max(0,(1-0.005*eventdata.VerticalScrollCount))+Txy(1);
yp=yp*max(0,(1-0.005*eventdata.VerticalScrollCount))+Txy(2);
try delete(handles.gedge);catch end
handles.gedge=plot(xp',yp','r-');
handles.gadjust.xo=xp;
handles.gadjust.yo=yp;
guidata(handles.figure1,handles);

function mesh_follow_mouse(hObject,eventdata)
handles=guidata(hObject);
po=handles.initial_point;
pt=getposition(handles);

xp=handles.gadjust.xo;
yp=handles.gadjust.yo;
get(handles.figure1,'selectiontype')
try delete(handles.gedge);catch end
switch get(handles.figure1,'selectiontype')
    case 'extend'
        xp=xp+pt(1,1)-po(1);
        yp=yp+pt(1,2)-po(2);
    case 'normal'
        po=po*[1;1i];
        pt=pt*[1;1i];
        Zcp=xp+1i*yp;
        pc=mean(Zcp(:));
        rot=((pt-pc)/abs(pt-pc)*abs(po-pc)/(po-pc));
        Zcp=(Zcp-pc)*rot+pc;
        xp=real(Zcp);
        yp=real(Zcp);
        
end
handles.gedge=plot(xp',yp','r-');
guidata(hObject,handles);

function handles=GenUnstructuredMesh(handles)
set(handles.message_text,'String','Unstructured mesh optimization.....')
if handles.ana==1
    if numel(handles.mvisu.entropy)>1
        mask=handles.mvisu.entropy>handles.fem_model.smin;
        nf=3;
        mask=filter2(ones(nf)/nf/nf,double(mask));
        mask=double(mask>0);
    else
        mask=[];
    end
else
    mask=[];
end
[xo,yo,conn,indc]=SmoothMesh(handles,-2*mask+1);
Nnodes=size(xo);
Nelems=[size(conn,1),1,1];
elt=3*ones(size(conn,1),1);
set(handles.message_text,'String','Unstructured mesh optimization.....done')
handles.mvisu.xo=xo;
handles.mvisu.yo=yo;
handles.mvisu.conn=conn;
handles.mvisu.Nnodes=Nnodes;
handles.mvisu.Nelems=Nelems;
handles.mvisu.elt=elt;
if ~isempty(indc)
    for iz=1:size(handles.fem_model.zone,2)
        
        if handles.fem_model.zone{4,iz}==5
            
            handles.fem_model.zone{10,iz}=indc(:,iz);
            
        end
    end
end


function [xo,yo,conn,indc]=SmoothMesh(handles,im)
%roi=[1,size(im,1),1,size(im,2)];
h=handles.fem_model.mesh_size;
roi=handles.param.roi;
roi=roi+0.5*[h(1),-h(1),h(2),-h(2)];
xe=zeros(prod(handles.mvisu.Nelems),1);
ye=zeros(prod(handles.mvisu.Nelems),1);

for ie=1:prod(handles.mvisu.Nelems)
    xe(ie)=mean(handles.mvisu.xo(handles.mvisu.conn(ie,1:handles.mvisu.elt(ie))));
    ye(ie)=mean(handles.mvisu.yo(handles.mvisu.conn(ie,1:handles.mvisu.elt(ie))));
end


xe=reshape(xe,handles.mvisu.Nelems);
ye=reshape(ye,handles.mvisu.Nelems);



% d=GetSignedDistanceToZone(handles,xe(:),ye(:));
%
%  if numel(d)>1
%  figure
%  imagesc(reshape(d,size(xe)))
%  hold on
%  contour(reshape(d,size(xe)),[0,0],'k')
%  axis equal
%  end
[d,xyfixo,indc]=GetMeshDensity(handles.fem_model,xe(:),ye(:));

%  if numel(d)>1
%     figure
%     imagesc(reshape(d,size(xe))')
%     colorbar
%     axis equal
%     hold on
%      if ~isempty(xyfix)
%              figure
%          plot(xyfix(:,2),xyfix(:,1),'rx')
%      end
% end
if ~isempty(im)
    im=double(im);
    if any(im>0)
        lvl=LSReinit(im,4*mean(h),mean(h));
        lvl=lvl+mean(h);
    else
        lvl=0*im-10000000000;
    end
    %      figure
    %      imagesc(lvl)
    %      hold on
    %      contour(lvl,[0,0],'k')
    %      axis equal
    
    ld=@(xy) max(mexInterpLVL7((xy(:,1)-min(xe(:)))/h(1)+1,(xy(:,2)-min(ye(:)))/h(2)+1,lvl),GetSignedDistanceToZone(handles.fem_model,roi,xy(:,1),xy(:,2)));
else
    ld=@(xy) GetSignedDistanceToZone(handles.fem_model,roi,xy(:,1),xy(:,2));
end
%hmin=(min(d)-0*std(d));
lh=@(xy) (GetMeshDensity(handles.fem_model,xy(:,1),xy(:,2)));
[xo,yo,conn,xyfix]=GenMeshFromLVL7(roi,mean(h)*min(d),ld,lh,xyfixo);
% couls='rbrbrbrbrbrbrbrb';
% simb='sosososososososo';
% figure
% hold on
% %triplot(conn,xo,yo)
% %plot(xyfixo*[1;1i],'kx')
% plot(xyfix*[1;1i],'kx')
if ~isempty(indc)
    out=feval(ld,xyfix)>=(0.001*mean(h)*min(d));
    %    out=feval(ld,xyfix)>0.75*2*mean(h)*hmin*feval(lh,xyfix);
    %    xyfix=xyfix(~out,:);
    new_id=0*out;
    new_id(~out)=1:sum(~out);
    %plot(xyfix*[1;1i],'ms')
    %plot([xo(1:sum(~out)),yo(1:sum(~out))]*[1;1i],'g*')
    
    
    duplicated=[];
    for iz=1:size(handles.fem_model.zone,2)
        if handles.fem_model.zone{4,iz}==5
            for id=1:size(indc,1)
                nodes=indc{id,iz}';
                if ~isempty(nodes)
                    for in=1:numel(nodes)
                        dist=abs((xyfix(nodes(in),1)-xo)+1i*(xyfix(nodes(in),2)-yo));
                        [dmin,idmin]=min(dist);
                        if dmin>(0.001*mean(h)*min(d))
                            nodes(in)=0;
                        else
                            nodes(in)=idmin;
                        end
                    end
                    %                    nodes=new_id(nodes);
                    if id==1,tipin=[nodes(1),nodes(end)]>0;end
                                        nodes(nodes==0)=[];
                    segc=[nodes(1:end-1),nodes(2:end)];
                    eneighboor=GetEltsFromNodes(conn,3*ones(size(conn,1),1),nodes);
                    neighboor=conn(eneighboor,:);
                    segs=reshape(neighboor(:,[1,2,1,3,2,3])',2,3*size(neighboor,1))';
                    segs=unique(sort(segs,2),'rows');
                    
                    [~,isc,~]=intersect(segc,segs,'rows');
                    segc(isc,:)=[];
                    toinsert=zeros(size(segc,1),1);
                    for ic=1:size(segc,1)
                        has1=neighboor((sum(neighboor==segc(ic,1),2)>0),:);
                        has2=neighboor((sum(neighboor==segc(ic,2),2)>0),:);
                        nodes12=intersect(has1(:),has2(:))';
                        switch numel(nodes12)
                            case 2
                                for ii=1:2
                                    elti=sort([segc(ic,:),nodes12(ii)],2);
                                    try
                                        [~,ie,~]=intersect(sort(neighboor,2),sort([nodes12,segc(ic,ii)],2),'rows');
                                    catch
                                        keyboard
                                    end
                                    detJ=((xo(elti(2))-xo(elti(1)))*(yo(elti(3))-yo(elti(1)))-(yo(elti(2))-yo(elti(1)))*(xo(elti(3))-xo(elti(1))))/2;
                                    if detJ<0
                                        elti([1,2])=elti([2,1]);
                                    end
                                    conn(eneighboor(ie),:)=elti;
                                end
                            case 1
                                toinsert(ic)=nodes12;
                        end
                    end
                    if any(toinsert)
                        nnodes=zeros(length(nodes)+sum(toinsert>0),1);
                        for ic=1:size(segc,1)
                            if toinsert(ic)>0
                                in=find(nodes==segc(ic,1));
                                nnodes(in+1)=toinsert(ic);
                            end
                        end
                        in=1;
                        for ii=1:length(nnodes)
                            if nnodes(ii)==0
                                nnodes(ii)=nodes(in);
                                in=in+1;
                            end
                        end
                        nodes=nnodes;
                    end
                    
                    
                    if id==1
                        %                        t=gradient(xo(nodes)+1i*yo(nodes));
                        nodes=repmat(nodes,1,2);
                        toadd=(1+tipin(1)):(size(nodes,1)-tipin(2)*handles.fem_model.zone{9,iz});
                        if ~isempty(duplicated)
                            ldup=[];
                            inn=numel(toadd)-1;
                            inod=nodes(toadd(inn),1);
                            elts=sum(conn==inod,2)>0;inods=conn(elts,:);
                            
                            for ii=1:numel(inods)
                                if any(duplicated(:,2)==inods(ii))
                                    try
                                        nodes(toadd(end),:)=duplicated(duplicated(:,1)==nodes(toadd(end),1),2);
                                    catch
                                    end
                                    break
                                end
                            end
                        end
                        nodes(toadd,2)=length(xo)+(1:length(toadd));
                        xo=[xo;xo(nodes(toadd,1))];
                        yo=[yo;yo(nodes(toadd,1))];
                        duplicated=[duplicated;nodes(toadd,:)];
                        %                       plot(xo(nodes(toadd,1)),yo(nodes(toadd,1)),[couls(iz),simb(iz)])
                        for ip=1:length(toadd)
                            %                             inods=nodes(toadd(ip),1);
                            %                             seg=-t(toadd(ip));
                            %                             mid=(xo(inods)+1i*yo(inods));
                            %                             elts=sum(conn==inods,2)>0;
                            %                             xc=sum(xo(conn(elts,:)),2)/3;
                            %                             yc=sum(yo(conn(elts,:)),2)/3;
                            %                             side=real(((xc+1i*yc)-mid)*(seg*exp(1i*pi/2))');
                            %                             elts=find(elts);
                            %
                            %
                            %                             for ii=1:length(side)
                            %                                 if side(ii)<0
                            %                                     nelt=conn(elts(ii),:);
                            %                                     nelt(nelt==inods)=nodes(toadd(ip),2);
                            %                                     conn(elts(ii),:)=nelt;
                            %                                 end
                            %
                            %                             end
                            
                            
                            inods=nodes(toadd(ip),1);
                            zn=(xo(inods)+1i*yo(inods));
                            elts=sum(conn==inods,2)>0;
                            elts=find(elts);
                            for ii=1:length(elts)
                                xc=sum(xo(conn(elts(ii),:)))/3;
                                yc=sum(yo(conn(elts(ii),:)))/3;
                                
                                if toadd(ip)==1
                                    behind=0;
                                    infront=1;
                                elseif toadd(ip)==size(nodes,1)
                                    behind=1;
                                    infront=0;
                                else
                                    t=diff(xo(nodes(toadd(ip)+(0:1),1))+1i*yo(nodes(toadd(ip)+(0:1),1)));
                                    infront=(real(((xc+1i*yc)-zn)*(t)'))>0;
                                    t=diff(xo(nodes(toadd(ip)+(-1:0),1))+1i*yo(nodes(toadd(ip)+(-1:0),1)));
                                    behind=(real(((xc+1i*yc)-zn)*(t)'))<0;
                                end
                                if infront
                                    t=diff(xo(nodes(toadd(ip)+(0:1),1))+1i*yo(nodes(toadd(ip)+(0:1),1)));
                                elseif behind
                                    t=diff(xo(nodes(toadd(ip)+(-1:0),1))+1i*yo(nodes(toadd(ip)+(-1:0),1)));
                                else
                                    t=diff(xo(nodes(toadd(ip)+(0:1),1))+1i*yo(nodes(toadd(ip)+(0:1),1)));
                                    t=t+diff(xo(nodes(toadd(ip)+(-1:0),1))+1i*yo(nodes(toadd(ip)+(-1:0),1)));
                                end
                                n=-t*exp(1i*pi/2);
                                side=real(((xc+1i*yc)-zn)*(n)');
                                if side<0
                                    %                                    'duplicated'
                                    nelt=conn(elts(ii),:);
                                    nelt(nelt==inods)=nodes(toadd(ip),2);
                                    conn(elts(ii),:)=nelt;
                                    %                                    plot(xo(nodes(toadd(ip),1)),yo(nodes(toadd(ip),1)),[couls(iz),'o'])
                                end
                                
                            end
                            
                        end
                    else
                        if ~isempty(nodes)
                            if tipin(mod(id,2)+1)
                                face_nodes=indc{1,iz};
                                tipid=(mod(id,2)==0)+size(face_nodes,1)*(mod(id,2)==1);
                                ztip=xo(face_nodes(tipid,1))+1i*yo(face_nodes(tipid,1));
                                zzone=xo(nodes)+1i*yo(nodes);
                                zface=xo(face_nodes(:,1))+1i*yo(face_nodes(:,1));
                                rtip=mean(abs(zzone-ztip));
                                [dmin,idmin]=min(abs(abs(zface-ztip)-rtip));
                                nodes(end+1)=face_nodes(idmin,1);
                                nodes(end+1)=face_nodes(idmin,2);
                            end
                        end
                    end
                    indc{id,iz}=nodes;
                end
            end
        end
    end
    
end


function mask=GetMask(handles)
mask=1;
if handles.ana==1
    switch handles.fbasis
        case 'fem'
            if (handles.fem_model.mesh_type==1)
                xo=handles.mvisu.xo;
                yo=handles.mvisu.yo;
                %                 Nnodes=handles.mvisu.Nnodes;
                %                 Nelems=handles.mvisu.Nelems;
                %                 conn=handles.mvisu.conn;
                %                elt=handles.mvisu.elt;
                if numel(handles.mvisu.entropy)>1
                    s=handles.eton*handles.mvisu.entropy(:);
                    
                    
                    mask=((handles.fem_model.mask)&(s(:)>=handles.mvisu.smin));
                else
                    mask=1;
                end
                %                 mask=reshape(mask,Nelems);
                %                 h=handles.fem_model.mesh_size;
                %                 nf=3;
                %                 mask=filter2(ones(nf)/nf/nf,double(mask));
                %                 mask=double(mask>0);
                %                 lvl=LSReinit(-2*mask+1,2*mean(h),mean(h));
                %                 mask=double(lvl(:)+mean(h)<0);
                %
                %                 xe=zeros(prod(Nelems),1);
                %                 ye=zeros(prod(Nelems),1);
                %
                %                 for ie=1:prod(Nelems)
                %                     xe(ie)=mean(xo(conn(ie,1:elt(ie))));
                %                     ye(ie)=mean(yo(conn(ie,1:elt(ie))));
                %                 end
                %
                %                 masko=GetSignedDistanceToZone(handles,xe(:),ye(:));
                %                 masko=double(masko+mean(handles.fem_model.mesh_size)<0);
                % masko=0.;
                %                                 for iz=1:size(handles.fem_model.zone,2)
                %                                     if handles.fem_model.zone{1,iz}
                %                                         zone=handles.fem_model.zone{2,iz};
                %                                         masko=masko|(inpolygon(xe,ye,zone(:,1),zone(:,2)));
                %                                     end
                %                                 end
                %                                 if numel(masko)==1,masko=1;end
                %                                 for iz=1:size(handles.fem_model.zone,2)
                %                                     if ~handles.fem_model.zone{1,iz}
                %                                         zone=handles.fem_model.zone{2,iz};
                %                                         masko=masko&(~inpolygon(xe,ye,zone(:,1),zone(:,2)));
                %                                     end
                %                                 end
                masko=0;
                for iz=1:size(handles.fem_model.zone,2)
                    if handles.fem_model.zone{1,iz}>0
                        zone=handles.fem_model.zone{2,iz};
                        masko=masko|(inpolygon(xo,yo,zone(:,1),zone(:,2)));
                    end
                end
                if numel(masko)==1,masko=1;end
                for iz=1:size(handles.fem_model.zone,2)
                    if handles.fem_model.zone{1,iz}==0
                        zone=handles.fem_model.zone{2,iz};
                        masko=masko&(~inpolygon(xo,yo,zone(:,1),zone(:,2)));
                    end
                end
                
                mask=mask&masko;
                
            elseif (handles.fem_model.mesh_type==3)
                xo=handles.mvisu.xo;
                yo=handles.mvisu.yo;
                mask=0;
                for iz=1:size(handles.fem_model.zone,2)
                    if handles.fem_model.zone{1,iz}>0
                        zone=handles.fem_model.zone{2,iz};
                        mask=mask|(inpolygon(xo,yo,zone(:,1),zone(:,2)));
                    end
                end
                if numel(mask)==1,mask=1;end
                for iz=1:size(handles.fem_model.zone,2)
                    if handles.fem_model.zone{1,iz}==0
                        zone=handles.fem_model.zone{2,iz};
                        mask=mask&(~inpolygon(xo,yo,zone(:,1),zone(:,2)));
                    end
                end
                
            else
                mask=1;
            end
        case {'none','uni'}
            
            if handles.animation.iim&&handles.preview
                
                
                mask=(handles.eton*handles.mvisu.entropy(:)>=handles.mvisu.smin);
            else
                mask=zeros(prod(handles.mvisu.Nnodes),1);
            end
            
    end
end


function [Uxy]=GetUxy(handles)
iim=handles.animation.iim;
if handles.ana==1
    Ux=handles.upreview((1:prod(handles.mpreview.Nnodes)),iim);
    Uy=handles.upreview(prod(handles.mpreview.Nnodes)+(1:prod(handles.mpreview.Nnodes)),iim);
    coords.xi=handles.mvisu.xo(:);
    coords.yi=handles.mvisu.yo(:);
    Uxy=interpMesh(handles.mpreview,[Ux,Uy],coords);
    Uxy=[Uxy,zeros(prod(handles.mvisu.Nnodes),1)];
else
    if iscell(handles.uvisu)
        Uxy=[reshape(handles.uvisu{iim},prod(handles.mvisu.Nnodes),2),zeros(prod(handles.mvisu.Nnodes),1)];
    else
        Uxy=handles.uvisu((1:prod(handles.mvisu.Nnodes)),iim);
        Uxy=[Uxy,handles.uvisu(prod(handles.mvisu.Nnodes)+(1:prod(handles.mvisu.Nnodes)),iim)];
        Uxy=[Uxy,zeros(prod(handles.mvisu.Nnodes),1)];
        if handles.stereo
            for iz=0:2
                Uxy=[Uxy,handles.uxyz(iz*prod(handles.mvisu.Nnodes)+(1:prod(handles.mvisu.Nnodes)),iim)];
            end
        end
        if handles.param.thermo
            Uxy=[Uxy,handles.uvisu(2*prod(handles.mvisu.Nnodes)+(1:prod(handles.mvisu.Nnodes)),iim)];
        end
    end
end
function handles=ComputeStrainGage(handles,gage,id)
if handles.stereo
    U=handles.uxyz;
else
    U=handles.uvisu;
end
roi=handles.param.roi;

zone= inpolygon(handles.mvisu.xo+roi(1)-1,handles.mvisu.yo+roi(3)-1,gage(:,1),gage(:,2));
if isempty(zone)
    zone=abs(handles.mvisu.xo+roi(1)-1+1i*(handles.mvisu.yo+roi(3)-1)-mean(gage*[1;1i]));
    zone=zone==min(zone);
end
zone=find(zone);
if ~exist(fullfile('TMP','ufreckles.res'),'file')
    xo=handles.mvisu.xo;
    yo=handles.mvisu.yo;
    Nnodes=handles.mvisu.Nnodes;
    Nelems=handles.mvisu.Nelems;
    conn=handles.mvisu.conn;
    if size(conn,2)==3
        conn=[conn,zeros(size(conn,1),1)];
    end
    elt=handles.mvisu.elt;
    rint=handles.sonelt;
    conn=max(conn,1);
    save(fullfile('TMP','ufreckles.res'),'xo','yo','conn','elt','Nnodes','Nelems','rint');
    if handles.stereo
        Xo=handles.mvisu.Xo;
        Yo=handles.mvisu.Yo;
        Zo=handles.mvisu.Zo;
        save(fullfile('TMP','ufreckles.res'),'Xo','Yo','Zo','-append');
    end
end
if handles.stereo
    [dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis25D(fullfile('TMP','ufreckles.res'),handles.sizeim,1,zone,'Gauss_points',true);
else
    [dphidx,dphidy]=CreateGradFiniteElementBasis(fullfile('TMP','ufreckles.res'),handles.sizeim,1,zone,'Gauss_points',true);
end
load(fullfile('TMP','ufreckles.res'),'-mat','Nnodes')
ng=sum(sum(abs(dphidx),2)>0);
Exx=sum(dphidx*U((1:prod(Nnodes)),:),1)/ng;
Eyy=sum(dphidy*U(prod(Nnodes)+(1:prod(Nnodes)),:),1)/ng;
dudy=sum(dphidy*U((1:prod(Nnodes)),:),1)/ng;
dvdx=sum(dphidx*U(prod(Nnodes)+(1:prod(Nnodes)),:),1)/ng;
if handles.fstrain==0
    Exy=0.5*(dudy+dvdx);
    E=[Exx;Eyy;Exy];
else
    Exx=Exx+1;Eyy=Eyy+1;
    FTFxx=Exx.*Exx+dvdx.*dvdx;
    FTFxy=Exx.*dudy+dvdx.*Eyy;
    FTFyy=dudy.*dudy+Eyy.*Eyy;
    E=0.5*[(FTFxx-1);(FTFyy-1);FTFxy];
end
if handles.stereo
    Ezx=sum(dphidx*U(2*prod(Nnodes)+(1:prod(Nnodes)),:),1)/ng;
    Ezy=sum(dphidy*U(2*prod(Nnodes)+(1:prod(Nnodes)),:),1)/ng;
    E=[E;Ezx;Ezy];
end
handles.fem_model.zone{7,id}=E;




function handles=ComputeStrain(handles,restart)
if nargin<2,restart=0;end
if ~exist(fullfile('TMP','ufreckles.res'),'file')
    xo=handles.mvisu.xo;
    yo=handles.mvisu.yo;
    Nnodes=handles.mvisu.Nnodes;
    Nelems=handles.mvisu.Nelems;
    conn=handles.mvisu.conn;
    if size(conn,2)==3
        conn=[conn,zeros(size(conn,1),1)];
    end
    elt=handles.mvisu.elt;
    %    rint=1;
    rint=handles.sonelt;
    conn=max(conn,1);
    save(fullfile('TMP','ufreckles.res'),'xo','yo','conn','elt','Nnodes','Nelems','rint');
    if handles.stereo
        Xo=handles.mvisu.Xo;
        Yo=handles.mvisu.Yo;
        Zo=handles.mvisu.Zo;
        save(fullfile('TMP','ufreckles.res'),'Xo','Yo','Zo','-append');
    end
    restart=1;
end
persistent ijm dphidx dphidy eton gp2cell fstrain
go=0;
if isempty(ijm),go=1;ijm=-1;fstrain=0;end
if handles.animation.iim>0
    if handles.stereo
        U=handles.uxyz(:,handles.animation.iim);
    else
        U=handles.uvisu;
        if iscell(U)
            go=1;
            U=U{handles.animation.iim};
        else
            U=U(:,handles.animation.iim);
        end
    end
    if go||restart
        if handles.stereo
            [dphidx,dphidy,dphidz]=CreateGradFiniteElementBasis25D(fullfile('TMP','ufreckles.res'),handles.sizeim,1,[],'Gauss_points');
        else
            [dphidx,dphidy]=CreateGradFiniteElementBasis(fullfile('TMP','ufreckles.res'),handles.sizeim,1,[],'Gauss_points');
        end
        load(fullfile('TMP','ufreckles.res'),'-mat','xo','yo','rint','elt','Nnodes','Nelems','conn')
        if handles.sonelt
            gp2cell=1;
            eton=1;
            clear MedianFilterCell
        else
            
            val=double(conn>0);
            conn=max(conn,1);
            ie=repmat((1:length(elt))',1,size(conn,2));
            eton=sparse(conn,ie,val,prod(Nnodes),prod(Nelems));
            eton=diag(sparse(1./sum(eton,2)))*eton;
            ngq=4;
            if rint,ngq=1;end
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
            gp2cell=sparse(indi,indj,val);
        end
        
        handles.evisu.xg=zeros(numel(elt),1);
        handles.evisu.yg=zeros(numel(elt),1);
        handles.evisu.xg(elt==3)=mean(xo(conn(elt==3,1:3)),2);
        handles.evisu.xg(elt==4)=mean(xo(conn(elt==4,1:4)),2);
        handles.evisu.yg(elt==3)=mean(yo(conn(elt==3,1:3)),2);
        handles.evisu.yg(elt==4)=mean(yo(conn(elt==4,1:4)),2);
    end
    if  (~(ijm==handles.animation.iim))||(~(fstrain==handles.fstrain))||restart
        load(fullfile('TMP','ufreckles.res'),'-mat','Nnodes')
        ijm=handles.animation.iim;
        fstrain=handles.fstrain;
        if handles.fstrain==1
            Fxx=1+eton*(gp2cell*(dphidx*U((1:prod(Nnodes)),:)));
            Fyy=1+eton*(gp2cell*(dphidy*U(prod(Nnodes)+(1:prod(Nnodes)),:)));
            Fxy=eton*(gp2cell*(dphidy*U((1:prod(Nnodes)),:)));
            Fyx=eton*(gp2cell*(dphidx*U(prod(Nnodes)+(1:prod(Nnodes)),:)));
            tmp=0.5*(Fxx.*Fxx+Fyx.*Fyx-1);
            if handles.medcell&&handles.sonelt
                tmp=MedianFilterCell(handles.mvisu.conn,tmp);
            end
            handles.evisu.xx=tmp;
            tmp=0.5*(Fxx.*Fxy+Fyx.*Fyy);
            if handles.medcell&&handles.sonelt
                tmp=MedianFilterCell(handles.mvisu.conn,tmp);
            end
            handles.evisu.xy=tmp;
            handles.evisu.yx=handles.evisu.xy;
            tmp=0.5*(Fxy.*Fxy+Fyy.*Fyy-1);
            if handles.medcell&&handles.sonelt
                tmp=MedianFilterCell(handles.mvisu.conn,tmp);
            end
            handles.evisu.yy=tmp;
        else
            strain=eton*(gp2cell*(dphidx*U((1:prod(Nnodes)),:)));
            if handles.medcell&&handles.sonelt
                strain=MedianFilterCell(handles.mvisu.conn,strain);
            end
            handles.evisu.xx=strain;
            strain=eton*(gp2cell*(dphidy*U((1:prod(Nnodes)),:)));
            if handles.medcell&&handles.sonelt
                strain=MedianFilterCell(handles.mvisu.conn,strain);
            end
            handles.evisu.xy=strain;
            strain=eton*(gp2cell*(dphidx*U(prod(Nnodes)+(1:prod(Nnodes)),:)));
            if handles.medcell&&handles.sonelt
                strain=MedianFilterCell(handles.mvisu.conn,strain);
            end
            handles.evisu.yx=strain;
            strain=eton*(gp2cell*(dphidy*U(prod(Nnodes)+(1:prod(Nnodes)),:)));
            if handles.medcell&&handles.sonelt
                strain=MedianFilterCell(handles.mvisu.conn,strain);
            end
            handles.evisu.yy=strain;
        end
    end
end
function [xo,yo,conn,elt,Uxy]=MaskMesh(handles)
Uxy=[0,0,0,0,0,0];
getdisp=(nargout>4)&&handles.preview&&handles.animation.iim;
if getdisp, [Uxy]=GetUxy(handles);end
xo=handles.mvisu.xo;
yo=handles.mvisu.yo;
Nnodes=handles.mvisu.Nnodes;
conn=handles.mvisu.conn;
elt=handles.mvisu.elt;
mask=GetMask(handles);
if numel(mask)>1
    mask=find(mask);
    ielt=GetEltsFromNodes(conn,elt,mask,1);
    elt=elt(ielt);
    conn=conn(ielt,:);
    xo=xo(mask);
    yo=yo(mask);
    newids=zeros(prod(Nnodes)+1,1);
    newids(mask)=1:length(mask);
    conn(conn==0)=prod(Nnodes)+1;
    conn=newids(conn);
    
    
    keep=zeros(length(xo),1);
    for in=1:length(xo)
        if ~any(conn(:)==(in))
            keep(in)=1;
        end
    end
    keep=~keep;
    xo=xo(keep);
    yo=yo(keep);
    newids=zeros(length(keep)+1,1);
    newids(keep)=1:sum(keep);
    conn(conn==0)=length(keep)+1;
    conn=newids(conn);
    mask=mask(keep);
    if getdisp
        Uxy=Uxy(mask,:);
    end
end

function handles=plot_mesh(handles)

set(0,'CurrentFigure',handles.figure1);

set(handles.cmenu_disp_data_mag,'Checked','off')
set(handles.cmenu_disp_data_ux,'Checked','off')
set(handles.cmenu_disp_data_uy,'Checked','off')
set(handles.cmenu_disp_data_uz,'Checked','off')
set(handles.cmenu_strain_data_mag,'Checked','off')
set(handles.cmenu_strain_data_exx,'Checked','off')
set(handles.cmenu_strain_data_eyy,'Checked','off')
set(handles.cmenu_strain_data_exy,'Checked','off')
set(handles.cmenu_strain_data_emax,'Checked','off')
set(handles.cmenu_strain_data_emin,'Checked','off')
set(handles.cmenu_strain_data_tmax,'Checked','off')
set(handles.cmenu_field_ondefimage,'Checked','off')
set(handles.cmenu_field_data_axialstrain,'Checked','off')
set(handles.cmenu_field_data_shearstrain,'Checked','off')
set(handles.cmenu_field_data_error,'Checked','off')
set(handles.cmenu_erroronelt_elt,'Checked','off')
set(handles.cmenu_erroronelt_pixel,'Checked','off')
set(handles.cmenu_field_data_topo,'Checked','off')
set(handles.cmenu_field_data_temp,'Checked','off')
try
    for iz=1:size(handles.calculator,2)
        set((handles.calculator{3,iz}),'Checked','off');
    end
catch
end
buf=[sprintf('Step %d/%d',handles.animation.iim,handles.animation.nbstep)];

sizeim=handles.sizeim(1:2);
try delete(handles.gmesh); catch, end
try delete(handles.gedge); catch, end
handles.gmesh=0;
roi=handles.param.roi;
Uxoo=0;Uyoo=0;
Uxc=0;Uyc=0;
if isfield(handles.param,'first_module')
    load(handles.param.result_file,'-mat','Uc')
end
[az,el] = view;
dview=~((az==0)&&(el==90));

[xo,yo,conn,elt,Uxyz]=MaskMesh(handles);
zo=0*xo;
iim=0;nU=0;
decu=3*handles.stereo*handles.pframe*dview;
flatinter='inter';
if handles.preview
    iim=handles.animation.iim;
    if iim||(handles.field==4)
        if isfield(handles.param,'first_module')
            load(handles.param.result_file,'-mat','Uc')
            Uxc=Uc(1:numel(xo),iim);
            Uyc=Uc(numel(xo)+(1:numel(xo)),iim);
        end
        
        if isfield(handles.fem_model,'cn')
            conn=conn(handles.fem_model.cn(:,iim)==0,:);
        end
        buf=[buf,' | '];
        switch handles.field
            case 1
                comp=handles.ucomp;
                if handles.rmrbm
                    if handles.stereo
                        Xo=handles.mvisu.Xo;
                        Yo=handles.mvisu.Yo;
                        Zo=handles.mvisu.Zo;
                        L=[1+0*Xo,0*Xo,0*Xo,-Yo,Zo,0*Xo;...
                            0*Yo,1+0*Yo,0*Yo,Xo,0*Yo,-Zo;...
                            0*Zo,0*Zo,1+0*Zo,0*Zo,-Xo,Yo];
                        Uxo=L*(L\[Uxyz(:,4);Uxyz(:,5);Uxyz(:,6)]);
                        Uyo=Uxo(length(xo)+(1:length(xo)));
                        Uzo=Uxo(2*length(xo)+(1:length(xo)));
                        Uxo=Uxo((1:length(xo)));
                    else
                        L=[1+0*xo,0*xo,-yo;...
                            0*yo,1+0*yo,xo];
                        Uxo=L*(L\[Uxyz(:,1);Uxyz(:,2)]);
                        Uyo=Uxo(length(xo)+(1:length(xo)));
                        Uxo=Uxo((1:length(xo)));
                        Uzo=0;
                        
                        %                         filk=sprintf('%s-crack-%02d-sif.res',strrep(handles.param.result_file,'.res',''),1);
                        %                         Step=handles.fem_model.zone{5,1}.steps;
                        %                         Urbt=0*Step;
                        %                         load(filk,'-mat','xytips','Urbt')
                        %                         Uxo=real(Urbt(Step==iim));
                        %                         Uyo=imag(Urbt(Step==iim));
                        
                        if ~(handles.wfac==1)||~(handles.showim)
                            Uxyz(:,1)=Uxyz(:,1)-Uxo;
                            Uxyz(:,2)=Uxyz(:,2)-Uyo;
                            Uxoo=Uxo;Uyoo=Uyo;
                            Uxo=0;Uyo=0;
                        end
                    end
                    buf=[buf,'wo RBM '];
                else
                    Uxo=0;Uyo=0;Uzo=0;
                end
                if handles.vref==0
                    Uxo=handles.uvisu((1:prod(handles.mvisu.Nnodes)),1);
                    Uyo=handles.uvisu(prod(handles.mvisu.Nnodes)+(1:prod(handles.mvisu.Nnodes)),1);
                    Uzo=0;
                    if handles.rmrbm
                        Uxyo=L*(L\[Uxo;Uyo]);
                        Uxo=Uxo-Uxyo(length(xo)+(1:length(xo)));
                        Uyo=Uyo-Uxyo((1:length(xo)));
                    end
                end
                switch handles.ucomp
                    case 0
                        nU=((Uxyz(:,1+3*handles.stereo)-Uxo).^2+(Uxyz(:,2+3*handles.stereo)-Uyo).^2);
                        if handles.stereo
                            nU=nU+(Uxyz(:,3+3*handles.stereo)-Uzo).^2;
                        end
                        nU=sqrt(nU);
                        set(handles.cmenu_disp_data_mag,'Checked','on')
                        buf=[buf,'displacement magnitude in '];
                    case 1
                        set(handles.cmenu_disp_data_ux,'Checked','on')
                        nU=Uxyz(:,1+3*handles.stereo)-Uxo;
                        buf=[buf,'Ux displacement in '];
                    case 2
                        set(handles.cmenu_disp_data_uy,'Checked','on')
                        nU=Uxyz(:,2+3*handles.stereo)-Uyo;
                        buf=[buf,'Uy displacement in '];
                    case 3
                        set(handles.cmenu_disp_data_uz,'Checked','on')
                        nU=Uxyz(:,3+3*handles.stereo)-Uzo;
                        buf=[buf,'Uz displacement in '];
                end
                nU=nU*handles.param.pixel_size;
                if handles.param.pixel_size==1&&~(handles.stereo)
                    buf=[buf,'pixel'];
                else
                    if max(abs(nU))<1e-3
                        nU=nU*1.e6;
                        buf=[buf,'micron'];
                    elseif max(abs(nU))<1e0
                        nU=nU*1.e3;
                        buf=[buf,'mm'];
                    else
                        buf=[buf,'m'];
                    end
                end
            case 2
                if handles.sonelt
                    flatinter='flat';
                end
                handles=ComputeStrain(handles);
                dudx=handles.evisu.xx;
                dudy=handles.evisu.xy;
                dvdy=handles.evisu.yy;
                dvdx=handles.evisu.yx;
                %                 dudx=handles.evisu.xx(:,(iim-1)*(~iscell(handles.uvisu))+1);
                %                 dudy=handles.evisu.xy(:,(iim-1)*(~iscell(handles.uvisu))+1);
                %                 dvdy=handles.evisu.yy(:,(iim-1)*(~iscell(handles.uvisu))+1);
                %                 dvdx=handles.evisu.yx(:,(iim-1)*(~iscell(handles.uvisu))+1);
                comp=handles.ecomp;
                fbuf=[];
                if handles.fstrain==1,fbuf='GL ';end
                switch handles.ecomp
                    case 0
                        set(handles.cmenu_strain_data_mag,'Checked','on')
                        nU=sqrt(dudx.^2+dvdy.^2+0.5*(dudy+dvdx).^2);
                        buf=[buf,fbuf,'Strain magnitude in %'];
                    case 1
                        set(handles.cmenu_strain_data_exx,'Checked','on')
                        nU=dudx;
                        buf=[buf,'xx ',fbuf,'strain in %'];
                    case 2
                        set(handles.cmenu_strain_data_eyy,'Checked','on')
                        nU=dvdy;
                        buf=[buf,'yy ',fbuf,'strain in %'];
                    case 3
                        set(handles.cmenu_strain_data_exy,'Checked','on')
                        nU=0.5*(dudy+dvdx);
                        buf=[buf,'xy ',fbuf,'strain in %'];
                    case 4
                        set(handles.cmenu_field_data_axialstrain,'Checked','on')
                        nU=handles.beam_model.Eax(:,iim);
                        buf=[buf,'axial ',fbuf,'strain in %'];
                    case 5
                        set(handles.cmenu_strain_data_emax,'Checked','on')
                        nU=0.5*(dudx+dvdy+sqrt((dudx-dvdy).^2+(dudy+dvdx).^2));
                        buf=[buf,'max ',fbuf,'strain in %'];
                    case 6
                        set(handles.cmenu_strain_data_emin,'Checked','on')
                        nU=0.5*(dudx+dvdy-sqrt((dudx-dvdy).^2+(dudy+dvdx).^2));
                        buf=[buf,'min ',fbuf,'strain in %'];
                    case 7
                        set(handles.cmenu_strain_data_tmax,'Checked','on')
                        nU=sqrt((dudx-dvdy).^2+(dudy+dvdx).^2);
                        buf=[buf,'max ',fbuf,'shear in %'];
                    case 8
                        set(handles.cmenu_field_data_shearstrain,'Checked','on')
                        nU=handles.beam_model.Esh(:,iim);
                        buf=[buf,'shear ',fbuf,'strain in %'];
                end
                nU=100*nU;
            case 3
                comp=handles.ercomp;
                if handles.cracked==2
                    erroronelt=0;
                else
                    buf=[buf,'correlation error in %'];
                    if strcmp(handles.fbasis,'uni')
                        fide=fopen(strrep(handles.param.result_file,'.res','-01-error.res'),'r');
                    else
                        fide=fopen(strrep(handles.param.result_file,'.res','-error.res'),'r');
                    end
                    erroronelt=fread(fide,1);
                    dynamic=fread(fide,1);
                end
                if isfield(handles.fem_model,'degree')
                    erroronelt=1;
                end
                if erroronelt
                    switch comp
                        case 0
                            set(handles.cmenu_erroronelt_pixel,'Checked','on')
                            zone=handles.erzone{2,1};
                            zone=round([min(zone(:,1))+5,max(zone(:,1))-5,min(zone(:,2))+5,max(zone(:,2))-5]);
                            nU=GetDiscrepancy2D(handles.uvisu,1,zone,iim);
                            sizeim=size(nU);
                        case 1
                            set(handles.cmenu_erroronelt_elt,'Checked','on')
                            fseek(fide,(iim-1+handles.stereo*(iim+1))*prod(handles.mvisu.Nelems)+2,'bof');
                            nU=fread(fide,prod(handles.mvisu.Nelems));
                            fclose(fide);
                            nU=handles.eton*nU;
                            nU=100*nU/dynamic;
                            flatinter='flat';
                    end
                else
                    set(handles.cmenu_field_data_error,'Checked','on')
                    if handles.cracked==2
                        load(handles.param.result_file,'-mat','Ut')
                        Uw=[reshape(Ut{iim},prod(handles.mvisu.Nnodes),2),zeros(prod(handles.mvisu.Nnodes),1)];
                        
                        switch handles.ucomp
                            case 0
                                nU=((Uxyz(:,1+3*handles.stereo)-Uw(:,1+3*handles.stereo)).^2+(Uxyz(:,2+3*handles.stereo)-Uw(:,2+3*handles.stereo)).^2);
                                nU=sqrt(nU);
                                set(handles.cmenu_disp_data_mag,'Checked','on')
                                buf=[buf,'Fit error magnitude in '];
                            case 1
                                set(handles.cmenu_disp_data_ux,'Checked','on')
                                nU=Uxyz(:,1+3*handles.stereo)-Uw(:,1+3*handles.stereo);
                                buf=[buf,'Fit error along x in '];
                            case 2
                                set(handles.cmenu_disp_data_uy,'Checked','on')
                                nU=Uxyz(:,2+3*handles.stereo)-Uw(:,2+3*handles.stereo);
                                buf=[buf,'Fit error along y in '];
                        end
                        nU=nU*handles.param.pixel_size;
                        if handles.param.pixel_size==1&&~(handles.stereo)
                            buf=[buf,'pixel'];
                        else
                            if max(abs(nU))<1e-3
                                nU=nU*1.e6;
                                buf=[buf,'micron'];
                            elseif max(abs(nU))<1e0
                                nU=nU*1.e3;
                                buf=[buf,'mm'];
                            else
                                buf=[buf,'m'];
                            end
                        end
                        
                    else
                        sizeim=[roi(2)-roi(1),roi(4)-roi(3)]+1;
                        zone=roi;
                        if strcmp(handles.fbasis,'uni')
                            gzone=handles.uni_model.zone{2,1};
                            xp=gzone(1:4,1)+roi(1)-1;
                            yp=gzone(1:4,2)+roi(3)-1;
                            zone(1:2) = [(floor(min(xp))-1),(ceil(max(xp))+1)];
                            zone(3:4) = [(floor(min(yp))-1),(ceil(max(yp))+1)];
                            sizeim=[zone(2)-zone(1),zone(4)-zone(3)]+1;
                            
                        end
                        fseek(fide,(iim-1+handles.stereo*(iim+1))*prod(sizeim)+2,'bof');
                        nU=fread(fide,prod(sizeim));
                        fclose(fide);
                        nU=100*nU/dynamic;
                    end
                end
            case 4
                
                set(handles.cmenu_field_data_topo,'Checked','on')
                nU=handles.mvisu.Zo;
                if handles.rmrbm
                    LL=[1+0*nU,handles.mvisu.Xo,handles.mvisu.Yo];
                    AA=LL\nU;
                    nn=[AA(2);AA(3);1];
                    nn=nn/norm(nn);
                    Zp=LL*AA;
                    nU=(nU-Zp)*nn(3);
                    buf=[buf,'Topography / mean plane in '];
                else
                    buf=[buf,'Topography in '];
                end
                if max(abs(nU))<1e-3
                    nU=nU*1.e6;
                    buf=[buf,'micron'];
                elseif max(abs(nU))<1e0
                    nU=nU*1.e3;
                    buf=[buf,'mm'];
                else
                    buf=[buf,'m'];
                end
                
                comp=0;
            case 5
                comp=handles.ccomp;
                if handles.stereo
                    X=handles.mvisu.Xo;
                    Y=handles.mvisu.Yo;
                    Z=handles.mvisu.Zo;
                else
                    X=xo;
                    Y=yo;
                    Z=zo;
                end
                if handles.rmrbm
                    if handles.stereo
                        L=[1+0*X,0*X,0*X,-Y,Z,0*X;...
                            0*Y,1+0*Y,0*Y,X,0*Y,-Z;...
                            0*Z,0*Z,1+0*Z,0*Z,-X,Y];
                        Uxo=L*(L\[Uxyz(:,4);Uxyz(:,5);Uxyz(:,6)]);
                        Uyo=Uxo(length(X)+(1:length(X)));
                        Uzo=Uxo(2*length(X)+(1:length(X)));
                        Uxo=Uxo((1:length(X)));
                    else
                        L=[1+0*X,0*X,-Y;...
                            0*Y,1+0*Y,X];
                        Uxo=L*(L\[Uxyz(:,1);Uxyz(:,2)]);
                        Uyo=Uxo(length(X)+(1:length(X)));
                        Uxo=Uxo((1:length(X)));
                        Uzo=0;
                    end
                    buf=[buf,'wo RBM '];
                else
                    Uxo=0;Uyo=0;Uzo=0;
                end
                
                Ux=Uxyz(:,1+decu)-Uxo;
                Uy=Uxyz(:,2+decu)-Uyo;
                Uz=Uxyz(:,3+decu)-Uzo;
                if handles.param.thermo==1
                    T=Uxyz(:,4);
                end
                handles=ComputeStrain(handles);
                Exx=handles.evisu.xx;Uxx=Exx;
                Eyy=handles.evisu.yy;Uyy=Eyy;
                Uxy=handles.evisu.xy;
                Uyx=handles.evisu.yx;
                Exy=0.5*(Uxy+Uyx);
                form=strrep(handles.calculator{2,comp+1},'*','.*');
                form=strrep(form,'/','./');
                form=strrep(form,'^','.^');
                nU=eval(form);
                if (numel(nU)==numel(Exx))&& handles.sonelt
                    flatinter='flat';
                end
                buf=[buf,handles.calculator{1,comp+1}];
                if ~isempty(handles.calculator{4,comp+1})
                    buf=[buf,[' in ',handles.calculator{4,comp+1}]];
                end
                set(handles.calculator{3,comp+1},'Checked','on');
            case 6
                
                set(handles.cmenu_field_data_temp,'Checked','on');
                nU=Uxyz(:,4);
                buf=[buf,'T elevation in DL'];
                comp=0;
        end
        cticks=(0:8)/8;
        switch handles.scale_mode(handles.field,comp+1)
            case 1
                umin=min(nU(:));umax=max(nU(:));
                if (handles.field==3)&&(umin>0), umin=0;end
            case 2
                mm=handles.minmax{handles.field,comp+1};
                umin=mm(1);
                umax=mm(2);
            case 3
                umin=min(nU(:));umax=max(nU(:));
                mm=handles.ominmax{handles.field,comp+1};
                umin=min(mm(1),umin);
                umax=max(mm(2),umax);
                handles.ominmax{handles.field,comp+1}=[umin,umax];
        end
        ticks=(umax-umin)*cticks+umin;
        nU=255*(nU-umin)/(max(eps,umax-umin));
        if strcmp(handles.fbasis,'uni')&&handles.ana==1,Ux=0;Uy=0;end
        for it=1:length(ticks)
            sticks{it}=sprintf('%3.2f',(ticks(it)));
        end
        if handles.showcb
            set(handles.gcbar,'Visible','on')
            set(handles.gcbar,'Ytick',255*cticks,'YTickLabel',sticks)
        else
            set(handles.gcbar,'Visible','off')
            set(handles.gcbar,'Ytick',255*cticks,'YTickLabel',sticks)
        end
        
    else
        set(handles.gcbar,'Visible','off')
    end
    if dview
        dx=abs(diff(roi(1:2)));
        dy=abs(diff(roi(3:4)));
        dxy=min(dx,dy);
        %    dxy=1;
        if handles.stereo&&handles.pframe
            xo=handles.mvisu.Xo;
            yo=handles.mvisu.Yo;
            zo=handles.mvisu.Zo;
            scale=min((max(xo)-min(xo)),(max(yo)-min(yo)));
            xo=dxy*(xo-min(xo))/scale+min(handles.mvisu.xo);
            yo=dxy*(yo-min(yo))/scale+min(handles.mvisu.yo);
            zo=dxy*(zo-min(zo))/scale;
            Uxyz=Uxyz*diag([1,1,1,dxy/scale,dxy/scale,dxy/scale]);
        else
            if iim||(handles.field==4)
                zo=handles.zfac*0.5*dxy*nU/255;
            end
        end
    end
else
    set(handles.gcbar,'Visible','off')
end
if size(conn,2)==3
    conn=[conn,zeros(numel(elt),1)];
end
if iim|| handles.field==4
    if handles.field==3&&comp==0&&(handles.cracked<2)
        gmesh=imagesc(zone(1),zone(3),reshape(nU,sizeim)');
    else
        if handles.fem_model.mesh_type==3
            conn=[conn;circshift(conn(elt==4,:),[0,2])];
            conn=conn(:,1:3);
        else
            if any(conn(:,4)==0)
                conn=conn(:,1:3);
            end
        end
        %        gmesh=trimesh(conn,min(max(xo+roi(1)-1+Uxc*handles.ondefimage*(handles.wfac==1)+handles.wfac*handles.ondefimage*Uxyz(:,1+decu),1-1000000*(dview||~handles.showim)),sizeim(1)+1000000*(dview||~handles.showim)),...
        %            min(max(yo+roi(3)-1+Uyc*handles.ondefimage*(handles.wfac==1)+handles.wfac*handles.ondefimage*Uxyz(:,2+decu),1-1000000*(dview||~handles.showim)),sizeim(2)+1000000*(dview||~handles.showim)),...
        %            zo+handles.wfac*handles.ondefimage*Uxyz(:,3+decu),nU);
        
        gmesh=patch('Faces',conn,'Vertices',[min(max(xo+roi(1)-1+Uxc*handles.ondefimage*(handles.wfac==1)+handles.wfac*handles.ondefimage*Uxyz(:,1+decu),1-1000000*(dview||~handles.showim)),sizeim(1)+1000000*(dview||~handles.showim)),...
            min(max(yo+roi(3)-1+Uyc*handles.ondefimage*(handles.wfac==1)+handles.wfac*handles.ondefimage*Uxyz(:,2+decu),1-1000000*(dview||~handles.showim)),sizeim(2)+1000000*(dview||~handles.showim)),...
            zo+handles.wfac*handles.ondefimage*Uxyz(:,3+decu)],'FaceVertexCData',nU);
        
        
        %        gmesh=trimesh(conn(:,1:3),xo+roi(1)-1+handles.wfac*handles.ondefimage*Uxyz(:,1+decu),...
        %            yo+roi(3)-1+handles.wfac*handles.ondefimage*Uxyz(:,2+decu),...
        %            zo+handles.wfac*handles.ondefimage*Uxyz(:,3+decu),nU);
        set(gmesh,'EdgeColor','none',...
            'FaceColor',flatinter,...
            'Marker','none','FaceAlpha',0.75);
    end
    handles.gmesh=gmesh;
    if (handles.ana==2)
        set(gmesh,'uicontextmenu',handles.cmenu_mesh);
        %        set(gmesh,'uicontextmenu',handles.cmenu_export_beam);
    else
        if strcmp(handles.fbasis,'fem'),set(gmesh,'uicontextmenu',handles.cmenu_mesh);end
    end
end
if strcmp(handles.fbasis,'fem')&&(iim==0||handles.showedge)&&(~(handles.field==4))%&&(~isfield(handles.fem_model,'degree'))
    if handles.fem_model.mesh_type==3
        if size(conn,2)==3
            conn=[conn,zeros(size(conn,1),1)];
        end
        seg3=reshape(conn(elt==3,[1,2,1,3,2,3])',2,3*sum(elt==3))';
        seg4=reshape(conn(elt==4,[1,2,2,3,3,4,4,1])',2,4*sum(elt==4))';
        
        seg=unique(sort([seg3;seg4],2),'rows');
        
        xo=xo(seg);
        yo=yo(seg);
        zo=zo(seg);
        if size(Uxyz,1)>1
            Ux=Uxyz(:,1+decu);
            Uy=Uxyz(:,2+decu);
            Uz=Uxyz(:,3+decu);
            Ux=Ux(seg);
            Uy=Uy(seg);
            Uz=Uz(seg);
        else
            Ux=0;Uy=0;Uz=0;
        end
        gedge=plot3(xo'+roi(1)-1+handles.wfac*handles.ondefimage*Ux',yo'+roi(3)-1+handles.wfac*handles.ondefimage*Uy',zo'+handles.wfac*handles.ondefimage*Uz','b-');
    else
        try
            set(handles.gmesh,'EdgeColor','b')
            gedge=handles.gmesh;
        catch
            if any(conn(:,4)==0)
                conn=conn(:,1:3);
            end
            
            gedge=trimesh(conn,min(max(xo+roi(1)-1+handles.wfac*handles.ondefimage*Uxyz(:,1+decu),1-1000000*(dview||~handles.showim)),sizeim(1)+1000000*(dview||~handles.showim)),...
                min(max(yo+roi(3)-1+handles.wfac*handles.ondefimage*Uxyz(:,2+decu),1-1000000*(dview||~handles.showim)),sizeim(2)+1000000*(dview||~handles.showim)),...
                zo+handles.wfac*handles.ondefimage*Uxyz(:,3+decu),nU);
            set(gedge,'EdgeColor','b',...
                'FaceColor','none');
            
        end
    end
    handles.gedge=gedge;
    set(gedge,'uicontextmenu',handles.cmenu_mesh);
end
% else
%     seg3=reshape(conn(elt==3,[1,2,1,3,2,3])',2,3*sum(elt==3))';
%     seg4=reshape(conn(elt==4,[1,2,2,3,3,4,4,1])',2,4*sum(elt==4))';
%
%     seg=unique(sort([seg3;seg4],2),'rows');
%
%     xo=xo(seg);
%     yo=yo(seg);
%
%     gmesh=plot(xo'+roi(1)-1,yo'+roi(3)-1,'b-');
%     set(handles.cmenu_mesh_save_mesh,'Visible','on')
% end

if handles.ondefimage
    set(handles.cmenu_field_ondefimage,'Checked','on');
    %    buf=[buf,' on current frame'];
else
    %    buf=[buf,' on initial frame'];
end

set(handles.field_text,'String',buf);
switch handles.fbasis
    case 'fem'
        for iz=1:size(handles.fem_model.zone,2)
            try delete((handles.fem_model.zone{3,iz}));catch end
            zone=handles.fem_model.zone(:,iz);
            xyp=zone{2};
            gz=plot(xyp(:,1),xyp(:,2),'LineWidth',2);
            gz.Tag=num2str(rand(1));
            if zone{4}==6
                xyp=zone{2};
                in=inpolygon(xo+roi(1)-1,yo+roi(3)-1,xyp(:,1),xyp(:,2));
                in=find(in);
                xon=xo(in)+roi(1)-1;
                yon=yo(in)+roi(3)-1;
                xyon=[xon,yon];
                handles.fem_model.zone{5,iz}=xyon;
                set(gz,'Xdata',xyon(:,1),'Ydata',xyon(:,2))
                set(gz,'Color',[1,0,0],'LineWidth',1,'LineStyle','none','MarkerSize',5);
                loads=zone{7};
                if handles.ana==2
                    if any(abs(loads(1:2,:))>0)
                        set(gz,'Marker','^');
                    else
                        set(gz,'Marker','s');
                    end
                    set(gz,'uicontextmenu',handles.cmenu_export_bcs);
                else
                    if any(abs(loads(:,2))>0)
                        set(gz,'Marker','^');
                    else
                        set(gz,'Marker','s');
                    end
                    set(gz,'uicontextmenu',handles.cmenu_bcs_zone,'ButtonDownFcn','');
                end
            else
                switch handles.ana
                    case {1,3}
                        set(gz,'ButtonDownFcn',@zone_adjust)
                        if zone{1}<0
                            if zone{4}==5
                                set(gz,'Color',[1,0,0]);
                                set(gz,'uicontextmenu',handles.cmenu_crack);
                            else
                                set(gz,'Color',[0,0,0.5]);
                                set(gz,'uicontextmenu',handles.cmenu_attractor);
                            end
                        else
                            if zone{1}
                                set(gz,'Color',[0,0.5,0]);
                            else
                                set(gz,'Color',[0.5,0.,0]);
                            end
                            set(gz,'uicontextmenu',handles.cmenu_zone);
                        end
                    case 2
                        if zone{4}==5
                            zsize=zone{8};
                            if abs(zsize(1))>0
                                if zsize(1)<0
                                    if ~iscell(handles.uvisu)||handles.cracked==2
                                        set(gz,'ButtonDownFcn',@zone_adjust)
                                        try
                                            delete(handles.fem_model.zone{11,iz});
                                        catch
                                        end
                                        Step=zone{5}.steps;
                                        if handles.cracked==2
                                            filk=handles.param.result_file;
                                            Step=1:length(Step);
                                        else
                                            filk=sprintf('%s-crack-%02d-sif.res',strrep(handles.param.result_file,'.res',''),iz);
                                        end
                                        Urbt=0*Step;
                                        load(filk,'-mat','xytips','Urbt')
                                        if any(Step==iim)
                                            urbt=handles.wfac*(Urbt(Step==iim))*handles.ondefimage;
                                            if handles.rmrbm
                                                if ~(handles.wfac==1)||~(handles.showim)
                                                    urbt=handles.wfac*(Urbt(Step==iim)-(Uxoo+1i*Uyoo));
                                                end
                                            end
                                            gzt=plot(xytips(Step==iim)+urbt,'rx','LineWidth',2);
                                            gzt.Tag=num2str(rand(1));
                                            gzt.MarkerSize=10;
                                            handles.fem_model.zone{11,iz}=gzt;
                                            
                                        end
                                    end
                                end
                            end
                            set(gz,'Color',[1,0,0]);
                            set(gz,'uicontextmenu',handles.cmenu_export_crack);
                        else
                            set(gz,'ButtonDownFcn',@zone_adjust)
                            set(gz,'Color',[0,0,0.5]);
                            set(gz,'uicontextmenu',handles.cmenu_export_gage)
                        end
                end
            end
            handles.fem_model.zone{3,iz}=gz;
            
        end
    case 'uni'
        couls='rbmycg';
        for iz=1:size(handles.uni_model.zone,2)
            delete((handles.uni_model.zone{3,iz}));
            zone=handles.uni_model.zone(:,iz);
            xyp=zone{2};
            gz=plot(xyp(:,1),xyp(:,2),'LineWidth',2);
            gz.Tag=num2str(rand(1));
            switch handles.ana
                case 1
                    set(gz,'ButtonDownFcn',@zone_adjust)
                    set(gz,'Color',[0.5,0,0.75]);
                case 2
                    set(gz,'Color',couls(iz));
            end
            set(gz,'uicontextmenu',handles.cmenu_gage);
            handles.uni_model.zone{3,iz}=gz;
        end
    case 'beam'
        for iz=1:size(handles.beam_model.zone,2)
            delete((handles.beam_model.zone{3,iz}));
            zone=handles.beam_model.zone(:,iz);
            xyp=zone{2};
            gz=plot(xyp(:,1),xyp(:,2),'LineWidth',2);
            gz.Tag=num2str(rand(1));
            switch handles.ana
                case 1
                    set(gz,'ButtonDownFcn',@zone_adjust)
                    set(gz,'uicontextmenu',handles.cmenu_beam);
                case 2
                    set(gz,'uicontextmenu',handles.cmenu_export_beam);
            end
            set(gz,'Color',[0.75,0,0]);
            handles.beam_model.zone{3,iz}=gz;
        end
end
guidata(handles.figure1,handles)
if handles.ana==2&&strcmp(handles.fbasis,'fem')&&~(handles.animation.playing)&&handles.showplot
    for iz=1:size(handles.fem_model.zone,2)
        zone=handles.fem_model.zone(:,iz);
        switch zone{4}
            case 4
                switch handles.field
                    case 1
                        %0
                        try delete(100+handles.fgage(iz));catch, end
                        try delete(200+handles.fgage(iz));catch, end
                        try delete(300+handles.fgage(iz));catch, end
                        try delete(500+handles.fgage(iz));catch, end
                        try delete(600+handles.fgage(iz));catch, end
                    case 2
                        %100
                        try delete(handles.fgage(iz));catch, end
                        try delete(200+handles.fgage(iz));catch, end
                        try delete(300+handles.fgage(iz));catch, end
                        try delete(500+handles.fgage(iz));catch, end
                        try delete(600+handles.fgage(iz));catch, end
                    case 3
                        %300
                        try delete(handles.fgage(iz));catch, end
                        try delete(100+handles.fgage(iz));catch, end
                        try delete(200+handles.fgage(iz));catch, end
                        try delete(500+handles.fgage(iz));catch, end
                        try delete(600+handles.fgage(iz));catch, end
                    case 4
                        %200
                        try delete(handles.fgage(iz));catch, end
                        try delete(100+handles.fgage(iz));catch, end
                        try delete(300+handles.fgage(iz));catch, end
                        try delete(500+handles.fgage(iz));catch, end
                        try delete(600+handles.fgage(iz));catch, end
                    case 5
                        %500
                        try delete(handles.fgage(iz));catch, end
                        try delete(100+handles.fgage(iz));catch, end
                        try delete(200+handles.fgage(iz));catch, end
                        try delete(300+handles.fgage(iz));catch, end
                        try delete(600+handles.fgage(iz));catch, end
                    case 6
                        %600
                        try delete(handles.fgage(iz));catch, end
                        try delete(100+handles.fgage(iz));catch, end
                        try delete(200+handles.fgage(iz));catch, end
                        try delete(300+handles.fgage(iz));catch, end
                        try delete(500+handles.fgage(iz));catch, end
                end
                plot_gage_data(handles,zone{2},iz);
                handles=guidata(handles.figure1);
            case 5
                if handles.cracked>0
                    plot_crack_data(handles,iz,0,0);
                    handles=guidata(handles.figure1);
                end
                
        end
    end
end
guidata(handles.figure1,handles)
%set(handles.axes1,'Clim',[0,255]);
set(handles.axes1,'ClimMode','manual');
%keyboard
function handles=display_frame(handles)
if handles.ana==0,return;end
set(0,'CurrentFigure',handles.figure1)
iim=handles.animation.iim;
set(handles.axes1,'Clim',[0,255]);
set(handles.axes1,'ClimMode','manual');
if handles.showim
    if handles.ondefimage&&iim>0
        frame=handles.animation.frames(iim);
        if isfield(handles,'reader')
            im00=readim(handles.reader,frame);
        else
            if iscell(handles.param.deformed_image)
                filim=handles.param.deformed_image{1,frame+handles.stereo};
            else
                filim=handles.param.deformed_image;
                assert(frame<2);
            end
            im00=readim(filim);
        end
    else
        frame=1;
        if isfield(handles,'reader')
            im00=readim(handles.reader,1);
        else
            filim=handles.param.reference_image;
            im00=readim(filim);
        end
    end
    if size(im00,3)==1%Cas d'une image NB
        rgb='MN';
    else
        rgb='RGB';
    end
    if size(im00,3)==1%Cas d'une image NB
        im00=repmat(im00,[1,1,3]);
    end
end
if handles.preview||strcmp(handles.fbasis,'fem')
    if isfield(handles.mvisu,'xos')
        if iim==0
            ijm=handles.animation.nbstep;
        else
            ijm=iim;
        end
        xo=handles.mvisu.xos{ijm};
        yo=handles.mvisu.yos{ijm};
        conn=handles.mvisu.conns{ijm};
        if size(conn,2)==3
            conn=[conn,zeros(size(conn,1),1)];
        end
        Nnodes=size(xo);
        Nelems=[size(conn,1),1,1];
        elt=3*ones(prod(Nelems),1);
        handles.mvisu.xo=xo;
        handles.mvisu.yo=yo;
        handles.mvisu.Nnodes=Nnodes;
        handles.mvisu.Nelems=Nelems;
        handles.mvisu.conn=conn;
        handles.mvisu.elt=elt;
        rint=false;
        save(fullfile('TMP','ufreckles.res'),'xo','yo','conn','elt','Nnodes','Nelems','rint');
        if iim ,handles=ComputeStrain(handles);end
    end
    
    handles=plot_mesh(handles);
end

if handles.showim
    if isfield(handles,'reader')
        img_msg=sprintf('%s %dx%d frame %d/%d %d-bit %s',handles.param.reference_image,handles.sizeim(1:2),frame,handles.param.number_of_frames,8*(1+double(max(im00(:))>255)),rgb);
    else
        img_msg=sprintf('%s %dx%d %d-bit %s',filim,handles.sizeim(1:2),8*(1+double(max(im00(:))>255)),rgb);
    end
    set(handles.image_info_text,'String',img_msg);
    im00=permute(im00,[2,1,3]);
    if any(im00(:)>255)
        im00=double(im00);
        im00=255*(im00-min(im00(:)))/(max(im00(:))-min(im00(:)));
        %     if any(im00(:)>1023)
        %         if any(im00(:)>4095)
        %             im00=im00*(255/65535);
        %         else
        %             im00=im00*(255/4095);
        %         end
        %     else
        %         im00=im00*(255/1023);
        %     end
        im00=uint8(floor(im00));
    end
    set(handles.gim,'CData',im00)
else
    set(handles.image_info_text,'String','no image on display');
    set(handles.gim,'Visible','off')
end
handles=set_param(handles);

% --------------------------------------------------------------------
function save_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.ana
    case {1,3}
        filename=handles.param.result_file;
        
        filename=strrep(filename,'.res','');
        if strcmp(filename,'ufreckles')
            pathres='';
            while ~strcmp(pathres,eval('cd'))
                [filename,pathres,FilterIndex]=uiputfile({'*.dat','Ufreckles Input file (*.dat)'},'Save as a new data set');
                if ~FilterIndex, return; end
                if strcmp(handles.param.analysis,'correlation')
                    if~strcmp(pathres,cd)
                        set(handles.message_text,'String','You must write the result in the same folder as the images !')
                        
                    end
                else
                    cd(pathres);
                end
                pathres(end)=[];
            end
            set(handles.message_text,'String','Data file selected.....')
            
            set(handles.figure1,'Name',sprintf('UFreckles: %s',filename))
            filename=strrep(filename,'.dat','');
            
            handles.param.result_file=[filename,'.res'];
            if handles.stereo==1
                %            if isfield(handles.param,'deformed_image')
                %                if size(handles.param.deformed_image,1)>1
                [fil0,path0,FilterIndex]=uigetfile({...
                    '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2;*.raw','Image files';...
                    '*.cal;','Ufreckles calibration files (*.cal)'...
                    },'Calibration file(s)','MultiSelect', 'on');
                switch FilterIndex
                    case 1
                        pini=cd;
                        filc=Ucalib(path0,fil0);
                        filc=fullfile(path0,filc);
                        cd(pini);
                        LoadParameters(handles.param);
                        ReferenceImage();
                    case 2
                        filc=fullfile(path0,fil0);
                end
                load(filc,'-mat','T');
                for icam=1:size(handles.param.deformed_image,1)
                    handles.param.calibration_data{icam}=T(:,icam);
                end
                
                %                end
            end
            guidata(hObject,handles)
        end
        switch handles.fbasis
            case 'fem'
                save_fem_model(handles);
            case 'uni'
                
                param=handles.param;
                param.psample=1;
                param.umaster=handles.figure1.Tag;
                model=handles.uni_model;
                model.nscale=handles.nscale;
                save([filename,'.dat'],'param','model','-v7.3')
            case 'beam'
                param=handles.param;
                param.psample=1;
                param.umaster=handles.figure1.Tag;
                model=handles.beam_model;
                save([filename,'.dat'],'param','model','-v7.3')
            case 'vic'
                param=handles.param;
                param.psample=1;
                param.umaster=handles.figure1.Tag;
                param.sizeim=handles.sizeim;
                model=handles.vic_model;
                model.nscale=handles.nscale;
                
                for iz=1:size(model.zone,2)
                    if model.zone{4,iz}==1
                        if size(model.zone,1)==10
                            point1=model.zone{10,iz};
                        else
                            set(handles.message_text,'String','Pick up the starting point.....')
                            waitforbuttonpress
                            if strcmp(get(handles.figure1,'selectiontype'),'normal')
                                point1=getposition(handles);
                                point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
                                
                            else
                                set(handles.message_text,'String','Save result file.....cancelled')
                                return
                            end
                        end
                        setok=0;
                        while ~setok
                            gzp=plot(point1(1),point1(2),'o','Color',[1,0.5,0.],'MarkerSize',10,'LineStyle','-','LineWidth',2,'ButtonDownFcn',@zone_adjust);
                            set(handles.message_text,'String','Change starting point ? (right click to set position).....')
                            waitforbuttonpress
                            if strcmp(get(handles.figure1,'selectiontype'),'normal')
                                delete(gzp);
                                point1=getposition(handles);
                                point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
                                
                            else
                                setok=1;
                                model.zone{10,iz}=point1;
                                handles.vic_model.zone{10,iz}=point1;
                                guidata(handles.figure1,handles)
                                set(handles.message_text,'String','Starting point position set.....')
                                pause(1)
                                delete(gzp);
                                
                            end
                        end
                    end
                end
                
                
                
                save([filename,'.dat'],'param','model','-v7.3')
        end
        pause(1)
        
        set(handles.message_text,'String',sprintf('Saving %s.....done',[filename,'.dat']))
    case 2
        switch handles.fbasis
            case 'fem'
                model=handles.fem_model;
                visu.ominmax=handles.ominmax;
                visu.minmax=handles.minmax;
                visu.scale_mode=handles.scale_mode;
                visu.field=handles.field;
                visu.ucomp=handles.ucomp;
                visu.ecomp=handles.ecomp;
                visu.ercomp=handles.ercomp;
                visu.ucomp=handles.ucomp;
                visu.rmrbm=handles.rmrbm;
                visu.calculator=handles.calculator;
                visu.showedge=handles.showedge;
                visu.ondefimage=handles.ondefimage;
                visu.iim=handles.animation.iim;
                [az,el] = view;
                visu.view= [az,el];
                model.visu=visu;
                save(handles.param.result_file,'model','-append');
                set(handles.message_text,'String',sprintf('Result file %s.....saved',handles.param.result_file))
        end
        
end
function [handles]=select_nodes(handles)
ok=1;
while ok
    set(handles.message_text,'String','Click opposite corners to add a node set with Neumann conditions.....')
    [handles,ok]=bcs_zone_set(handles);
    
end
set(handles.message_text,'String','Neumann conditions.....done')




function [handles,ok]=bcs_zone_set(handles)
ok=1;
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    
    handles.initial_point=point1;
    guidata(handles.figure1,handles);
    xyp=[point1;point1];
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',1);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_rect);
    
    waitforbuttonpress
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    
    point2=getposition(handles);
    point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
    ly=abs((point1(1)-point2(1)));
    lx=abs((point1(2)-point2(2)));
    if strcmp(get(handles.figure1,'selectiontype'),'normal')&&(~(max(lx,ly)<1))
        delete(gz);
        roi=[max(1,min(handles.sizeim(1),sort(round([point1(1),point2(1)])))),max(1,min(handles.sizeim(2),sort(round([point1(2),point2(2)]))))];
        xyp=[roi([1,2,2,1,1])',roi([3,3,4,4,3])'];
        switch handles.fbasis
            case 'fem'
                roi=handles.param.roi;
                if (handles.fem_model.mesh_type)==1
                    [xo,yo,conn,elt]=MaskMesh(handles);
                else
                    xo=handles.mvisu.xo;
                    yo=handles.mvisu.yo;
                end
                xo=xo+roi(1)-1;
                yo=yo+roi(3)-1;
                in=inpolygon(xo,yo,xyp(:,1),xyp(:,2));
                in=find(in);
                xon=xo(in);
                yon=yo(in);
                gz=plot(xon,yon,'r^','MarkerSize',5,'LineWidth',1);
                gz.Tag=num2str(rand(1));
                handles.fem_model.zone{2,end+1}=xyp;
                handles.fem_model.zone{3,end}=gz;
                handles.fem_model.zone{4,end}=6;
                handles.fem_model.zone{5,end}=[xon,yon];
                handles.fem_model.zone{6,end}=in;
                handles.fem_model.zone{7,end}=[[NaN(2,1),ones(2,1)];0,0];
                set(gz,'uicontextmenu',handles.cmenu_bcs_zone);
                
                handles=set_param(handles);
        end
    else
        delete(gz);
        set(handles.message_text,'String','Zone definition.....cancelled')
        ok=0;
    end
else
    set(handles.message_text,'String','Zone definition.....cancelled')
    ok=0;
end



function save_fem_model(handles)
param=handles.param;
param.umaster=handles.figure1.Tag;
filename=strrep(param.result_file,'.res','');
roi=param.roi;
%save([filename,'.dat'],'param','model','-v7.3')
if handles.fem_model.mesh_type==1
    [xo,yo,conn,elt]=MaskMesh(handles);
    
else
    xo=handles.mvisu.xo;
    yo=handles.mvisu.yo;
    conn=handles.mvisu.conn;
    elt=handles.mvisu.elt;
    
end
if param.regularization_parameter==0
    param.regularization_type='none';
end
if strcmp(param.regularization_type,'equilibrium_gap')
    [handles]=select_nodes(handles);
end
model=handles.fem_model;
model.nscale=handles.nscale;


if ~(handles.fem_model.mesh_type==3)
    model.mesh_file=[filename,'.vtk'];
    model.gluing_parameters{1}={'translate',[0;0;0]};
end
if handles.ana==1
    param.roi=round([max(roi(1),min(xo)-5),min(roi(2),max(xo)+5),max(roi(3),min(yo)-5),min(roi(4),max(yo)+5)]);
end
xo=xo-roi(1)+1;
yo=yo-roi(3)+1;
Nnodes=[length(xo),1,1];
Nelems=[length(elt),1,1];
if size(conn,2)==3
    conn=[conn,zeros(numel(elt),1)];
end
if strcmp(param.regularization_type,'equilibrium_gap')
    selected=ones(Nnodes);
    for iz=1:size(model.zone,2)
        zone=model.zone(:,iz);
        if zone{4}==6
            if any(sum(isnan(zone{7})))
                xyp=zone{2};
                in=inpolygon(xo+roi(1)-1,yo+roi(3)-1,xyp(:,1),xyp(:,2));
                selected(in)=0;
            end
        end
    end
else
    selected=[];
end
save([filename,'.dat'],'param','model','xo','yo','Nnodes','Nelems','conn','elt','selected','-v7.3')
writeVTKmesh([filename,'.dat']);

writeINPFile([filename,'.ufr'],param,model)

function save_mesh(handles)

[file,pp,filterindex]=uiputfile({'*.vtk','VTK mesh file (*.vtk)';'*.inp','Abaqus mesh format (*.inp)'},'Save mesh as...');
[~,file,ext]=fileparts(file);
filename='tmp';
param=handles.param;
model=handles.fem_model;

if handles.fem_model.mesh_type==1
    [xo,yo,conn,elt]=MaskMesh(handles);
    Nnodes=[length(xo),1,1];
else
    xo=handles.mvisu.xo;
    yo=handles.mvisu.yo;
    Nnodes=handles.mvisu.Nnodes;
    conn=handles.mvisu.conn;
    elt=handles.mvisu.elt;
end
Nelems=[length(elt),1,1];

if size(conn,2)==3
    conn=[conn,zeros(numel(elt),1)];
end

save([filename,'.dat'],'param','model','xo','yo','Nnodes','Nelems','conn','elt','-v7.3')
switch filterindex
    case 1
        writeVTKmesh([filename,'.dat']);
        movefile([filename,'.vtk'],fullfile(pp,[file,'.vtk']));
    case 2
        writeINP([filename,'.dat']);
        movefile([filename,'.inp'],fullfile(pp,[file,'.inp']));
end
delete([filename,'.dat']);



% --------------------------------------------------------------------
function save_as_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save_as_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.ana
    case {1,3}
        pathres='';
        while ~strcmp(pathres,cd)
            [filename,pathres,FilterIndex]=uiputfile('*.dat','Save as a new data set');
            if ~FilterIndex, return; end
            if strcmp(handles.param.analysis,'correlation')
                if~strcmp(pathres,cd)
                    set(handles.message_text,'String','You must write the result in the same folder as the images !')
                    
                end
            else
                cd(pathres);
            end
            pathres(end)=[];
        end
        
        set(handles.message_text,'String','Result file selected.....')
        
        set(handles.figure1,'Name',sprintf('UFreckles: %s',filename))
        filename=strrep(filename,'.dat','');
        
        handles.param.result_file=[filename,'.res'];
        if size(handles.param.deformed_image,1)>1
            [fil0,path0,FilterIndex]=uigetfile({...
                '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2','Image files';...
                '*.cal;','Ufreckles calibration files (*.cal)'...
                },'Calibration file(s)','MultiSelect', 'on');
            switch FilterIndex
                case 1
                    pini=cd;
                    filc=Ucalib(path0,fil0);
                    filc=fullfile(path0,filc);
                    cd(pini);
                    LoadParameters(handles.param);
                    ReferenceImage();
                case 2
                    filc=fullfile(path0,fil0);
            end
            load(filc,'-mat','T');
            for icam=1:size(handles.param.deformed_image,1)
                handles.param.calibration_data{icam}=T(:,icam);
            end
            
        end
        guidata(hObject,handles)
        switch handles.fbasis
            case 'fem'
                save_fem_model(handles);
            case 'uni'
                param=handles.param;
                param.psample=1;
                model=handles.uni_model;
                model.nscale=handles.nscale;
                save([filename,'.dat'],'param','model','-v7.3')
            case 'beam'
                param=handles.param;
                param.psample=1;
                
                model=handles.beam_model;
                model.nscale=handles.nscale;
                save([filename,'.dat'],'param','model','-v7.3')
            case 'vic'
                param=handles.param;
                param.psample=1;
                
                model=handles.vic_model;
                model.nscale=handles.nscale;
                param.sizeim=handles.sizeim;
                save([filename,'.dat'],'param','model','-v7.3')
        end
        pause(1)
        
        set(handles.message_text,'String',sprintf('Saving %s.....done',[filename,'.dat']))
    case 2
        if isempty(handles.fem_model.zone)
            save_inp(handles);
        else
            if ~any(cell2mat(handles.fem_model.zone(4,:))==5)
                save_inp(handles);
            else
                load(strrep(handles.param.result_file,'.res','.dat'),'-mat','model');
                cracko=[];
                try cracko=model.zone(:,cell2mat(model.zone(4,:))==5);catch; end
                crack=handles.fem_model.zone(:,cell2mat(handles.fem_model.zone(4,:))==5);
                crack=crack(:,(size(cracko,2)+1):end);
                for ic=1:size(crack,2)
                    crack{8,ic}=0;
                end
                pathres='';
                xfem=menu('Crack model','FEM','X-FEM');
                filres=strrep(handles.param.result_file,'.res','-crack.dat');
                while ~strcmp(pathres,cd)
                    [filename,pathres,FilterIndex]=uiputfile('*.dat','Save as a new data set',filres);
                    if ~FilterIndex, return; end
                    if~strcmp(pathres,cd)
                        set(handles.message_text,'String','You must write the result in the same folder as the images !')
                    end
                    pathres(end)=[];
                end
                
                set(handles.message_text,'String','Result file selected.....')
                load(strrep(handles.param.result_file,'.res','.dat'),'-mat','param','model');
                nl=size(model.zone,1);
                
                model.zone=[model.zone;cell(size(crack,1)-nl,size(model.zone,2))];
                model.zone=[model.zone,crack];
                model.mesh_type=1;
                model.phantom_nodes=xfem-1;
                param.result_file=strrep(filename,'.dat','.res');
                
                save(filename,'param','model','-v7.3');
                handles=reset_param(handles);
                path0=cd;
                loading_mat_file(handles.figure1,handles,filename,path0,2);
            end
        end
end

function save_inp(handles)
[filename,pathres,FilterIndex]=uiputfile('*.inp','Abaqus Input File (*.inp)',strrep(handles.param.result_file,'.res','-fea.inp'));
if FilterIndex
    filename=strrep(filename,'.inp','');
    handles=reference_frame(handles);
    roi=handles.param.roi;
    handles=guidata(handles.figure1);
    xo=handles.mvisu.xo+roi(1)-1;
    yo=handles.mvisu.yo+roi(3)-1;
    Nnodes=handles.mvisu.Nnodes;
    No=Nnodes;
    gpoint=handles.gpoint;
    if isempty(gpoint)
        LoadParameters(handles.param)
        LoadParameters(handles.fem_model,1)
        ReferenceImage();
        LoadMeshes(1);
        LoadMask(1);
        set(handles.message_text,'String','Click and drag to select the Dirichlet boundary nodes.....')
    else
        set(handles.message_text,'String','Click and drag to add Dirichlet boundary nodes.....right click to end selection')
    end
    %selected=ones(Nnodes);
    load(fullfile('TMP','1_mesh_0'),'selected')
    try set(gpoint,'Visible','on');catch, end
    set(0,'CurrentFigure',handles.figure1);
    waitforbuttonpress
    if strcmp(get(handles.figure1,'selectiontype'),'normal')
        point1=get(gca,'CurrentPoint');
        point1=point1(1,1:2);
        while strcmp(get(handles.figure1,'selectiontype'),'normal')
            rbbox;
            point2=get(gca,'CurrentPoint');
            point2=point2(1,1:2);
            ly=abs((point1(1)-point2(1)));
            lx=abs((point1(2)-point2(2)));
            if ~(max(lx,ly)<1)
                zone=[sort(round([point1(1),point2(1)])),sort(round([point1(2),point2(2)]))];
                selected=selected(:)&(~inpolygon(xo,yo,zone([1,2,2,1,1]),zone([3,3,4,4,3])));
                if isempty(gpoint)
                    gpoint=plot(xo(~selected),yo(~selected),'ro','MarkerSize',5,'LineWidth',2);
                else
                    set(gpoint,'Xdata',xo(~selected));
                    set(gpoint,'Ydata',yo(~selected));
                end
                set(handles.message_text,'String','Selecting nodes.....done')
            else
                set(handles.message_text,'String','No node selected.....')
            end
            pause(1.0)
            set(handles.message_text,'String','Click and drag to add Dirichlet boundary nodes.....right click to end selection')
            
            waitforbuttonpress
            point1=get(gca,'CurrentPoint');
            point1=point1(1,1:2);
        end
    else
        if  ~any(~selected(:))
            try delete(gpoint);catch, end
            handles.gpoint=[];
            handles=set_param(handles);
            set(handles.message_text,'String','Nodes selection.....cancelled')
            return
        end
    end
    if isempty(selected)
        try delete(gpoint);catch, end
        set(handles.message_text,'String','Nodes selection.....cancelled')
        return
    else
        set(gpoint,'uicontextmenu',handles.cmenu_selected_point)
        handles.gpoint=gpoint;
        handles=set_param(handles);
        set(handles.message_text,'String','Nodes selection.....done')
    end
    filres=handles.param.result_file;
    if ~strcmp(filres,[filename,'.res'])
        copyfile(handles.param.result_file,[filename,'.res']);
        copyfile(handles.param.result_file,[filename,'.dat']);
    end
    load([filename,'.res'],'-mat','param','model','U');
    param.analysis='mechanics';
    param.result_file=[filename,'.res'];
    set(handles.message_text,'String','FE model definition.....')
    [model,doextrusion]=Uaba(param,model);
    set(handles.message_text,'String','FE model definition.....done')
    
    dflag=isfield(handles.param,'calibration_parameters')||(length(handles.param.roi)==6);
    LoadParameters(param)
    model.material_parameters.nu=0.3;
    model.material_parameters.young=2.e11;
    model.material_parameters.KIc=1.e7;
    model.mesh_file=[filename,'.vtk'];
    
    LoadParameters(model,1);
    save([filename,'.res'],'selected','model','param','xo','yo','-append');
    save([filename,'.dat'],'selected','model','param','xo','yo','-append');
    save(fullfile('TMP','1_mesh_0'),'selected','-append')
    selectedo=selected;
    if model.rmesh
        U=RefineMesh(1,model.rmesh,U);
        load(fullfile('TMP','1_mesh_0'),'selected','xo','yo','Nnodes','Nelems','conn','elt')
        xo=xo+param.roi(1)-1;
        yo=yo+param.roi(3)-1;
        save([filename,'.res'],'U','selected','xo','yo','Nnodes','Nelems','conn','elt','-append')
        save([filename,'.dat'],'U','selected','xo','yo','Nnodes','Nelems','conn','elt','-append')
    end
    if doextrusion
        [U]=ExtrudeModel(1,U,0);
        load(fullfile('TMP','1_3d_mesh_0'),'selected','xo','yo','zo','Nnodes','Nelems','conn','elt')
        xo=xo+param.roi(1)-1;
        yo=yo+param.roi(3)-1;
        save([filename,'.res'],'U','selected','xo','yo','zo','Nnodes','Nelems','conn','elt','-append')
        save([filename,'.dat'],'U','selected','xo','yo','zo','Nnodes','Nelems','conn','elt','-append')
    end
    
    Nddl=(2+(dflag||doextrusion))*prod(Nnodes);
    Fos=zeros(Nddl,1);
    indio=find(~selected(:));
    indi=[indio;indio+prod(Nnodes)];
    if dflag||doextrusion
        if dflag
            indi=[indi;indio+2*prod(Nnodes)];
        else
            if strcmp(model.extrusion_parameters.type,'3d')
                indi=[indi;2*prod(Nnodes)+((prod(Nnodes)-prod(No)+1):prod(Nnodes))'];
            end
        end
        
    end
    Cs=sparse(indi,1:length(indi),1,Nddl,length(indi));
    Ups=U(indi,:);
    switch FilterIndex
        case 1 %ABAQUS
            genINPFile(1,Fos,Cs,Ups)
    end
    if dflag||doextrusion
        writeVTKmesh3D([filename,'.res']);
    else
        if dflag
            writeVTKmesh25D([filename,'.res']);
        else
            writeVTKmesh([filename,'.res']);
        end
    end
    if model.rmesh
        LoadParameters(handles.param)
        LoadParameters(handles.fem_model,1)
        LoadMeshes(1);
        LoadMask(1);
        selected=selectedo;
        save(fullfile('TMP','1_mesh_0'),'selected','-append')
    end
end



% --------------------------------------------------------------------
function preview_button_OffCallback(hObject, eventdata, handles)
% hObject    handle to preview_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(handles.figure1);
handles.preview=0;
handles.ondefimage=1;
handles.showim=1;
if handles.ana>0
    handles=plot_mesh(handles);
    handles=set_param(handles);
end

% --------------------------------------------------------------------
function preview_button_OnCallback(hObject, eventdata, handles)
% hObject    handle to preview_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.preview=1;
switch handles.ana
    case 1
        if isempty(handles.upreview)
            handles=run_preview(handles);
            handles=set_param(handles);
        else
            handles=plot_mesh(handles);
            handles=set_param(handles);
            
        end
    case 2
        handles=plot_mesh(handles);
        handles=set_param(handles);
end

% --------------------------------------------------------------------
function run_analysis_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to run_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Urun;

% --- Executes during object creation, after setting all properties.
function nbs_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbs_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function cmenu_export_gage_export_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_gage_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
gage=handles.fem_model.zone{2,iz};
plot_gage_data(handles,gage,iz,1)
% --------------------------------------------------------------------
function cmenu_export_gage_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_gage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_field_ondefimage_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_ondefimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ondefimage=~handles.ondefimage;
handles=display_frame(handles);
handles=set_param(handles);
% --------------------------------------------------------------------
function cmenu_field_data_error_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
if isfield(handles.fem_model,'degree')
    if isempty(handles.erzone)
        LoadParameters(handles.param);
        LoadParameters(handles.fem_model,1);
        ReferenceImage(1);
        LoadMeshes(1);
        LoadMask(1);
        zone_er_set(handles);
    else
        [az,el] = view;
        dview=(az==0)&&(el==90);
        if ~dview
            handles.view=[az,el];
            view([0,90]);
        end
    end
    handles=guidata(handles.figure1);
    if ~isempty(handles.erzone)
        handles.field=3;
        handles.ercomp=0;
        handles=plot_mesh(handles);
        handles=set_param(handles);
    end
else
    handles.field=3;
    handles.ercomp=0;
    handles=reference_frame(handles);
    handles=plot_mesh(handles);
    handles=set_param(handles);
end

% --------------------------------------------------------------------
function cmenu_field_data_erroronelt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_erroronelt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_erroronelt_elt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_erroronelt_elt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=3;
handles.ercomp=1;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_erroronelt_pixel_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_erroronelt_pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.erzone)
    LoadParameters(handles.param);
    LoadParameters(handles.fem_model,1);
    ReferenceImage(1);
    LoadMeshes(1);
    LoadMask(1);
    zone_er_set(handles);
else
    [az,el] = view;
    dview=(az==0)&&(el==90);
    if ~dview
        handles.view=[az,el];
        view([0,90]);
    end
end
handles=guidata(handles.figure1);
handles.field=3;
handles.ercomp=0;
handles=plot_mesh(handles);
handles=set_param(handles);


function cmenu_field_data_topo_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_topo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=4;
handles.animation.iim=0;
handles=display_frame(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_field_data_temp_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=6;
handles=plot_mesh(handles);
handles=set_param(handles);


function cmenu_field_data_disp_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function cmenu_disp_data_mag_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_disp_data_mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=1;
handles.ucomp=0;
handles=plot_mesh(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_disp_data_ux_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_disp_data_ux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=1;
handles.ucomp=1;
handles=plot_mesh(handles);
handles=set_param(handles);



% --------------------------------------------------------------------
function cmenu_disp_data_uy_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_disp_data_uy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=1;
handles.ucomp=2;
handles=plot_mesh(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_disp_data_uz_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_disp_data_uz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=1;
handles.ucomp=3;
handles=plot_mesh(handles);
handles=set_param(handles);




% --------------------------------------------------------------------
function cmenu_field_data_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_mesh_create_gage_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_create_gage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%p_create_gage_Callback(handles.p_create_gage, eventdata, handles)
zone_set(handles)



% --------------------------------------------------------------------
function cmenu_mesh_save_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_save_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_mesh(handles);

% --------------------------------------------------------------------
function cmenu_mesh_save_inp_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_save_inp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_inp(handles)

% --------------------------------------------------------------------
function cmenu_export_gage_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_gage_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.fem_model.zone(:,iz)=[];
try delete(handles.fgage(iz));catch , end
try delete(100+handles.fgage(iz));catch , end
try delete(200+handles.fgage(iz));catch , end
try delete(300+handles.fgage(iz));catch , end
try delete(400+handles.fgage(iz));catch , end
try delete(4000+handles.fgage(iz));catch , end
try delete(500+handles.fgage(iz));catch , end
try delete(600+handles.fgage(iz));catch , end
guidata(hObject,handles)
delete(gco)

% --------------------------------------------------------------------
function cmenu_field_data_strain_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_strain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
function cmenu_field_data_axialstrain_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_axialstrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
handles.field=2;
handles.ecomp=4;
handles=plot_mesh(handles);
handles=set_param(handles);

function cmenu_strain_data_mag_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_strain_data_mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=0;
handles=plot_mesh(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_strain_data_exx_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_strain_data_exx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=1;
handles=plot_mesh(handles);
handles=set_param(handles);



% --------------------------------------------------------------------
function cmenu_strain_data_eyy_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_strain_data_eyy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=2;
handles=plot_mesh(handles);
handles=set_param(handles);



% --------------------------------------------------------------------
function cmenu_strain_data_exy_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_strain_data_exy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=3;
handles=plot_mesh(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_strain_data_emax_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_strain_data_emax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=5;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_strain_data_emin_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_strain_data_emin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=6;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_strain_data_tmax_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_strain_data_tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=7;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_colorbar_auto_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_colorbar_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.field
    case 1
        comp=handles.ucomp;
    case 2
        comp=handles.ecomp;
    case 5
        comp=handles.ccomp;
    case 3
        comp=handles.ercomp;
    case {4,6}
        comp=0;
end
handles.scale_mode(handles.field,comp+1)=1;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_colorbar_minmax_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_colorbar_minmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.field
    case 1
        comp=handles.ucomp;
    case 2
        comp=handles.ecomp;
    case 5
        comp=handles.ccomp;
    case 3
        comp=handles.ercomp;
    case {4,6}
        comp=0;
end
handles.scale_mode(handles.field,comp+1)=3;
handles.ominmax{handles.field,comp+1}=[Inf,-Inf];
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_colorbar_user_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_colorbar_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.field
    case 1
        comp=handles.ucomp;
    case 2
        comp=handles.ecomp;
    case 5
        comp=handles.ccomp;
    case 3
        comp=handles.ercomp;
        
    case {4,6}
        comp=0;
end
handles.scale_mode(handles.field,comp+1)=2;
try
    mm=handles.minmax{handles.field,comp+1};
    if isempty(mm)
        ticks=get(handles.gcbar,'YTickLabel');
        mm=[min(eval(ticks{1})),max(eval(ticks{end}))];
    end
    answer=inputdlg({'Min value:','Max value:',},'Defined min/max field value',1,...
        {num2str(mm(1)),num2str(mm(2))});
    mm(1)=eval(answer{1});
    mm(2)=eval(answer{2});
    
catch
    ticks=get(handles.gcbar,'YTickLabel');
    mm=[min(eval(ticks{1})),max(eval(ticks{end}))];
    answer=inputdlg({'Min value:','Max value:',},'Defined min/max field value',1,...
        {num2str(mm(1)),num2str(mm(2))});
    mm(1)=eval(answer{1});
    mm(2)=eval(answer{2});
    
    
end
handles.minmax{handles.field,comp+1}=mm;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cmenu_colorbar_user,'Checked','off')
set(handles.cmenu_colorbar_minmax,'Checked','off')
set(handles.cmenu_colorbar_auto,'Checked','off')
switch handles.field
    case 1
        comp=handles.ucomp;
    case 2
        comp=handles.ecomp;
    case 5
        comp=handles.ccomp;
    case 3
        comp=handles.ercomp;
    case {4,6}
        comp=0;
end
switch handles.scale_mode(handles.field,comp+1)
    case 1
        set(handles.cmenu_colorbar_auto,'Checked','on')
    case 2
        set(handles.cmenu_colorbar_user,'Checked','on')
    case 3
        set(handles.cmenu_colorbar_minmax,'Checked','on')
end


% --------------------------------------------------------------------
function cmenu_zone_inout_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_zone_inout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
if (handles.fem_model.zone{1,iz})
    handles.fem_model.zone{1,iz}=0;
    set(gco,'Color',[0.5,0,0]);
else
    handles.fem_model.zone{1,iz}=1;
    set(gco,'Color',[0,0.5,0]);
end
if handles.fem_model.mesh_type==2
    handles=remesh(handles);
else
    handles=plot_mesh(handles);
end
% --------------------------------------------------------------------
function cmenu_zone_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_zone_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.fem_model.zone(:,iz)=[];
delete(gco)
if handles.fem_model.mesh_type==2
    handles=remesh(handles);
else
    handles=plot_mesh(handles);
end
% --------------------------------------------------------------------
function cmenu_mesh_zone_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_mesh_zone_rect_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_zone_rect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zone_set(handles)

function plot_gage_data(handles,gage,id,export)
if nargin<4,export=0;end
if handles.preview
    handles.showplot=1;
    handles.fgage(id)=id;
    guidata(handles.figure1,handles);
    if size(gage,1)>2 %cas d'une jauge
        roi=handles.param.roi;
        l1=abs(diff(gage(1:2,:)*[1;1i]));
        l2=abs(diff(gage(2:3,:)*[1;1i]));
        angl=angle(diff(gage(1:2,:)*[1;1i]));
        if l2<l1,angl=angl+pi/2;end
        angl=exp(1i*angl);
        R=[real(angl*exp(1i*pi/2)),real(angl);...
            imag(angl*exp(1i*pi/2)),imag(angl)];
        
        switch handles.ana
            case 1
                zone= inpolygon(handles.mpreview.xo,handles.mpreview.yo,gage(:,1),gage(:,2));
                if isempty(zone)
                    zone=abs(handles.mpreview.xo+1i*handles.mpreview.yo-mean(gage*[1;1i]));
                    zone=zone==min(zone);
                end
                Nnodes=handles.mpreview.Nnodes;
                Exx=mean(handles.epreview.xx(zone,:),1);
                dudy=mean(handles.epreview.xy(zone,:),1);
                dvdx=mean(handles.epreview.yx(zone,:),1);
                Eyy=mean(handles.epreview.yy(zone,:),1);
                Exy=0.5*(dudy+dvdx);
                
            case 2
                switch handles.fbasis
                    case 'uni'
                        Up=handles.uni_model.Up{id};
                        Exx=Up(end-2,:);
                        Eyy=Up(end-1,:);
                        Exy=Up(end,:);
                    case 'fem'
                        %                         zone= inpolygon(handles.mvisu.xo+roi(1)-1,handles.mvisu.yo+roi(3)-1,gage(:,1),gage(:,2));
                        %                         if isempty(zone)
                        %                             zone=abs(handles.mvisu.xo+roi(1)-1+1i*(handles.mvisu.yo+roi(3)-1)-mean(gage*[1;1i]));
                        %                             zone=zone==min(zone);
                        %                         end
                        %                         Nnodes=handles.mvisu.Nnodes;
                        %                         Exx=mean(handles.evisu.xx(zone,:),1);
                        %                         Eyy=mean(handles.evisu.yy(zone,:),1);
                        %                         dudy=mean(handles.evisu.xy(zone,:),1);
                        %                         dvdx=mean(handles.evisu.yx(zone,:),1);
                        %                         Exy=0.5*(dudy+dvdx);
                        Exx=handles.fem_model.zone{7,id};
                        Eyy=Exx(2,:);
                        Exy=Exx(3,:);
                        if handles.stereo
                            Ezx=Exx(4,:);
                            Ezy=Exx(5,:);
                        end
                        Exx=Exx(1,:);
                end
        end
        El=R(1,1)*(R(1,1)*Exx+R(2,1)*Exy)+R(2,1)*(R(1,1)*Exy+R(2,1)*Eyy);
        Et=R(1,2)*(R(1,2)*Exx+R(2,2)*Exy)+R(2,2)*(R(1,2)*Exy+R(2,2)*Eyy);
        Es=R(1,2)*(R(1,1)*Exx+R(2,1)*Exy)+R(2,2)*(R(1,1)*Exy+R(2,1)*Eyy);
        if handles.stereo
            Ewl=Ezx*R(1,1)+Ezy*R(2,1);
            Ewt=Ezx*R(1,2)+Ezy*R(2,2);
        end
        %try delete(id);catch, end
        fgage=figure(id);
        if ~isempty(findobj(gca,'DisplayName','\epsilon_{l}'))
            set(findobj(gca,'DisplayName','\epsilon_{l}'),'Ydata',El);
            set(findobj(gca,'DisplayName','\epsilon_{t}'),'Ydata',Et);
            set(findobj(gca,'DisplayName','\tau'),'Ydata',Es);
        else
            delete(gca)
            axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                'FontName','Times');
            set(fgage,'Name',sprintf('Gage #%d',id));
            set(fgage,'NumberTitle','off');
            box('on');
            hold('all');
            if numel(El)==1
                plot(El,'x','DisplayName','\epsilon_{l}','Parent',axes1,'LineWidth',2);
                plot(Et,'x','DisplayName','\epsilon_{t}','Parent',axes1,'LineWidth',2);
                plot(Es,'x','DisplayName','\tau','Parent',axes1,'LineWidth',2);
            else
                plot(El,'DisplayName','\epsilon_{l}','Parent',axes1,'LineWidth',2);
                plot(Et,'DisplayName','\epsilon_{t}','Parent',axes1,'LineWidth',2);
                plot(Es,'DisplayName','\tau','Parent',axes1,'LineWidth',2);
            end
            ylabel('Strain []','FontSize',20,'FontName','Times');
            xlabel('Step','FontSize',20,'FontName','Times');
            leg=legend(axes1,'show','Location','NorthWest');
            set(leg,'Box','off')
            xlim([0,(numel(El)+1)])
        end
        if export
            [pp,filename,ext]=fileparts(handles.param.result_file);
            filexp=sprintf('%s-gage-%02d.csv',filename,id);
            if exist(filexp,'file')
                set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save gage data...',filexp);
            end
            
            fid=fopen(filexp,'w');
            fprintf(fid,'Result file;%s\n',handles.param.result_file);
            if ~strcmp(handles.param.analysis,'mechanics')
                fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
            end
            fprintf(fid,'Strain gage;%d\n',id);
            if strcmp(handles.fbasis,'fem')
                if handles.fstrain==0
                    fprintf(fid,'Small strain\n');
                else
                    
                    fprintf(fid,'Finite strain\n');
                end
            end
            fprintf(fid,'Position X [pixel];Position Y [pixel]; Orientation; Length [pixel]; Width [pixel]\n');
            fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1)),mean(gage(1:4,2)),(angle(angl)+pi/2)*180/pi,max(l1,l2),min(l1,l2));
            if ~(handles.param.pixel_size==1)
                pix2m=handles.param.pixel_size;
                fprintf(fid,'Position X [m];Position Y [m]; Orientation; Length [m]; Width [m]\n');
                fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1))*pix2m,mean(gage(1:4,2))*pix2m,(angle(angl)+pi/2)*180/pi,max(l1,l2)*pix2m,min(l1,l2)*pix2m);
            end
            fprintf(fid,'Filename;Step;Longitudinal strain [];Transverse strain [];Shear strain []\n');
            for iim=1:length(El)
                if iscell(handles.param.deformed_image)
                    fprintf(fid,'"%s";%d;%.3e;%.3e;%.3e\n',handles.param.deformed_image{1,iim+handles.stereo},iim,El(iim),Et(iim),Es(iim));
                else
                    fprintf(fid,'"%s";%d;%.3e;%.3e;%.3e\n',handles.param.deformed_image,iim,El(iim),Et(iim),Es(iim));
                end
            end
            fclose(fid);
            
            
            %    dlmwrite(filexp,{'Step;Longitudinal strain;Transverse strain;Shear strain'},'delimiter','');
            %   dlmwrite(filexp,[(1:length(El))',El',Et',Es'],'delimiter',';','-append');
            param=handles.param;
            model=handles.fem_model;
            save(sprintf('%s-gage-%02d.res',filename,id),'El','Et','Es','gage','param','model','-v7.3')
            if handles.stereo
                save(sprintf('%s-gage-%02d.res',filename,id),'Ewl','Ewt','-append')
            end
            set(handles.message_text,'String','Exporting gage data.....done')
            
            
        end
    else %cas d'une ligne
        roi=handles.param.roi;
        t=(diff(gage(1:2,:)*[1;1i]));
        l1=abs(t);
        t=t/l1;
        angl=angle(diff(gage(1:2,:)*[1;1i]));
        angl=exp(1i*angl);
        R=[real(angl*exp(1i*pi/2)),real(angl);...
            imag(angl*exp(1i*pi/2)),imag(angl)];
        nn=round(l1/max(1,l1/100));
        dl=l1/nn;
        s=dl*(0:nn)';
        xyp=gage(1,:)*[1;1i]+s*t;
        coords.xi=real(xyp)-roi(1)+1;
        coords.yi=imag(xyp)-roi(3)+1;
        
        Nnodes=handles.mvisu.Nnodes;
        iim= handles.animation.iim;
        switch handles.field
            case 1
                if iim
                    if handles.stereo
                        U=handles.uxyz(:,iim);
                        U=[reshape(U,prod(Nnodes),3)];
                    else
                        U=handles.uvisu(:,iim);
                        U=reshape(U(1:2*prod(Nnodes)),prod(Nnodes),2);
                    end
                    if handles.rmrbm
                        if handles.stereo
                            Xo=handles.mvisu.Xo;
                            Yo=handles.mvisu.Yo;
                            Zo=handles.mvisu.Zo;
                            L=[1+0*Xo,0*Xo,0*Xo,-Yo,Zo,0*Xo;...
                                0*Yo,1+0*Yo,0*Yo,Xo,0*Yo,-Zo;...
                                0*Zo,0*Zo,1+0*Zo,0*Zo,-Xo,Yo];
                            Uxo=L*(L\[U(:,1);U(:,2);U(:,3)]);
                            Uyo=Uxo(length(Xo)+(1:length(Xo)));
                            Uzo=Uxo(2*length(Xo)+(1:length(Xo)));
                            Uxo=Uxo((1:length(Xo)));
                            U=U-[Uxo,Uyo,Uzo];
                        else
                            L=[ones(prod(Nnodes),1),zeros(prod(Nnodes),1),-handles.mvisu.yo;...
                                zeros(prod(Nnodes),1),ones(prod(Nnodes),1),handles.mvisu.xo];
                            Uxo=L*(L\[U(:,1);U(:,2)]);
                            Uyo=Uxo(prod(Nnodes)+(1:prod(Nnodes)));
                            Uxo=Uxo((1:prod(Nnodes)));
                            U=U-[Uxo,Uyo];
                        end
                    end
                    
                    Uxyz=interpMesh(handles.mvisu,U,coords,1);
                    
                    
                    
                else
                    Uxyz=zeros(length(s),3);
                end
                
                fgage=figure(id);
                if ~isempty(findobj(gca,'DisplayName','U_x'))
                    set(findobj(gca,'DisplayName','U_x'),'Ydata',Uxyz(:,1));
                    set(findobj(gca,'DisplayName','U_x'),'Xdata',s);
                    set(findobj(gca,'DisplayName','U_y'),'Ydata',Uxyz(:,2));
                    set(findobj(gca,'DisplayName','U_y'),'Xdata',s);
                    if ~isempty(findobj(gca,'DisplayName','U_z'))
                        set(findobj(gca,'DisplayName','U_z'),'Ydata',Uxyz(:,3));
                        set(findobj(gca,'DisplayName','U_z'),'Xdata',s);
                    end
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                else
                    delete(gca)
                    axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                        'FontName','Times');
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                    set(fgage,'NumberTitle','off');
                    box('on');
                    hold('all');
                    plot(s,Uxyz(:,1),'DisplayName','U_x','Parent',axes1,'LineWidth',2);
                    plot(s,Uxyz(:,2),'DisplayName','U_y','Parent',axes1,'LineWidth',2);
                    if handles.stereo
                        plot(s,Uxyz(:,3),'DisplayName','U_z','Parent',axes1,'LineWidth',2);
                        ylabel('Displacement [m]','FontSize',20,'FontName','Times');
                    else
                        ylabel('Displacement [pixel]','FontSize',20,'FontName','Times');
                    end
                    xlabel('s [pixel]','FontSize',20,'FontName','Times');
                    leg=legend(axes1,'show','Location','NorthWest');
                    set(leg,'Box','off')
                end
                if export
                    if handles.stereo
                        U=handles.uxyz;
                    else
                        U=handles.uvisu;
                    end
                    Nnodes=handles.mvisu.Nnodes;
                    Ux=interpMesh(handles.mvisu,U((1:prod(Nnodes)),:),coords);
                    Uy=interpMesh(handles.mvisu,U(prod(Nnodes)+(1:prod(Nnodes)),:),coords);
                    if handles.stereo
                        Uz=interpMesh(handles.mvisu,U(2*prod(Nnodes)+(1:prod(Nnodes)),:),coords);
                    end
                    [pp,filename,ext]=fileparts(handles.param.result_file);
                    filexp=sprintf('%s-line-%02d-disp.csv',filename,id);
                    if exist(filexp,'file')
                        set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                        [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save line data...',filexp);
                    end
                    
                    fid=fopen(filexp,'w');
                    fprintf(fid,'Result file;%s\n',handles.param.result_file);
                    if ~strcmp(handles.param.analysis,'mechanics')
                        fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
                    end
                    fprintf(fid,'Line plot;%d\n',id);
                    fprintf(fid,'X1 [pixel];Y1 [pixel];X2 [pixel];Y2 [pixel]; Orientation; Length [pixel];\n');
                    fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:),gage(2,:),(angle(angl))*180/pi,l1);
                    if ~(handles.param.pixel_size==1)
                        pix2m=handles.param.pixel_size;
                        fprintf(fid,'X1 [m];Y1 [m];X2 [m];Y2 [m]; Orientation; Length [m];\n');
                        fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:)*pix2m,gage(2,:)*pix2m,(angle(angl))*180/pi,l1*pix2m);
                        fprintf(fid,'Position X [m];Position Y [m]; Orientation; Length [m]; Width [m]\n');
                        fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1))*pix2m,mean(gage(1:4,2))*pix2m,(angle(angl)+pi/2)*180/pi,max(l1,l2)*pix2m,min(l1,l2)*pix2m);
                    end
                    if handles.stereo
                        
                        if isfield(handles.param,'deformed_image')
                            fprintf(fid,'Filename;;;');
                            for iim=1:size(Uy,2)
                                if iscell(handles.param.deformed_image)
                                    fprintf(fid,'"%s";;;',handles.param.deformed_image{1,iim+handles.stereo});
                                else
                                    fprintf(fid,'"%s";;;',handles.param.deformed_image);
                                end
                            end
                            fprintf(fid,'\n');
                        end
                        fprintf(fid,'Step;;;');
                        for iim=1:size(Uy,2)
                            fprintf(fid,'%d;;;',iim);
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'x [pixel];y [pixel];s [pixel];Ux[m];Uy[m];Uz[m]\n');
                        for ip=1:length(s)
                            fprintf(fid,'%.3f;%.3f;%.3f;',real(xyp(ip)),imag(xyp(ip)),s(ip));
                            for iim=1:size(Uy,2)
                                fprintf(fid,'%.3e;%.3e;%.3e;',Ux(ip,iim),Uy(ip,iim),Uz(ip,iim));
                            end
                            fprintf(fid,'\n');
                        end
                        
                    else
                        
                        if isfield(handles.param,'deformed_image')
                            fprintf(fid,'Filename;;;');
                            for iim=1:size(Uy,2)
                                if iscell(handles.param.deformed_image)
                                    fprintf(fid,'"%s";;',handles.param.deformed_image{1,iim+handles.stereo});
                                else
                                    fprintf(fid,'"%s";;',handles.param.deformed_image);
                                end
                            end
                            fprintf(fid,'\n');
                        end
                        fprintf(fid,'Step;;;');
                        for iim=1:size(Uy,2)
                            fprintf(fid,'%d;;',iim);
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'x [pixel];y [pixel];s [pixel];Ux[pixel];Uy[pixel];\n');
                        for ip=1:length(s)
                            fprintf(fid,'%.3f;%.3f;%.3f;',real(xyp(ip)),imag(xyp(ip)),s(ip));
                            for iim=1:size(Uy,2)
                                fprintf(fid,'%.3f;%.3f;',Ux(ip,iim),Uy(ip,iim));
                            end
                            fprintf(fid,'\n');
                        end
                    end
                    fclose(fid);
                    
                    
                    param=handles.param;
                    model=handles.fem_model;
                    save(sprintf('%s-line-%02d-disp.res',filename,id),'Ux','Uy','gage','param','model','-v7.3')
                    if handles.stereo
                        save(sprintf('%s-line-%02d-disp.res',filename,id),'Uz','-append')
                    end
                    set(handles.message_text,'String','Exporting line data.....done')
                    
                    
                end
            case 2
                if iim
                    %                    E=interpMesh(handles.mvisu,[handles.evisu.xx(:,iim),handles.evisu.yx(:,iim),handles.evisu.xy(:,iim),handles.evisu.yy(:,iim)],coords,1);
                    if handles.sonelt
                        E=griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xx,coords.xi,coords.yi,'nearest');
                        E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yx,coords.xi,coords.yi,'nearest')];
                        E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xy,coords.xi,coords.yi,'nearest')];
                        E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yy,coords.xi,coords.yi,'nearest')];
                    else
                        E=interpMesh(handles.mvisu,[handles.evisu.xx,handles.evisu.yx,handles.evisu.xy,handles.evisu.yy],coords,1);
                    end
                else
                    E=zeros(length(s),4);
                end
                fgage=figure(100+id);
                if ~isempty(findobj(gca,'DisplayName','\epsilon_{xx}'))
                    set(findobj(gca,'DisplayName','\epsilon_{xx}'),'Ydata',E(:,1));
                    set(findobj(gca,'DisplayName','\epsilon_{xx}'),'Xdata',s);
                    set(findobj(gca,'DisplayName','\epsilon_{yy}'),'Ydata',E(:,4));
                    set(findobj(gca,'DisplayName','\epsilon_{yy}'),'Xdata',s);
                    set(findobj(gca,'DisplayName','\epsilon_{xy}'),'Ydata',0.5*(E(:,2)+E(:,3)));
                    set(findobj(gca,'DisplayName','\epsilon_{xy}'),'Xdata',s);
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                    
                else
                    delete(gca)
                    axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                        'FontName','Times');
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                    set(fgage,'NumberTitle','off');
                    box('on');
                    hold('all');
                    plot(s,E(:,1),'DisplayName','\epsilon_{xx}','Parent',axes1,'LineWidth',2);
                    plot(s,E(:,4),'DisplayName','\epsilon_{yy}','Parent',axes1,'LineWidth',2);
                    plot(s,0.5*(E(:,2)+E(:,3)),'DisplayName','\epsilon_{xy}','Parent',axes1,'LineWidth',2);
                    ylabel('Strain []','FontSize',20,'FontName','Times');
                    xlabel('s [pixel]','FontSize',20,'FontName','Times');
                    leg=legend(axes1,'show','Location','NorthWest');
                    set(leg,'Box','off')
                end
                if export
                    [pp,filename,ext]=fileparts(handles.param.result_file);
                    filexp=sprintf('%s-line-%02d-strain.csv',filename,id);
                    if exist(filexp,'file')
                        set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                        [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save gage data...',filexp);
                    end
                    
                    fid=fopen(filexp,'w');
                    fprintf(fid,'Result file;%s\n',handles.param.result_file);
                    if ~strcmp(handles.param.analysis,'mechanics')
                        fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
                    end
                    fprintf(fid,'Line plot;%d\n',id);
                    fprintf(fid,'X1 [pixel];Y1 [pixel];X2 [pixel];Y2 [pixel]; Orientation; Length [pixel];\n');
                    fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:),gage(2,:),(angle(angl))*180/pi,l1);
                    if ~(handles.param.pixel_size==1)
                        pix2m=handles.param.pixel_size;
                        fprintf(fid,'X1 [m];Y1 [m];X2 [m];Y2 [m]; Orientation; Length [m];\n');
                        fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:)*pix2m,gage(2,:)*pix2m,(angle(angl))*180/pi,l1*pix2m);
                        fprintf(fid,'Position X [m];Position Y [m]; Orientation; Length [m]; Width [m]\n');
                        fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1))*pix2m,mean(gage(1:4,2))*pix2m,(angle(angl)+pi/2)*180/pi,max(l1,l2)*pix2m,min(l1,l2)*pix2m);
                    end
                    
                    if isfield(handles.param,'deformed_image')
                        fprintf(fid,'Filename;;;');
                        for iim=1:handles.animation.nbstep
                            
                            if iscell(handles.param.deformed_image)
                                fprintf(fid,'"%s";;;',handles.param.deformed_image{1,iim+handles.stereo});
                            else
                                fprintf(fid,'"%s";;;',handles.param.deformed_image);
                            end
                        end
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,'Step;;;');
                    for iim=1:handles.animation.nbstep
                        fprintf(fid,'%d;;;',iim);
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'x [pixel];y [pixel];s [pixel];Exx[];Eyy[];Exy[]\n');
                    
                    Exx=zeros(length(s),handles.animation.nbstep);
                    Eyy=zeros(length(s),handles.animation.nbstep);
                    Exy=zeros(length(s),handles.animation.nbstep);
                    for iim=1:handles.animation.nbstep
                        handles.animation.iim=iim;
                        handles=ComputeStrain(handles);
                        if handles.sonelt
                            E=griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xx,coords.xi,coords.yi,'nearest');
                            E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yx,coords.xi,coords.yi,'nearest')];
                            E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xy,coords.xi,coords.yi,'nearest')];
                            E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yy,coords.xi,coords.yi,'nearest')];
                        else
                            E=interpMesh(handles.mvisu,[handles.evisu.xx,handles.evisu.yx,handles.evisu.xy,handles.evisu.yy],coords,1);
                        end
                        Exx(:,iim)=E(:,1);
                        Eyy(:,iim)=E(:,4);
                        Exy(:,iim)=0.5*(E(:,2)+E(:,3));
                        
                    end
                    
                    for ip=1:length(s)
                        fprintf(fid,'%.3f;%.3f;%.3f;',real(xyp(ip)),imag(xyp(ip)),s(ip));
                        ijm=handles.animation.iim;
                        for iim=1:handles.animation.nbstep
                            fprintf(fid,'%.3e;%.3e;%.3e;',Exx(ip,iim),Eyy(ip,iim),Exy(ip,iim));
                        end
                        fprintf(fid,'\n');
                    end
                    
                    fclose(fid);
                    handles.animation.iim=ijm;
                    
                    handles=ComputeStrain(handles);
                    
                    param=handles.param;
                    model=handles.fem_model;
                    save(sprintf('%s-line-%02d-strain.res',filename,id),'Exx','Eyy','Exy','gage','param','model','-v7.3')
                    set(handles.message_text,'String','Exporting line data.....done')
                    
                    
                end
            case 3
                if iim
                    fide=fopen([strrep(handles.param.result_file,'.res',''),'-error.res'],'r');
                    erroronelt=fread(fide,1);
                    dynamic=fread(fide,1);
                    if erroronelt
                        fseek(fide,(iim-1+handles.stereo*(iim+1))*prod(handles.mvisu.Nelems)+2,'bof');
                        Eri=fread(fide,prod(handles.mvisu.Nelems));
                        fclose(fide);
                        %                      Eri=handles.eton*Eri;
                        Eri=100*Eri/dynamic;
                        if ~isfield(handles.evisu,'xg')
                            handles=ComputeStrain(handles);
                            guidata(handles.figure1,handles);
                        end
                        %                      Er=interpMesh(handles.mvisu,Eri,coords,1);
                        Er=griddata(handles.evisu.xg,handles.evisu.yg,Eri,coords.xi,coords.yi,'nearest');
                        
                    else
                        sizeim=[roi(2)-roi(1),roi(4)-roi(3)]+1;
                        fseek(fide,(iim-1+handles.stereo*(iim+1))*prod(sizeim)+2,'bof');
                        Eri=fread(fide,prod(sizeim));
                        fclose(fide);
                        Er=NaN(length(s),1);
                        in=inpolygon(coords.xi,coords.yi,roi([1,2,2,1,1,])-roi(1)+1,roi([3,3,4,4,3])-roi(3)+1);
                        ind=sub2ind(sizeim,round(coords.xi(in)),round(coords.yi(in)));
                        Er(in)=100*Eri(ind)/dynamic;
                    end
                    
                    
                else
                    Er=zeros(length(s),1);
                end
                fgage=figure(300+id);
                if ~isempty(findobj(gca,'DisplayName','Er'))
                    set(findobj(gca,'DisplayName','Er'),'Ydata',Er);
                    set(findobj(gca,'DisplayName','Er'),'Xdata',s);
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                else
                    delete(gca)
                    axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                        'FontName','Times');
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                    set(fgage,'NumberTitle','off');
                    box('on');
                    hold('all');
                    plot(s,Er,'DisplayName','Er','Parent',axes1,'LineWidth',2);
                    ylabel('Correlation error [% dyn]','FontSize',20,'FontName','Times');
                    xlabel('s [pixel]','FontSize',20,'FontName','Times');
                end
                
                if export
                    Er=NaN(length(s),size(handles.uvisu,2));
                    if ~erroronelt      sizeim=[roi(2)-roi(1),roi(4)-roi(3)]+1;
                        in=inpolygon(coords.xi,coords.yi,roi([1,2,2,1,1,])-roi(1)+1,roi([3,3,4,4,3])-roi(3)+1);
                        ind=sub2ind(sizeim,round(coords.xi(in)),round(coords.yi(in)));
                    end
                    fide=fopen([strrep(handles.param.result_file,'.res',''),'-error.res'],'r');
                    erroronelt=fread(fide,1);
                    dynamic=fread(fide,1);
                    for iim=1:size(Er,2)
                        if erroronelt
                            Eri=fread(fide,prod(handles.mvisu.Nelems));
                            %                            Eri=handles.eton*Eri;
                            Eri=100*Eri/dynamic;
                            Er(:,iim)=griddata(handles.evisu.xg,handles.evisu.yg,Eri,coords.xi,coords.yi,'nearest');
                            
                            % Er(:,iim)=interpMesh(handles.mvisu,Eri,coords,1);
                            
                        else
                            Eri=fread(fide,prod(sizeim));
                            Er(in,iim)=100*Eri(ind)/dynamic;
                        end
                    end
                    fclose(fide);
                    
                    
                    
                    
                    
                    [pp,filename,ext]=fileparts(handles.param.result_file);
                    filexp=sprintf('%s-line-%02d-err.csv',filename,id);
                    if exist(filexp,'file')
                        set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                        [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save gage data...',filexp);
                    end
                    
                    fid=fopen(filexp,'w');
                    fprintf(fid,'Result file;%s\n',handles.param.result_file);
                    if ~strcmp(handles.param.analysis,'mechanics')
                        fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
                    end
                    fprintf(fid,'Line plot;%d\n',id);
                    fprintf(fid,'X1 [pixel];Y1 [pixel];X2 [pixel];Y2 [pixel]; Orientation; Length [pixel];\n');
                    fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:),gage(2,:),(angle(angl))*180/pi,l1);
                    if ~(handles.param.pixel_size==1)
                        pix2m=handles.param.pixel_size;
                        fprintf(fid,'X1 [m];Y1 [m];X2 [m];Y2 [m]; Orientation; Length [m];\n');
                        fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:)*pix2m,gage(2,:)*pix2m,(angle(angl))*180/pi,l1*pix2m);
                        fprintf(fid,'Position X [m];Position Y [m]; Orientation; Length [m]; Width [m]\n');
                        fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1))*pix2m,mean(gage(1:4,2))*pix2m,(angle(angl)+pi/2)*180/pi,max(l1,l2)*pix2m,min(l1,l2)*pix2m);
                    end
                    
                    if isfield(handles.param,'deformed_image')
                        fprintf(fid,'Filename;;;');
                        for iim=1:size(Er,2)
                            if iscell(handles.param.deformed_image)
                                fprintf(fid,'"%s";',handles.param.deformed_image{1,iim+handles.stereo});
                            else
                                fprintf(fid,'"%s";',handles.param.deformed_image);
                            end
                        end
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,'Step;;;');
                    for iim=1:size(Er,2)
                        fprintf(fid,'%d;',iim);
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'x [pixel];y [pixel];s [pixel];Er[%%dyn]\n');
                    for ip=1:length(s)
                        fprintf(fid,'%.3f;%.3f;%.3f;',real(xyp(ip)),imag(xyp(ip)),s(ip));
                        for iim=1:size(Er,2)
                            fprintf(fid,'%.3e;',Er(ip,iim));
                        end
                        fprintf(fid,'\n');
                    end
                    fclose(fid);
                    
                    
                    param=handles.param;
                    model=handles.fem_model;
                    save(sprintf('%s-line-%02d-err.res',filename,id),'Er','gage','param','model','-v7.3')
                    set(handles.message_text,'String','Exporting line data.....done')
                    
                    
                end
            case 4
                if handles.stereo
                    Zz=interpMesh(handles.mvisu,handles.mvisu.Zo,coords,1);
                    
                    fgage=figure(200+id);
                    if ~isempty(findobj(gca,'DisplayName','Z'))
                        set(findobj(gca,'DisplayName','Z'),'Ydata',Uxyz(:,4));
                        set(findobj(gca,'DisplayName','Z'),'Xdata',s);
                        set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                    else
                        delete(gca)
                        axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                            'FontName','Times');
                        set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                        set(fgage,'NumberTitle','off');
                        box('on');
                        hold('all');
                        plot(s,Zz,'DisplayName','Z','Parent',axes1,'LineWidth',2);
                        ylabel('Topography [m]','FontSize',20,'FontName','Times');
                        xlabel('s [pixel]','FontSize',20,'FontName','Times');
                    end
                    if export
                        
                        
                        [pp,filename,ext]=fileparts(handles.param.result_file);
                        filexp=sprintf('%s-line-%02d-topo.csv',filename,id);
                        if exist(filexp,'file')
                            set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                            [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save gage data...',filexp);
                        end
                        
                        fid=fopen(filexp,'w');
                        fprintf(fid,'Result file;%s\n',handles.param.result_file);
                        if ~strcmp(handles.param.analysis,'mechanics')
                            fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
                        end
                        fprintf(fid,'Line plot;%d\n',id);
                        fprintf(fid,'X1 [pixel];Y1 [pixel];X2 [pixel];Y2 [pixel]; Orientation; Length [pixel];\n');
                        fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:),gage(2,:),(angle(angl))*180/pi,l1);
                        if ~(handles.param.pixel_size==1)
                            pix2m=handles.param.pixel_size;
                            fprintf(fid,'X1 [m];Y1 [m];X2 [m];Y2 [m]; Orientation; Length [m];\n');
                            fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:)*pix2m,gage(2,:)*pix2m,(angle(angl))*180/pi,l1*pix2m);
                            fprintf(fid,'Position X [m];Position Y [m]; Orientation; Length [m]; Width [m]\n');
                            fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1))*pix2m,mean(gage(1:4,2))*pix2m,(angle(angl)+pi/2)*180/pi,max(l1,l2)*pix2m,min(l1,l2)*pix2m);
                        end
                        fprintf(fid,'x [pixel];y [pixel];s [pixel];Z[m];\n');
                        for ip=1:length(s)
                            fprintf(fid,'%.3f;%.3f;%.3f;%.3e;',real(xyp(ip)),imag(xyp(ip)),s(ip),Zz(ip));
                            fprintf(fid,'\n');
                        end
                        
                        fclose(fid);
                        
                        
                        param=handles.param;
                        model=handles.fem_model;
                        save(sprintf('%s-line-%02d-topo.res',filename,id),'Zz','gage','param','model','-v7.3')
                        set(handles.message_text,'String','Exporting line data.....done')
                        
                        
                    end
                end
            case 5
                if iim
                    if handles.stereo
                        U=handles.uxyz(:,iim);
                        U=[reshape(U,prod(Nnodes),3)];
                    else
                        U=handles.uvisu(:,iim);
                        U=reshape(U(1:2*prod(Nnodes)),prod(Nnodes),2);
                    end
                    if handles.rmrbm
                        if handles.stereo
                            Xo=handles.mvisu.Xo;
                            Yo=handles.mvisu.Yo;
                            Zo=handles.mvisu.Zo;
                            L=[1+0*Xo,0*Xo,0*Xo,-Yo,Zo,0*Xo;...
                                0*Yo,1+0*Yo,0*Yo,Xo,0*Yo,-Zo;...
                                0*Zo,0*Zo,1+0*Zo,0*Zo,-Xo,Yo];
                            Uxo=L*(L\[U(:,1);U(:,2);U(:,3)]);
                            Uyo=Uxo(length(Xo)+(1:length(Xo)));
                            Uzo=Uxo(2*length(Xo)+(1:length(Xo)));
                            Uxo=Uxo((1:length(Xo)));
                            U=U-[Uxo,Uyo,Uzo];
                        else
                            L=[ones(prod(Nnodes),1),zeros(prod(Nnodes),1),-handles.mvisu.yo;...
                                zeros(prod(Nnodes),1),ones(prod(Nnodes),1),handles.mvisu.xo];
                            Uxo=L*(L\[U(:,1);U(:,2)]);
                            Uyo=Uxo(prod(Nnodes)+(1:prod(Nnodes)));
                            Uxo=Uxo((1:prod(Nnodes)));
                            U=U-[Uxo,Uyo];
                        end
                    end
                    
                    Uz=interpMesh(handles.mvisu,U,coords,1);
                    if handles.sonelt
                        Exy=griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xx,coords.xi,coords.yi,'nearest');
                        Exy=[Exy,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yx,coords.xi,coords.yi,'nearest')];
                        Exy=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xy,coords.xi,coords.yi,'nearest')];
                        Exy=[Exy,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yy,coords.xi,coords.yi,'nearest')];
                    else
                        Exy=interpMesh(handles.mvisu,[handles.evisu.xx,handles.evisu.yx,handles.evisu.xy,handles.evisu.yy],coords,1);
                    end
                else
                    Uxyz=zeros(length(s),3);
                    Exy=zeros(length(s),3);
                end
                Ux=Uz(:,1);
                Uy=Uz(:,2);
                if size(Uz,2)==3
                    Uz=Uz(:,3);
                else
                    Uz=0*Ux;
                end
                Exx=Exy(:,1);
                Eyy=Exy(:,2);
                Exy=0.5*Exy(:,3);
                form=strrep(handles.calculator{2,handles.ccomp+1},'*','.*');
                form=strrep(form,'/','./');
                form=strrep(form,'^','.^');
                nU=eval(form);
                tit=sprintf('%s',handles.calculator{1,handles.ccomp+1});
                fgage=figure(500+id);
                if ~isempty(findobj(gca,'DisplayName',tit))
                    set(findobj(gca,'DisplayName',tit),'Ydata',nU);
                    set(findobj(gca,'DisplayName',tit),'Xdata',s);
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                else
                    delete(gca)
                    axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                        'FontName','Times');
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                    set(fgage,'NumberTitle','off');
                    box('on');
                    hold('all');
                    plot(s,nU,'DisplayName',tit,'Parent',axes1,'LineWidth',2);
                    ylabel([tit,sprintf('[%s]',handles.calculator{4,handles.ccomp+1})],'FontSize',20,'FontName','Times');
                    xlabel('s [pixel]','FontSize',20,'FontName','Times');
                    leg=legend(axes1,'show','Location','NorthWest');
                    set(leg,'Box','off')
                end
                if export
                    if handles.stereo
                        U=handles.uxyz;
                    else
                        U=handles.uvisu;
                    end
                    Nnodes=handles.mvisu.Nnodes;
                    Ux=interpMesh(handles.mvisu,U((1:prod(Nnodes)),:),coords);
                    Uy=interpMesh(handles.mvisu,U(prod(Nnodes)+(1:prod(Nnodes)),:),coords);
                    if handles.stereo
                        Uz=interpMesh(handles.mvisu,U(2*prod(Nnodes)+(1:prod(Nnodes)),:),coords);
                    end
                    Exx=zeros(length(s),handles.animation.nbstep);
                    Eyy=zeros(length(s),handles.animation.nbstep);
                    Exy=zeros(length(s),handles.animation.nbstep);
                    for iim=1:handles.animation.nbstep
                        handles.animation.iim=iim;
                        handles=ComputeStrain(handles);
                        if handles.sonelt
                            E=griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xx,coords.xi,coords.yi,'nearest');
                            E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yx,coords.xi,coords.yi,'nearest')];
                            E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.xy,coords.xi,coords.yi,'nearest')];
                            E=[E,griddata(handles.evisu.xg,handles.evisu.yg,handles.evisu.yy,coords.xi,coords.yi,'nearest')];
                        else
                            E=interpMesh(handles.mvisu,[handles.evisu.xx,handles.evisu.yx,handles.evisu.xy,handles.evisu.yy],coords,1);
                        end
                        Exx(:,iim)=E(:,1);
                        Eyy(:,iim)=E(:,4);
                        Exy(:,iim)=0.5*(E(:,2)+E(:,3));
                        
                    end
                    data=eval(form);
                    
                    [pp,filename,ext]=fileparts(handles.param.result_file);
                    filexp=sprintf('%s-line-%02d-%s.csv',filename,id,tit);
                    if exist(filexp,'file')
                        set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                        [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save line data...',filexp);
                    end
                    
                    fid=fopen(filexp,'w');
                    fprintf(fid,'Result file;%s\n',handles.param.result_file);
                    if ~strcmp(handles.param.analysis,'mechanics')
                        fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
                    end
                    fprintf(fid,'Calculator;%s\n',handles.calculator{2,handles.ccomp+1});
                    fprintf(fid,'Line plot;%d\n',id);
                    fprintf(fid,'X1 [pixel];Y1 [pixel];X2 [pixel];Y2 [pixel]; Orientation; Length [pixel];\n');
                    fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:),gage(2,:),(angle(angl))*180/pi,l1);
                    if ~(handles.param.pixel_size==1)
                        pix2m=handles.param.pixel_size;
                        fprintf(fid,'X1 [m];Y1 [m];X2 [m];Y2 [m]; Orientation; Length [m];\n');
                        fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:)*pix2m,gage(2,:)*pix2m,(angle(angl))*180/pi,l1*pix2m);
                        fprintf(fid,'Position X [m];Position Y [m]; Orientation; Length [m]; Width [m]\n');
                        fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1))*pix2m,mean(gage(1:4,2))*pix2m,(angle(angl)+pi/2)*180/pi,max(l1,l2)*pix2m,min(l1,l2)*pix2m);
                    end
                    if isfield(handles.param,'deformed_image')
                        fprintf(fid,'Filename;;;');
                        for iim=1:size(Uy,2)
                            if iscell(handles.param.deformed_image)
                                fprintf(fid,'"%s";',handles.param.deformed_image{1,iim+handles.stereo});
                            else
                                fprintf(fid,'"%s";',handles.param.deformed_image);
                            end
                        end
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,'Step;;;');
                    for iim=1:size(Uy,2)
                        fprintf(fid,'%d;',iim);
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'x [pixel];y [pixel];s [pixel];data[%s];\n',handles.calculator{4,handles.ccomp+1});
                    for ip=1:length(s)
                        fprintf(fid,'%.3f;%.3f;%.3f;',real(xyp(ip)),imag(xyp(ip)),s(ip));
                        for iim=1:size(Uy,2)
                            fprintf(fid,'%.3e;',data(ip,iim));
                        end
                        fprintf(fid,'\n');
                    end
                    
                    fclose(fid);
                    
                    
                    param=handles.param;
                    model=handles.fem_model;
                    save(sprintf('%s-line-%02d-%s.res',filename,id,tit),'data','gage','param','model','-v7.3')
                    set(handles.message_text,'String','Exporting line data.....done')
                    
                    
                end
                
                
            case 6
                if iim
                    U=handles.uvisu(2*prod(Nnodes)+(1:prod(Nnodes)),iim);
                    U=interpMesh(handles.mvisu,U,coords,1);
                else
                    U=zeros(length(s),1);
                end
                
                fgage=figure(600+id);
                if ~isempty(findobj(gca,'DisplayName','T'))
                    set(findobj(gca,'DisplayName','T'),'Ydata',U);
                    set(findobj(gca,'DisplayName','T'),'Xdata',s);
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                else
                    delete(gca)
                    axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                        'FontName','Times');
                    set(fgage,'Name',sprintf('Line #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                    set(fgage,'NumberTitle','off');
                    box('on');
                    hold('all');
                    plot(s,U,'DisplayName','T','Parent',axes1,'LineWidth',2);
                    ylabel('T [DL]','FontSize',20,'FontName','Times');
                    xlabel('s [pixel]','FontSize',20,'FontName','Times');
                    leg=legend(axes1,'show','Location','NorthWest');
                    set(leg,'Box','off')
                end
                if export
                    U=handles.uvisu;
                    Nnodes=handles.mvisu.Nnodes;
                    U=interpMesh(handles.mvisu,U(2*prod(Nnodes)+(1:prod(Nnodes)),:),coords);
                    [pp,filename,ext]=fileparts(handles.param.result_file);
                    filexp=sprintf('%s-line-%02d-temp.csv',filename,id);
                    if exist(filexp,'file')
                        set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                        [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save line data...',filexp);
                    end
                    
                    fid=fopen(filexp,'w');
                    fprintf(fid,'Result file;%s\n',handles.param.result_file);
                    if ~strcmp(handles.param.analysis,'mechanics')
                        fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
                    end
                    fprintf(fid,'Line plot;%d\n',id);
                    fprintf(fid,'X1 [pixel];Y1 [pixel];X2 [pixel];Y2 [pixel]; Orientation; Length [pixel];\n');
                    fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:),gage(2,:),(angle(angl))*180/pi,l1);
                    if ~(handles.param.pixel_size==1)
                        pix2m=handles.param.pixel_size;
                        fprintf(fid,'X1 [m];Y1 [m];X2 [m];Y2 [m]; Orientation; Length [m];\n');
                        fprintf(fid,'%f;%f;%f;%f;%f;%f\n',gage(1,:)*pix2m,gage(2,:)*pix2m,(angle(angl))*180/pi,l1*pix2m);
                        fprintf(fid,'Position X [m];Position Y [m]; Orientation; Length [m]; Width [m]\n');
                        fprintf(fid,'%f;%f;%f;%f;%f\n',mean(gage(1:4,1))*pix2m,mean(gage(1:4,2))*pix2m,(angle(angl)+pi/2)*180/pi,max(l1,l2)*pix2m,min(l1,l2)*pix2m);
                    end
                    if isfield(handles.param,'deformed_image')
                        fprintf(fid,'Filename;;;');
                        for iim=1:size(U,2)
                            if iscell(handles.param.deformed_image)
                                fprintf(fid,'"%s";',handles.param.deformed_image{1,iim+handles.stereo});
                            else
                                fprintf(fid,'"%s";',handles.param.deformed_image);
                            end
                        end
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,'Step;;;');
                    for iim=1:size(U,2)
                        fprintf(fid,'%d;',iim);
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'x [pixel];y [pixel];s [pixel];T[DL];\n');
                    for ip=1:length(s)
                        fprintf(fid,'%.3f;%.3f;%.3f;',real(xyp(ip)),imag(xyp(ip)),s(ip));
                        for iim=1:size(U,2)
                            fprintf(fid,'%.3e;',U(ip,iim));
                        end
                        fprintf(fid,'\n');
                    end
                    fclose(fid);
                    
                    T=U;
                    param=handles.param;
                    model=handles.fem_model;
                    save(sprintf('%s-line-%02d-temp.res',filename,id),'T','gage','param','model','-v7.3')
                    set(handles.message_text,'String','Exporting line data.....done')
                    
                    
                end
        end
    end
    set(0,'CurrentFigure',handles.figure1);
end


function zone_set(handles,attractor)
if nargin<2,attractor=0;end
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Click opposite corners of the rectangular zone.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    
    handles.initial_point=point1;
    guidata(handles.figure1,handles);
    xyp=[point1;point1];
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_rect);
    
    waitforbuttonpress
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    
    point2=getposition(handles);
    point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
    ly=abs((point1(1)-point2(1)));
    lx=abs((point1(2)-point2(2)));
    if strcmp(get(handles.figure1,'selectiontype'),'normal')&&(~(max(lx,ly)<1))
        delete(gz);
        roi=[max(1,min(handles.sizeim(1),sort(round([point1(1),point2(1)])))),max(1,min(handles.sizeim(2),sort(round([point1(2),point2(2)]))))];
        xyp=[roi([1,2,2,1,1])',roi([3,3,4,4,3])'];
        switch handles.fbasis
            case 'fem'
                gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','-','LineWidth',2,'ButtonDownFcn',@zone_adjust);
                gz.Tag=num2str(rand(1));
                handles.fem_model.zone{2,end+1}=xyp;
                handles.fem_model.zone{3,end}=gz;
                handles.fem_model.zone{4,end}=1;
                if ~attractor
                    handles.fem_model.zone{1,end}=1;
                    handles.fem_model.zone{6,end}=0;
                else
                    handles.fem_model.zone{1,end}=-attractor;
                    handles.fem_model.zone{6,end}=mean(handles.fem_model.mesh_size);
                end
                
                switch handles.ana
                    case {1,3}
                        if ~attractor
                            set(gz,'uicontextmenu',handles.cmenu_zone);
                        else
                            set(gz,'uicontextmenu',handles.cmenu_attractor);
                            set(gz,'Color',[0,0,0.5]);
                        end
                        if handles.fem_model.mesh_type==2
                            handles=remesh(handles);
                        else
                            handles=plot_mesh(handles);
                        end
                        handles=set_param(handles);
                        
                        set(handles.message_text,'String','Creating rectangular zone.....done')
                    case 2
                        
                        handles=set_param(handles);
                        handles=ComputeStrainGage(handles,xyp,size(handles.fem_model.zone,2));
                        plot_gage_data(handles,xyp,size(handles.fem_model.zone,2))
                        set(gz,'uicontextmenu',handles.cmenu_export_gage)
                        set(handles.message_text,'String','Creating strain gage.....done')
                end
                
            case 'uni'
                gz=plot(xyp(:,1),xyp(:,2),'Color',[0.5,0,0.75],'LineStyle','-','LineWidth',2,'uicontextmenu',handles.cmenu_gage,'ButtonDownFcn',@zone_adjust);
                gz.Tag=num2str(rand(1));
                handles.uni_model.zone{2,end+1}=xyp;
                handles.uni_model.zone{3,end}=gz;
                handles.uni_model.zone{1,end}=1;
                handles.uni_model.zone{4,end}=1;
                handles=set_param(handles);
                plot_gage_data(handles,xyp,size(handles.uni_model.zone,2))
                set(handles.message_text,'String','Creating strain gage.....done')
            case 'beam'
                gz=plot(xyp(:,1),xyp(:,2),'Color',[0.75,0,0],'LineStyle','-','LineWidth',2,'uicontextmenu',handles.cmenu_beam,'ButtonDownFcn',@zone_adjust);
                gz.Tag=num2str(rand(1));
                handles.beam_model.zone{2,end+1}=xyp;
                handles.beam_model.zone{3,end}=gz;
                handles.beam_model.zone{1,end}=1;
                handles.beam_model.zone{4,end}=1;
                handles=set_param(handles);
                set(handles.message_text,'String','Creating beam.....done')
                
        end
    else
        delete(gz);
        set(handles.message_text,'String','Zone definition.....cancelled')
    end
else
    set(handles.message_text,'String','Zone definition.....cancelled')
end



function zone_adjust(hObject,eventdata)
if ~strcmp(get(gcf,'selectiontype'),'alt')
    handles=guidata(hObject);
    pt=getposition(handles);
    handles.initial_point=pt;
    guidata(hObject,handles);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonUpFcn',@stop_follow);
    set(gcf,'WindowButtonMotionFcn',@follow_mouse);
end
function follow_mouse(hObject,eventdata)
handles=guidata(hObject);
po=handles.initial_point;
pt=getposition(handles);


try
    switch handles.fbasis
        case 'fem'
            iz=findg((handles.fem_model.zone(3,:)),gco);
            xyp=handles.fem_model.zone{2,iz};
            typ=handles.fem_model.zone{4,iz};
        case 'uni'
            iz=findg((handles.uni_model.zone(3,:)),gco);
            xyp=handles.uni_model.zone{2,iz};
            typ=handles.uni_model.zone{4,iz};
        case 'beam'
            iz=findg((handles.beam_model.zone(3,:)),gco);
            xyp=handles.beam_model.zone{2,iz};
            typ=handles.beam_model.zone{4,iz};
        case 'vic'
            iz=findg((handles.vic_model.zone(3,:)),gco);
            xyp=handles.vic_model.zone{2,iz};
            typ=handles.vic_model.zone{4,iz};
    end
catch
    assert(handles.field==3);
    iz=1;
    xyp=handles.erzone{2,1};
    typ=handles.erzone{4,1};
    
end
np=size(xyp,1);
Zcp=(xyp*[1;1i]);
switch get(handles.figure1,'selectiontype')
    case 'extend'
        set(gco,'Xdata',xyp(:,1)+pt(1,1)-po(1),'Ydata',xyp(:,2)+pt(1,2)-po(2))
    case 'normal'
        po=po*[1;1i];
        pt=pt*[1;1i];
        switch typ
            case {1,6} %rectangle
                Zcp=[Zcp(1:np);0.5*(Zcp(1:np-1)+Zcp(2:np))];
                [dmin,icp]=min(abs(Zcp-po));
                icp=icp(1);
                po=Zcp(icp);
                if icp<=np
                    if ~(handles.field==3)
                        pc=mean(Zcp(1:np-1));
                        rot=((pt-pc)/abs(pt-pc)*abs(po-pc)/(po-pc));
                        Zcp=(Zcp-pc)*rot+pc;
                    end
                else
                    d=pt-po;
                    seg=diff(Zcp(icp-np+(0:1)));
                    seg=(-1i*seg)/abs(seg);
                    d=real(d'*seg)*seg;
                    Zcp(icp-np+(0:1))=Zcp(icp-np+(0:1))+d;
                    if icp-np==1
                        Zcp(np)=Zcp(1);
                    elseif icp==length(Zcp)
                        Zcp(1)=Zcp(np);
                    end
                end
                
                Zcp=Zcp(1:np);
            case 2 %poylgone
                Zcp=Zcp(1:end-1);
                [dmin,icp]=min(abs(Zcp-po));
                icp=icp(1);
                po=Zcp(icp);
                d=pt-po;
                Zcp(icp)=Zcp(icp)+d;
                Zcp(end+1)=Zcp(1);
            case 4 %line
                [dmin,icp]=min(abs(Zcp-po));
                icp=icp(1);
                po=Zcp(icp);
                d=pt-po;
                Zcp(icp)=Zcp(icp)+d;
            case 3 %circle
                [dmin,icmin]=min(abs(Zcp-po));
                [dmax,icmax]=max(abs(Zcp-po));
                p1=Zcp(icmin(1));
                p2=Zcp(icmax(1));
                seg=p2-p1;
                seg=seg/abs(seg);
                d=real(seg'*(pt-po));
                p1=p1+d*seg;
                xyp=0.5*[real(p1+p2),imag(p1+p2),norm(p2-p1)];
                
                nface=round(2*pi*xyp(2)/5);
                Zcp=xyp(1:2)*[1;1i]+xyp(3)*exp(1i*[(0:2*pi/nface:2*pi)';2*pi]);
            case 5 %crack
                [dmin,icp]=min(abs(Zcp-po));
                icp=icp(1);
                po=Zcp(icp);
                d=pt-po;
                Zcp(icp)=Zcp(icp)+d;
                
        end
        set(gco,'Xdata',real(Zcp),'Ydata',imag(Zcp))
end
function iz=findg(cellh,gz)
for iz=1:numel(cellh)
    if strcmp(get(cellh{iz},'Tag'),get(gz,'Tag')), return; end
end
error('Failed at finding the graphical object')

function stop_follow(hObject,eventdata)
handles=guidata(hObject);
set(gcf,'WindowButtonUpFcn','');
set(handles.figure1,'WindowButtonMotionFcn',@show_position);
xp=get(gco,'Xdata');
yp=get(gco,'Ydata');

try
    switch handles.fbasis
        case 'fem'
            doremesh=1;
            iz=findg((handles.fem_model.zone(3,:)),gco);
            handles.fem_model.zone{2,iz}=[xp(:),yp(:)];
            if (handles.fem_model.zone{4,iz}==3)
                p1=[xp(1),yp(1)];
                p2=[xp(round(numel(xp)/2)),yp(round(numel(xp)/2))];
                xyp=0.5*[(p1+p2),norm(p2-p1)];
                handles.fem_model.zone{5,iz}=xyp;
            elseif handles.fem_model.zone{4,iz}==6
                handles=UpdateBCSZone(handles);
                doremesh=0;
            end
            switch handles.ana
                case {1,3}
                    if handles.fem_model.mesh_type==2&&doremesh
                        handles=remesh(handles);
                    else
                        handles=plot_mesh(handles);
                    end
                    handles=set_param(handles);
                case 2
                    iz=findg((handles.fem_model.zone(3,:)),gco);
                    handles.fem_model.zone{2,iz}=[xp(:),yp(:)];
                    switch handles.fem_model.zone{4,iz}
                        case {1,4}
                            handles=ComputeStrainGage(handles,handles.fem_model.zone{2,iz},iz);
                            guidata(hObject,handles);
                            plot_gage_data(handles,[xp(:),yp(:)],iz);
                        case 5
                            guidata(hObject,handles);
                            if handles.cracked
                                plot_crack_data(handles,iz);
                            end
                    end
            end
        case 'uni'
            iz=findg((handles.uni_model.zone(3,:)),gco);
            handles.uni_model.zone{2,iz}=[xp(:),yp(:)];
            guidata(hObject,handles);
            plot_gage_data(handles,[xp(:),yp(:)],iz);
        case 'beam'
            iz=findg((handles.beam_model.zone(3,:)),gco);
            handles.beam_model.zone{2,iz}=[xp(:),yp(:)];
            guidata(hObject,handles);
        case 'vic'
            iz=findg((handles.vic_model.zone(3,:)),gco);
            handles.vic_model.zone{2,iz}=[xp(:),yp(:)];
            guidata(hObject,handles);
    end
catch
    handles.erzone{2,1}=[xp(:),yp(:)];
    handles=plot_mesh(handles);
    handles=set_param(handles);
    
end
set(gcf,'WindowButtonDownFcn',@double_clic);

function line_set(handles,attractor)
if nargin<2,attractor=0;end
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Click both ends of the line.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    xyp=[point1;point1];
    point2=-1000;
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_poly);
    
    waitforbuttonpress
    if strcmp(get(handles.figure1,'selectiontype'),'normal')
        point2=getposition(handles);
        point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
        xyp=[xyp;point2;point2];
        set(gz,'Xdata',xyp(:,1));
        set(gz,'Ydata',xyp(:,2));
    else
        delete(gz);
        set(gcf,'WindowButtonDownFcn',@double_clic);
        set(gcf,'WindowButtonMotionFcn',@show_position);
        return;
        set(handles.message_text,'String','Zone definition.....cancelled')
    end
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    xyp=xyp(1:2:end,:);
    gz.Tag=num2str(rand(1));
    handles.fem_model.zone{2,end+1}=xyp;
    handles.fem_model.zone{3,end}=gz;
    handles.fem_model.zone{4,end}=4;
    if ~attractor
        handles.fem_model.zone{1,end}=1;
        handles.fem_model.zone{6,end}=0;
    else
        handles.fem_model.zone{1,end}=-1;
        handles.fem_model.zone{6,end}=mean(handles.fem_model.mesh_size);
    end
    
    set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
    set(gz,'Color',[0,0,0.5]);
    switch handles.ana
        case {1,3}
            set(gz,'uicontextmenu',handles.cmenu_attractor);
            
            if handles.fem_model.mesh_type==2
                handles=remesh(handles);
            else
                handles=plot_mesh(handles);
            end
            handles=set_param(handles);
        case 2
            handles=set_param(handles);
            plot_gage_data(handles,xyp,size(handles.fem_model.zone,2))
            set(gz,'uicontextmenu',handles.cmenu_export_gage)
            handles=guidata(handles.figure1);
            handles=set_param(handles);
            
            
    end
    set(handles.message_text,'String','Creating line.....done')
    
else
    set(handles.message_text,'String','Line definition.....cancelled')
end
set(gcf,'WindowButtonDownFcn',@double_clic);
set(gcf,'WindowButtonMotionFcn',@show_position);

% --------------------------------------------------------------------
function cmenu_mesh_zone_poly_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_zone_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
poly_set(handles);

function poly_set(handles,attractor)
if nargin<2,attractor=0;end
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Click and close a polygone to define a zone.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    xyp=[point1;point1];
    point2=-10000000000;
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_poly);
    
    while abs((point2-xyp(1,:))*[1;1i])>0.01*min(handles.sizeim(1:2))
        waitforbuttonpress
        if strcmp(get(handles.figure1,'selectiontype'),'normal')
            point2=getposition(handles);
            point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
            xyp=[xyp;point2;point2];
            set(gz,'Xdata',xyp(:,1));
            set(gz,'Ydata',xyp(:,2));
        else
            delete(gz);
            set(gcf,'WindowButtonDownFcn',@double_clic);
            set(gcf,'WindowButtonMotionFcn',@show_position);
            return;
            set(handles.message_text,'String','Zone definition.....cancelled')
        end
        
    end
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    xyp=xyp(1:2:end,:);
    %xyp=unique(xyp,'rows','stable');
    xyp(end,:)=xyp(1,:);
    gz.Tag=num2str(rand(1));
    handles.fem_model.zone{2,end+1}=xyp;
    handles.fem_model.zone{3,end}=gz;
    handles.fem_model.zone{4,end}=2;
    if ~attractor
        handles.fem_model.zone{1,end}=1;
        handles.fem_model.zone{6,end}=0;
    else
        handles.fem_model.zone{1,end}=-attractor;
        handles.fem_model.zone{6,end}=mean(handles.fem_model.mesh_size);
    end
    set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
    if ~attractor
        set(gz,'uicontextmenu',handles.cmenu_zone);
    else
        set(gz,'Color',[0,0,0.5]);
        set(gz,'uicontextmenu',handles.cmenu_attractor);
    end
    if handles.fem_model.mesh_type==2
        handles=remesh(handles);
    else
        handles=plot_mesh(handles);
    end
    handles=set_param(handles);
    set(handles.message_text,'String','Creating polygonal zone.....done')
    
else
    set(handles.message_text,'String','Zone definition.....cancelled')
end
set(gcf,'WindowButtonDownFcn',@double_clic);
set(gcf,'WindowButtonMotionFcn',@show_position);

function follow_mouse_poly(hObject,eventdata)
handles=guidata(hObject);
pt=getposition(handles);

gz=findobj(gca,'Color',[0,0.5,0],'LineStyle','--');
xp=get(gz,'Xdata');
yp=get(gz,'Ydata');
np=length(xp);
xp=xp(1:max(1,np-1));
yp=yp(1:max(1,np-1));
set(gz,'Xdata',[xp(:);pt(1)],'Ydata',[yp(:);pt(2)])


% --------------------------------------------------------------------
function cmenu_mesh_zone_circle_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_zone_circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
circle_set(handles)
function handles=reference_frame(handles)
[az,el] = view;
dview=(az==0)&&(el==90);
if ~dview
    handles.view=[az,el];
    view([0,90]);
end
switch handles.ana
    case 1
        if handles.preview
            if handles.ondefimage
                handles.ondefimage=0;
                handles=display_frame(handles);
            end
        else
            if handles.animation.iim
                handles.animation.iim=0;
                handles=display_frame(handles);
                
            end
        end
    case 2
        if handles.ondefimage
            handles.ondefimage=0;
            handles=display_frame(handles);
        end
        
end

function circle_set(handles,attractor)
if nargin<2,attractor=0;end
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Click and drag to define a circular zone.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    
    handles.initial_point=point1;
    guidata(handles.figure1,handles);
    xyp=[point1;point1];
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_circle);
    
    waitforbuttonpress
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    if strcmp(get(handles.figure1,'selectiontype'),'normal')
        point2=getposition(handles);
        point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
        xyp=0.5*[(point1+point2),norm(point2-point1)];
        nface=round(2*pi*xyp(2)/5);
        Zcp=xyp(1:2)*[1;1i]+xyp(3)*exp(1i*[(0:2*pi/nface:2*pi)';2*pi]);
        Zcp=[real(Zcp),imag(Zcp)];
        set(gz,'Xdata',Zcp(:,1));
        set(gz,'Ydata',Zcp(:,2));
    else
        delete(gz);
        return;
        set(handles.message_text,'String','Zone definition.....cancelled')
    end
    gz.Tag=num2str(rand(1));
    handles.fem_model.zone{2,end+1}=Zcp;
    handles.fem_model.zone{3,end}=gz;
    handles.fem_model.zone{4,end}=3;
    handles.fem_model.zone{5,end}=xyp;
    if ~attractor
        handles.fem_model.zone{1,end}=1;
        handles.fem_model.zone{6,end}=0;
    else
        handles.fem_model.zone{1,end}=-attractor;
        handles.fem_model.zone{6,end}=mean(handles.fem_model.mesh_size);
    end
    set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
    if ~attractor
        set(gz,'uicontextmenu',handles.cmenu_zone);
    else
        set(gz,'Color',[0,0,0.5]);
        set(gz,'uicontextmenu',handles.cmenu_attractor);
    end
    
    if handles.fem_model.mesh_type==2
        handles=remesh(handles);
    else
        handles=plot_mesh(handles);
    end
    handles=set_param(handles);
    set(handles.message_text,'String','Creating circular zone.....done')
    
else
    set(handles.message_text,'String','Zone definition.....cancelled')
end
set(gcf,'WindowButtonDownFcn',@double_clic);
set(gcf,'WindowButtonMotionFcn',@show_position);

function follow_mouse_rect(hObject,eventdata)
handles=guidata(hObject);
point1=handles.initial_point;
point2=getposition(handles);

gz=findobj(gca,'Color',[0,0.5,0],'LineStyle','--');
roi=[max(1,min(handles.sizeim(1),sort(round([point1(1),point2(1)])))),max(1,min(handles.sizeim(2),sort(round([point1(2),point2(2)]))))];
xyp=[roi([1,2,2,1,1])',roi([3,3,4,4,3])'];
set(gz,'Xdata',xyp(:,1));
set(gz,'Ydata',xyp(:,2));

function follow_mouse_circle(hObject,eventdata)
handles=guidata(hObject);
po=handles.initial_point;
pt=getposition(handles);

gz=findobj(gca,'Color',[0,0.5,0],'LineStyle','--');
xyp=0.5*[(po+pt),norm(pt-po)];
nface=round(2*pi*xyp(2)/5);

Zcp=xyp(1:2)*[1;1i]+xyp(3)*exp(1i*(0:2*pi/nface:2*pi));
set(gz,'Xdata',real(Zcp));
set(gz,'Ydata',imag(Zcp));


% --------------------------------------------------------------------
function cmenu_zone_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --------------------------------------------------------------------
function cmenu_animation_goto_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_goto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_mesh_elt_type_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_elt_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function cmenu_mesh_elt_type_adapt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_elt_type_adapt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fem_model.element_type=4;
handles.fem_model.mesh_type=2;
handles=remesh(handles);
handles=set_param(handles);



% --------------------------------------------------------------------
function cmenu_mesh_elt_type_quad_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_elt_type_quad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fem_model.element_type=4;
handles.fem_model.mesh_type=1;
handles=remesh(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_mesh_elt_type_tri_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_elt_type_tri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fem_model.element_type=3;
handles.fem_model.mesh_type=1;
handles=remesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_mesh_smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function cmenu_mesh_smoothing_lc_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_smoothing_lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='lc';
guidata(hObject,handles);
if strcmp(handles.param.regularization_type,'median')
    set(handles.message_text,'String','Scroll mouse wheel to adjust number of neighboors.....')
else
    set(handles.message_text,'String','Scroll mouse wheel to adjust cut-off wave length.....')
end
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);

% --------------------------------------------------------------------
function cmenu_mesh_load_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_load_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function cmenu_mesh_load_mesh_load_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_load_mesh_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fem_model.mesh_type=3;
handles=remesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_mesh_load_mesh_unload_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_load_mesh_unload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mesh_type=1;
handles=remesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
function cmenu_mesh_show_edge_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_show_edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
handles.showedge=~handles.showedge;
handles=plot_mesh(handles);
handles=set_param(handles);

function cmenu_mesh_mesh_size_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_mesh_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function scroll_adjust(hObject,eventdata)
handles=guidata(hObject);
switch handles.scroll_adjusted
    case 'step'
        iim=handles.animation.iim;
        iim=min(handles.animation.nbstep,max(0,iim-1*eventdata.VerticalScrollCount));
        handles.animation.iim=iim;
        handles=display_frame(handles);
    case 'lc'
        h=eval(get(handles.cmenu_mesh_smoothing_lc,'Label'));
        hmin=0;
        if isfield(handles.fem_model,'degree')
            hmin=round(handles.fem_model.mesh_size(1)/4);
        end
        h=max(hmin,h-1*eventdata.VerticalScrollCount);
        set(handles.cmenu_mesh_smoothing_lc,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('Cut-off: %d',h))
    case 'degree'
        h=eval(get(handles.cmenu_mesh_degree_text,'Label'));
        h=max(1,h-eventdata.VerticalScrollCount);
        set(handles.cmenu_mesh_degree_text,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('NURBS degree: %g',h))
    case 'mesh_size'
        h=eval(get(handles.cmenu_mesh_mesh_size_text,'Label'));
        h=max(0,h-(0.1+0.9*(~(h<1)))*eventdata.VerticalScrollCount);
        set(handles.cmenu_mesh_mesh_size_text,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('Mesh size: %g',h))
    case 'tip_mesh_size'
        h=handles.scroll_value;
        h=max(0,h-(0.1+0.9*(~(h<1)))*eventdata.VerticalScrollCount);
        handles.scroll_value=h;
        guidata(hObject,handles);
        set(handles.message_text,'String',sprintf('Crack tip mesh size: %g',h))
    case 'attractor_mesh_size'
        h=real(handles.scroll_value);
        h=max(0,h-(0.1+0.9*(~(h<1)))*eventdata.VerticalScrollCount);
        handles.scroll_value=h+imag(handles.scroll_value);
        guidata(hObject,handles);
        set(handles.message_text,'String',sprintf('Attractor mesh size: %g',real(h)))
    case 'face_mesh_size'
        h=handles.scroll_value;
        h=max(0,h-(0.1+0.9*(~(h<1)))*eventdata.VerticalScrollCount);
        handles.scroll_value=h;
        guidata(hObject,handles);
        set(handles.message_text,'String',sprintf('Crack face mesh size: %g',h))
    case 'sif_zone_size'
        h=handles.scroll_value;
        h=max(0,h-1*eventdata.VerticalScrollCount);
        handles.scroll_value=h;
        guidata(hObject,handles);
        set(handles.message_text,'String',sprintf('Extraction zone size: %d',h))
    case 'entropy'
        smin=handles.fem_model.smin;
        ds=0.01;
        smin=min(1,max(0,smin-ds(1)*eventdata.VerticalScrollCount));
        handles.mvisu.smin=smin;
        handles.fem_model.smin=smin;
        if handles.fem_model.mesh_type==2
            handles=remesh(handles);
        end
        handles=plot_mesh(handles);
        handles=set_param(handles);
        set(handles.message_text,'String',sprintf('Threshold: %0.2f',smin))
    case 'nbelt'
        h=eval(get(handles.cmenu_beam_nbelt_edit,'Label'));
        h=max(1,h-1*eventdata.VerticalScrollCount);
        set(handles.cmenu_beam_nbelt_edit,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('Nbr of Elt: %d',h))
    case 'vic_nbelt'
        h=eval(get(handles.cmenu_vic_nbelt_edit,'Label'));
        h=max(1,h-1*eventdata.VerticalScrollCount);
        set(handles.cmenu_vic_nbelt_edit,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('Nbr of Elt: %d',h))
    case 'vic_deg'
        h=eval(get(handles.cmenu_vic_deg_edit,'Label'));
        h=max(1,h-eventdata.VerticalScrollCount);
        set(handles.cmenu_vic_deg_edit,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('B-Spline degree: %d',h))
    case 'vic_ep'
        h=eval(get(handles.cmenu_vic_ep_edit,'Label'));
        h=max(2,h-1*eventdata.VerticalScrollCount);
        set(handles.cmenu_vic_ep_edit,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('Transition length: %d pixels',h))
    case 'vic_t'
        h=eval(get(handles.cmenu_vic_t_edit,'Label'));
        h=max(0,h-1*eventdata.VerticalScrollCount);
        set(handles.cmenu_vic_t_edit,'Label',num2str(h));
        set(handles.message_text,'String',sprintf('Thickness: %d pixels',h))
    case 'warp'
        fac=handles.wfac;
        fac=max(0,fac-(0.1+0.9*(~(fac<1)))*eventdata.VerticalScrollCount);
        handles.wfac=fac;
        handles=plot_mesh(handles);
        handles=set_param(handles);
        set(handles.message_text,'String',sprintf('Warp factor: %0.1f',fac))
    case 'elevation'
        fac=handles.zfac;
        fac=max(0,fac-(0.1+0.9*(~(fac<1)))*eventdata.VerticalScrollCount);
        handles.zfac=fac;
        handles=plot_mesh(handles);
        handles=set_param(handles);
end
function stop_adjust(hObject, eventdata)
handles=guidata(hObject);
set(handles.figure1,'WindowScrollWheelFcn','');
if strcmp(handles.scroll_adjusted,'step')
    set(handles.message_text,'String','')
else
    if strcmp(get(handles.figure1,'selectiontype'),'normal')
        
        set(handles.message_text,'String','Adjustment.....done')
        switch handles.scroll_adjusted
            case 'step'
                handles=display_frame(handles);
            case 'lc'
                h=eval(get(handles.cmenu_mesh_smoothing_lc,'Label'));
                handles.param.regularization_parameter=h;
            case 'degree'
                h=eval(get(handles.cmenu_mesh_degree_text,'Label'));
                handles.fem_model.degree=h*ones(2,1);
            case 'mesh_size'
                h=eval(get(handles.cmenu_mesh_mesh_size_text,'Label'));
                handles.fem_model.mesh_size=h*ones(2,1);
                for iz=1:size(handles.fem_model.zone,2)
                    if handles.fem_model.zone{4,iz}==5
                        handles.fem_model.zone{7,iz}=min(h,handles.fem_model.zone{7,iz});
                        handles.fem_model.zone{6,iz}=min(h,real(handles.fem_model.zone{6,iz}))+1i;
                    end
                end
                handles=remesh(handles);
            case 'attractor_mesh_size'
                h=handles.scroll_value;
                handles.fem_model.zone{6,handles.active_zone}=h+1i*imag(handles.fem_model.zone{6,handles.active_zone});
                if ~(handles.ana==2),handles=remesh(handles);end
            case 'face_mesh_size'
                h=handles.scroll_value;
                handles.fem_model.zone{6,handles.active_zone}=h+1i*imag(handles.fem_model.zone{6,handles.active_zone});
                htip=min(real(handles.fem_model.zone{7,handles.active_zone}),h);
                handles.fem_model.zone{7,handles.active_zone}=htip;
                if ~(handles.ana==2),handles=remesh(handles);end
            case 'tip_mesh_size'
                h=handles.scroll_value;
                handles.fem_model.zone{7,handles.active_zone}=h;
                hface=max(real(handles.fem_model.zone{6,handles.active_zone}),h)+1i;
                handles.fem_model.zone{6,handles.active_zone}=hface;
                if ~(handles.ana==2),handles=remesh(handles);end
            case 'sif_zone_size'
                h=handles.scroll_value;
                rc= handles.fem_model.zone{8,handles.active_zone};
                rc(1)=h;
                if h>0
                    if ~isfield(handles.fem_model,'material_parameters')
                        handles.fem_model.material_model='elastic_homogeneous_isotropic';
                        matmod.young=200e9;
                        matmod.nu=0.3;
                        handles.fem_model.material_parameters=matmod;
                        cmenu_mesh_material_Callback(handles.figure1,[], handles);
                        handles=guidata(handles.figure1);
                    end
                end
                handles.fem_model.zone{8,handles.active_zone}=rc;
                if ~(handles.ana==2),handles=remesh(handles);end
            case 'nbelt'
                h=eval(get(handles.cmenu_beam_nbelt_edit,'Label'));
                handles.beam_model.nb_element=h;
            case 'vic_nbelt'
                h=eval(get(handles.cmenu_vic_nbelt_edit,'Label'));
                handles.vic_model.zone{6,handles.active_zone}=h;
            case 'vic_deg'
                h=eval(get(handles.cmenu_vic_deg_edit,'Label'));
                handles.vic_model.zone{7,handles.active_zone}=h;
            case 'vic_ep'
                h=eval(get(handles.cmenu_vic_ep_edit,'Label'));
                handles.vic_model.zone{8,handles.active_zone}=h;
            case 'vic_t'
                h=eval(get(handles.cmenu_vic_t_edit,'Label'));
                handles.vic_model.zone{9,handles.active_zone}=h;
        end
        handles=set_param(handles);
    else
        set(handles.message_text,'String','Adjustment.....cancelled')
        switch handles.scroll_adjusted
            case 'lc'
                h=handles.param.regularization_parameter;
                set(handles.cmenu_mesh_smoothing_lc,'Label',num2str(h));
            case 'mesh_size'
                h=handles.fem_model.mesh_size(1);
                set(handles.cmenu_mesh_mesh_size_text,'Label',num2str(h));
            case 'nbelt'
                h=handles.beam_model.nb_elt;
                set(handles.cmenu_beam_nbelt_edit,'Label',num2str(h));
            case 'warp'
                handles.wfac=1;
                handles=plot_mesh(handles);
                handles=set_param(handles);
            case 'elevation'
                handles.zfac=1;
                handles=plot_mesh(handles);
                handles=set_param(handles);
        end
    end
end
set(handles.figure1,'WindowButtonDownFcn',@double_clic);

% --------------------------------------------------------------------
function cmenu_mesh_degree_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function cmenu_mesh_degree_text_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_degree_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='degree';
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust NURBS degree.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_mesh_mesh_size_text_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_mesh_size_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='mesh_size';
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust mesh size.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);

function cmenu_mesh_entropy_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_entropy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='entropy';
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust entropy erosion threshold.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);

% --------------------------------------------------------------------
function cmenu_animation_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function cmenu_animation_goto_nb_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_goto_nb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if 0
    handles.scroll_adjusted='step';
    guidata(hObject,handles);
    set(handles.message_text,'String','Scroll mouse wheel to select frame.....')
    set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
    set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);
else
    iim=handles.animation.iim;
    answer=inputdlg({'Frame'},sprintf('Go to frame...'),1,...
        {num2str(iim)});
    if ~isempty(answer)
        handles.animation.iim=min(max(eval(answer{1}),0),handles.animation.nbstep);
        handles=display_frame(handles);
    end
end

% --------------------------------------------------------------------
function cmenu_animation_goto_begin_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_goto_begin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.animation.iim=0;
handles=display_frame(handles);


% --------------------------------------------------------------------
function cmenu_animation_goto_previous_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_goto_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iim=handles.animation.iim;
handles.animation.iim=max(iim-1,0);
handles=display_frame(handles);

% --------------------------------------------------------------------
function cmenu_animation_goto_next_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_goto_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iim=handles.animation.iim;
handles.animation.iim=min(iim+1,handles.animation.nbstep);
handles=display_frame(handles);

% --------------------------------------------------------------------
function cmenu_animation_goto_end_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_goto_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.animation.iim=handles.animation.nbstep;
handles=display_frame(handles);

% --------------------------------------------------------------------
function cmenu_animation_play_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if handles.animation.playing
%     set(handles.cmenu_animation_play,'Label','Pause');
%     pause_button_Callback(handles.pause_button,[],handles);
% else
if handles.animation.playing
    pause_animation(handles);
else
    play_animation(handles);
end

function pause_animation(handles)
set(handles.cmenu_animation_play,'Label','Play');
handles.animation.playing=0;
guidata(handles.figure1,handles)
%handles=display_frame(handles);


function play_animation(handles)
set(handles.cmenu_animation_play,'Label','Pause');
iim=handles.animation.iim;
nbs=handles.animation.nbstep;
if iim==nbs, iim=0;end
handles.animation.playing=1;
while (iim<=nbs)&&handles.animation.playing
    handles.animation.iim=iim;
    handles=display_frame(handles);
    pause(min(10/nbs,1))
    handles=guidata(handles.figure1);
    iim=iim+1;
end
handles.animation.playing=0;
set(handles.cmenu_animation_play,'Label','Play');
guidata(handles.figure1,handles)
if handles.ana==2&&strcmp(handles.fbasis,'fem')
    for iz=1:size(handles.fem_model.zone,2)
        zone=handles.fem_model.zone(:,iz);
        if zone{4}==4
            plot_gage_data(handles,zone{2},iz);
        end
    end
end
% --------------------------------------------------------------------
function cmenu_param_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
function cmenu_param_model_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cmenu_param_model_nurbs,'Checked','off')
set(handles.cmenu_param_model_fem,'Checked','off')
set(handles.cmenu_param_model_uni,'Checked','off')
set(handles.cmenu_param_model_beam,'Checked','off')

switch handles.fbasis
    case 'fem'
        if isfield(handles.fem_model,'degree')
            set(handles.cmenu_param_model_nurbs,'Checked','on')
        else
            set(handles.cmenu_param_model_fem,'Checked','on')
        end
    case 'uni'
        set(handles.cmenu_param_model_uni,'Checked','on')
    case 'beam'
        set(handles.cmenu_param_model_beam,'Checked','on')
end

% --------------------------------------------------------------------
function cmenu_param_model_nurbs_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_model_nurbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fbasis='fem';
handles.fem_model.degree=1+0*handles.fem_model.mesh_size;
handles=basis_functions(handles);


% --------------------------------------------------------------------
function cmenu_param_model_fem_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_model_fem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fbasis='fem';
if isfield(handles.fem_model,'degree')
    handles.fem_model=rmfield(handles.fem_model,'degree');
end
handles=basis_functions(handles);

% --------------------------------------------------------------------
function cmenu_param_model_uni_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_model_uni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fbasis='uni';
handles=basis_functions(handles);
if isempty(handles.uni_model.zone)
    zone_set(handles)
end

% --------------------------------------------------------------------
function cmenu_param_model_beam_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_model_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fbasis='beam';
handles=basis_functions(handles);
if isempty(handles.beam_model.zone)
    zone_set(handles)
end



% --- Executes on button press in field_text.
function field_text_Callback(hObject, eventdata, handles)
% hObject    handle to field_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over field_text.
function field_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to field_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_gage_plot_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_gage_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.uni_model.zone(3,:)),gco);
gage=handles.uni_model.zone{2,iz};
plot_gage_data(handles,gage,iz)


% --------------------------------------------------------------------
function cmenu_gage_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_gage_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.uni_model.zone(3,:)),gco);
handles.uni_model.zone(:,iz)=[];
try delete(handles.fgage(iz));catch , end
try delete(100+handles.fgage(iz));catch , end
try delete(200+handles.fgage(iz));catch , end
try delete(300+handles.fgage(iz));catch , end
try delete(400+handles.fgage(iz));catch , end
try delete(4000+handles.fgage(iz));catch , end
try delete(500+handles.fgage(iz));catch , end
try delete(600+handles.fgage(iz));catch , end
guidata(hObject,handles)
delete(gco)

% --------------------------------------------------------------------
function cmenu_gage_add_gage_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_gage_add_gage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zone_set(handles);

% --------------------------------------------------------------------
function cmenu_gage_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_gage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_export_gage_plot_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_gage_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
gage=handles.fem_model.zone{2,iz};
plot_gage_data(handles,gage,iz)


% --------------------------------------------------------------------
function cmenu_field_showim_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_showim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showim=~handles.showim;
if handles.showim
    handles.wfac=1;
end
handles=display_frame(handles);


% --------------------------------------------------------------------
function cmenu_field_rmrbm_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_rmrbm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.rmrbm=~handles.rmrbm;
handles=plot_mesh(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_field_finite_strain_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_finite_strain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fstrain=~handles.fstrain;
handles=plot_mesh(handles);
handles=set_param(handles);



% --------------------------------------------------------------------
function cmenu_animation_savescreenshot_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_savescreenshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,pp,filterindex]=uiputfile({'*.png','Image file (*.png)'},'Save screenshot...');
if filterindex>0
    file=fullfile(pp,file);
    im = getframe(handles.axes1);
    imwrite(im.cdata,file);
end


% --------------------------------------------------------------------
function cmenu_animation_save_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,pp,filterindex]=uiputfile({'*.png','Image file (*.png)';'*.avi','Video file (*.avi)'},'Save animation...');

switch filterindex
    case 1
        [~,file,ext]=fileparts(file);
        nbs=handles.animation.nbstep;
        iim=0;
        handles.animation.playing=1;
        while (iim<=nbs)
            handles.animation.iim=iim;
            handles=display_frame(handles);
            pause(0.1);
            im = getframe(handles.axes1);
            filei=fullfile(pp,sprintf('%s-%04d%s',file,iim,ext));
            imwrite(im.cdata,filei);
            handles=guidata(handles.figure1);
            iim=iim+1;
        end
        handles.animation.playing=0;
        handles=guidata(handles.figure1);
        
        
    case 2
        
        nbs=handles.animation.nbstep;
        writerObj = VideoWriter(fullfile(pp,file));
        writerObj.FrameRate=max(1,round(nbs/10));
        open(writerObj);
        iim=0;
        handles.animation.playing=1;
        while (iim<=nbs)
            handles.animation.iim=iim;
            handles=display_frame(handles);
            im = getframe(handles.axes1);
            writeVideo(writerObj,im);
            handles=guidata(handles.figure1);
            iim=iim+1;
            
        end
        
        close(writerObj);
end


% --------------------------------------------------------------------
function cmenu_param_param_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if~strcmp(handles.fbasis,'none')
    handles.fsolver=USolver(handles.figure1);
    usolver=guidata(handles.fsolver);
    set(usolver.pix2m_edit,'String',num2str(handles.param.pixel_size))
    set(usolver.maxiter_edit,'String',num2str(handles.param.iter_max))
    set(usolver.convergence_edit,'String',num2str(handles.param.convergance_limit))
    set(usolver.nscale_edit,'String',num2str(handles.nscale))
    set(usolver.psample_edit,'String',num2str(handles.param.psample))
    set(usolver.do_pgd_edit,'Value',1+handles.param.do_pgd_prediction)
    set(usolver.restart_edit,'Value',~handles.param.restart)
    set(usolver.normalize_grey_levels_edit,'Value',handles.param.normalize_grey_level)
    handles=set_param(handles);
end

% --------------------------------------------------------------------
function cmenu_field_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showcb=~handles.showcb;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_beam_type_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_beam_nbelt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_nbelt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_beam_deg_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_deg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_beam_exx_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_exx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.exx=~handles.beam_model.exx;
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_beam_deg_2_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_deg_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.degree=2;
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_beam_deg_3_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_deg_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.degree=3;
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_beam_deg_4_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_deg_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.degree=4;
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_beam_deg_5_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_deg_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.degree=5;
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_beam_type_beam_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_type_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.basis='beam';
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_beam_type_beam_nurbs_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_type_beam_nurbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.basis='beam-nurbs';
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_beam_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_beam_nbelt_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_nbelt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='nbelt';
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust the number of elements.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);

% --------------------------------------------------------------------
function cmenu_export_beam_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function cmenu_export_beam_plot_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_beam_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.preview
    handles.showplot=1;
    handles.fbeam=1:4;
    guidata(handles.figure1,handles);
    
    pix2m=handles.param.pixel_size;
    load(handles.param.result_file,'-mat','s','fleche','curv','rot')
    
    dt=ceil(size(fleche,2)/7);
    id=1;
    fbeam(id)=figure(id);
    delete(gca)
    axes1 = axes('Parent',fbeam(id),'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    set(fbeam(id),'Name',sprintf('Beam deflection'));
    set(fbeam(id),'NumberTitle','off');
    box('on');
    hold('all');
    
    for iim=1:dt:size(fleche,2)
        plot(s'*pix2m,fleche(:,iim)*pix2m,'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
    end
    if pix2m==1
        ylabel('Deflection [pixel]','FontSize',20,'FontName','Times');
        xlabel('s [pixel]','FontSize',20,'FontName','Times');
    else
        ylabel('Deflection [m]','FontSize',20,'FontName','Times');
        xlabel('s [m]','FontSize',20,'FontName','Times');
    end
    %    if size(fleche,2)<6
    leg=legend(axes1,'show','Location','NorthEast');
    set(leg,'Box','off')
    %    end
    
    
    if strcmp(handles.beam_model.basis,'beam')
        id=2;
        fbeam(id)=figure(id);
        delete(gca)
        axes1 = axes('Parent',fbeam(id),'LineWidth',2,'FontSize',16,...
            'FontName','Times');
        set(fbeam(id),'Name',sprintf('Beam curvature'));
        set(fbeam(id),'NumberTitle','off');
        box('on');
        hold('all');
        
        plot(curv/pix2m,'Parent',axes1,'LineWidth',2);
        if pix2m==1
            ylabel('Curvature [1/pixel]','FontSize',20,'FontName','Times');
            xlabel('Step','FontSize',20,'FontName','Times');
        else
            ylabel('Curvature [1/m]','FontSize',20,'FontName','Times');
            xlabel('Step','FontSize',20,'FontName','Times');
        end
        xlim([0,(numel(curv)+1)])
        
        load(handles.param.result_file,'-mat','strain')
        TT=-1:0.01:1;
        id=3;
        fbeam(id)=figure(id);
        delete(gca)
        axes1 = axes('Parent',fbeam(id),'LineWidth',2,'FontSize',16,...
            'FontName','Times');
        set(fbeam(id),'Name',sprintf('Axial strain'));
        set(fbeam(id),'NumberTitle','off');
        box('on');
        hold('all');
        for iim=1:dt:size(fleche,2)
            plot(strain(:,iim),0.5*TT,'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
        end
        ylabel('Position t/h []','FontSize',20,'FontName','Times');
        xlabel('Longitudinal strain []','FontSize',20,'FontName','Times');
        %        if size(fleche,2)<6
        leg=legend(axes1,'show','Location','NorthEast');
        set(leg,'Box','off')
        %        end
        if handles.beam_model.exx
            load(handles.param.result_file,'-mat','naxis')
            
            id=4;
            fbeam(id)=figure(id);
            delete(gca)
            
            axes1 = axes('Parent',fbeam(id),'LineWidth',2,'FontSize',16,...
                'FontName','Times');
            set(fbeam(id),'Name',sprintf('Neutral axis position'));
            set(fbeam(id),'NumberTitle','off');
            box('on');
            hold('all');
            
            plot(naxis,'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
            
            ylabel('Neutral axis position t/h []','FontSize',20,'FontName','Times');
            xlabel('Step','FontSize',20,'FontName','Times');
            ylim([-0.5,0.5])
            xlim([0,(numel(naxis)+1)])
            
        end
    else
        id=3;
        fbeam(id)=figure(id);
        delete(gca)
        axes1 = axes('Parent',fbeam(id),'LineWidth',2,'FontSize',16,...
            'FontName','Times');
        set(fbeam(id),'Name',sprintf('Beam rotation'));
        set(fbeam(id),'NumberTitle','off');
        box('on');
        hold('all');
        for iim=1:dt:size(fleche,2)
            plot(s'*pix2m,rot(:,iim),'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
        end
        ylabel('Rotation []','FontSize',20,'FontName','Times');
        if pix2m==1
            xlabel('s [pixel]','FontSize',20,'FontName','Times');
        else
            xlabel('s [m]','FontSize',20,'FontName','Times');
        end
        %        if size(fleche,2)<6
        leg=legend(axes1,'show','Location','NorthEast');
        set(leg,'Box','off')
        %        end
        id=2;
        fbeam(id)=figure(id);
        delete(gca)
        axes1 = axes('Parent',fbeam(id),'LineWidth',2,'FontSize',16,...
            'FontName','Times');
        set(fbeam(id),'Name',sprintf('Beam curvature'));
        set(fbeam(id),'NumberTitle','off');
        box('on');
        hold('all');
        for iim=1:dt:size(fleche,2)
            plot(s'*pix2m,curv(:,iim)/pix2m,'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
        end
        if pix2m==1
            ylabel('Curvature [1/pixel]','FontSize',20,'FontName','Times');
            xlabel('s [pixel]','FontSize',20,'FontName','Times');
        else
            ylabel('Curvature [1/m]','FontSize',20,'FontName','Times');
            xlabel('s [m]','FontSize',20,'FontName','Times');
        end
        %        if size(fleche,2)<6
        leg=legend(axes1,'show','Location','NorthEast');
        set(leg,'Box','off')
        %        end
        if handles.beam_model.exx
            load(handles.param.result_file,'-mat','naxis')
            
            id=4;
            fbeam(id)=figure(id);
            delete(gca)
            
            axes1 = axes('Parent',fbeam(id),'LineWidth',2,'FontSize',16,...
                'FontName','Times');
            set(fbeam(id),'Name',sprintf('Neutral axis position'));
            set(fbeam(id),'NumberTitle','off');
            box('on');
            hold('all');
            
            for iim=1:dt:size(fleche,2)
                plot(s'*pix2m,naxis(:,iim),'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
            end
            
            ylabel('Neutral axis position t/h []','FontSize',20,'FontName','Times');
            ylim([-0.5,0.5])
            if pix2m==1
                xlabel('s [pixel]','FontSize',20,'FontName','Times');
            else
                xlabel('s [m]','FontSize',20,'FontName','Times');
            end
            %        if size(fleche,2)<6
            leg=legend(axes1,'show','Location','NorthEast');
            set(leg,'Box','off')
            
        end
        
    end
    set(0,'CurrentFigure',handles.figure1);
    
end

% --------------------------------------------------------------------
function cmenu_selected_point_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_selected_point_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try delete(handles.gpoint); catch, end
handles.gpoint=0;
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_selected_point_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_selected_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function zone_er_set(handles)
% hObject    handle to inclusion_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Click opposite corners of the rectangular zone.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    handles.initial_point=point1;
    guidata(handles.figure1,handles);
    xyp=[point1;point1];
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_rect);
    
    waitforbuttonpress
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    
    point2=getposition(handles);
    point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
    ly=abs((point1(1)-point2(1)));
    lx=abs((point1(2)-point2(2)));
    if strcmp(get(handles.figure1,'selectiontype'),'normal')&&(~(max(lx,ly)<1))
        delete(gz);
        roi=[max(1,min(handles.sizeim(1),sort(round([point1(1),point2(1)])))),max(1,min(handles.sizeim(2),sort(round([point1(2),point2(2)]))))];
        xyp=[roi([1,2,2,1,1])',roi([3,3,4,4,3])'];
        gz=plot(xyp(:,1),xyp(:,2),'Color',[0.9,0.9,0.2],'LineStyle','-','LineWidth',2,'ButtonDownFcn',@zone_adjust);
        gz.Tag=num2str(rand(1));
        handles.erzone{2,1}=xyp;
        handles.erzone{3,1}=gz;
        handles.erzone{1,1}=2;
        handles.erzone{4,1}=1;
        handles=plot_mesh(handles);
        handles=set_param(handles);
        set(handles.message_text,'String','Creating rectangular zone.....done')
    else
        delete(gz);
        set(handles.message_text,'String','Zone definition.....cancelled')
    end
else
    set(handles.message_text,'String','Zone definition.....cancelled')
end


% --------------------------------------------------------------------
function cmenu_export_error_delete_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_error_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try delete(handles.erzone{3,1}); catch, end
handles.erzone={};
handles.ercomp=1;
handles=plot_mesh(handles);
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_export_error_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
close all


% --------------------------------------------------------------------
% --------------------------------------------------------------------
function cmenu_zone_attractor_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_zone_attractor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='attractor_mesh_size';
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.active_zone=iz;
if (handles.fem_model.zone{6,iz})>0
    handles.scroll_value=handles.fem_model.zone{6,iz};
else
    handles.scroll_value=mean(handles.fem_model.mesh_size);
end
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust attractor mesh size.....0 to disable attraction')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);

function cmenu_mesh_attractors_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function cmenu_mesh_attractors_line_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractors_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
line_set(handles,1)



% --------------------------------------------------------------------
function cmenu_mesh_attractors_rectangle_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractors_rectangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zone_set(handles,1)


% --------------------------------------------------------------------
function cmenu_mesh_attractors_poly_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractors_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
poly_set(handles,1);


% --------------------------------------------------------------------
function cmenu_mesh_attractors_circle_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractors_circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
circle_set(handles,1)


% --------------------------------------------------------------------
function cmenu_attractor_define_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_attractor_define (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
set(handles.cmenu_attractor_define_txt,'Label',num2str(real(handles.fem_model.zone{6,iz})));

% --------------------------------------------------------------------
function cmenu_attractor_define_txt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_attractor_define_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='attractor_mesh_size';
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.active_zone=iz;
handles.scroll_value=handles.fem_model.zone{6,iz};
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust attractor mesh size.....0 to disable attraction')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);



% --------------------------------------------------------------------
function cmenu_attractor_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_attractor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_attractor_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_attractor_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.fem_model.zone(:,iz)=[];
delete(gco)
handles=remesh(handles);


% --------------------------------------------------------------------
function cmenu_attractor_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_attractor_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
if imag(handles.fem_model.zone{6,iz})
    handles.fem_model.zone{6,iz}=real(handles.fem_model.zone{6,iz});
else
    handles.fem_model.zone{6,iz}=handles.fem_model.zone{6,iz}+1i;
end
handles=remesh(handles);


% --------------------------------------------------------------------
function cmenu_mesh_attractor_rect_area_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractor_rect_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zone_set(handles,2)


% --------------------------------------------------------------------
function cmenu_mesh_attractor_poly_area_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractor_poly_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
poly_set(handles,2)


% --------------------------------------------------------------------
function cmenu_mesh_attractor_disc_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_attractor_disc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
circle_set(handles,2)


% --------------------------------------------------------------------
function cmenu_mesh_post_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_post (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_mesh_plot_line_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_plot_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
line_set(handles)


% --------------------------------------------------------------------
function switch_view_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to switch_view_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[az,el] = view;
dview=(az==0)&&(el==90);
if dview
    view(handles.view);
else
    handles.view=[az,el];
    view([0,90]);
end
handles=plot_mesh(handles);
handles=set_param(handles);

function update_display(hObject,eventdata)
handles=guidata(hObject);
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function rotate_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to rotate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'State')
    case 'on'
        h=rotate3d(handles.figure1);
        set(h,'ActionPostCallback',@update_display)
        set(h,'Enable','on')
    case 'off'
        rotate3d off
end


% --------------------------------------------------------------------
function switch_frame_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to switch_frame_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.pframe=~handles.pframe;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_mesh_factors_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_factors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_mesh_warp_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_warp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='warp';
handles.showim=0;
handles.ondefimage=1;
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust the amplification factor.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_mesh_elevation_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='elevation';
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust the amplification factor.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_field_data_calculator_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_calculator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_calculator_data_new_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_calculator_data_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer=inputdlg({'Label','Formula','Unit'},'Calculator',1,...
    {'Result','',''});
if~isempty(answer)
    handles.calculator{1,end+1}=answer{1};
    handles.calculator{2,end}=answer{2};
    newitem=uimenu(handles.cmenu_field_data_calculator,'Label',answer{1},'Callback',@calculator_callback);
    newitem.Tag=num2str(rand(1));
    handles.calculator{3,end}=newitem;
    handles.calculator{4,end}=answer{3};
    handles.field=5;
    handles.ccomp=size(handles.calculator,2)-1;
    handles=plot_mesh(handles);
    handles=set_param(handles);
end

function calculator_callback(hObject, eventdata)
handles=guidata(gcf);
iz=findg((handles.calculator(3,:)),hObject);
handles.field=5;
handles.ccomp=iz-1;
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function last_frame_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to last_frame_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.animation.iim=handles.animation.nbstep;
handles=display_frame(handles);


% --------------------------------------------------------------------
function next_frame_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to next_frame_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iim=handles.animation.iim;
handles.animation.iim=min(iim+1,handles.animation.nbstep);
handles=display_frame(handles);


% --------------------------------------------------------------------
function play_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to play_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.animation.playing
    pause_animation(handles);
else
    play_animation(handles);
end

% --------------------------------------------------------------------
function previous_frame_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to previous_frame_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iim=handles.animation.iim;
handles.animation.iim=max(iim-1,0);
handles=display_frame(handles);


% --------------------------------------------------------------------
function pprevious_frame_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pprevious_frame_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iim=handles.animation.iim;
dt=round(handles.animation.nbstep/10);
handles.animation.iim=max(iim-dt,0);
handles=display_frame(handles);


% --------------------------------------------------------------------
function nnext_frame_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to nnext_frame_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iim=handles.animation.iim;
dt=round(handles.animation.nbstep/10);
handles.animation.iim=min(iim+dt,handles.animation.nbstep);
handles=display_frame(handles);

% --------------------------------------------------------------------
function first_frame_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to first_frame_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.animation.iim=0;
handles=display_frame(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over image_info_text.
function image_info_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to image_info_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_mesh_crack_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_crack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=crack_set(handles);


% --------------------------------------------------------------------
function cmenu_mesh_crack_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_crack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function handles=crack_set(handles)
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Start defining the crack by picking its tip.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    xyp=[point1;point1];
    point2=-10000000000;
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(handles.message_text,'String','.....and now pick up points on its path, right click to end.....')
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_poly);
    
    while abs((point2-xyp(1,:))*[1;1i])>0.01*min(handles.sizeim(1:2))
        waitforbuttonpress
        if strcmp(get(handles.figure1,'selectiontype'),'normal')
            point2=getposition(handles);
            point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
            xyp=[xyp;point2;point2];
            set(gz,'Xdata',xyp(:,1));
            set(gz,'Ydata',xyp(:,2));
        else
            if size(xyp,1)<3
                delete(gz);
                set(gcf,'WindowButtonDownFcn',@double_clic);
                set(gcf,'WindowButtonMotionFcn',@show_position);
                return;
                set(handles.message_text,'String','Crack definition.....cancelled')
            else
                break;
            end
        end
        
    end
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    xyp=xyp(1:2:end,:);
    set(gz,'Xdata',xyp(:,1));
    set(gz,'Ydata',xyp(:,2));
    gz.Tag=num2str(rand(1));
    handles.fem_model.zone{2,end+1}=xyp;
    handles.fem_model.zone{3,end}=gz;
    handles.fem_model.zone{4,end}=5;
    handles.fem_model.zone{1,end}=-1;
    handles.fem_model.zone{6,end}=mean(handles.fem_model.mesh_size)+1i;
    handles.fem_model.zone{7,end}=mean(handles.fem_model.mesh_size);
    handles.fem_model.zone{9,end}=0;
    
    set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
    set(gz,'Color',[1,0,0]);
    if ~(handles.ana==2)
        set(gz,'uicontextmenu',handles.cmenu_crack);
        handles.fem_model.zone{8,end}=0;
        handles.fem_model.mesh_type=2;
        handles=remesh(handles);
    else
        handles.fem_model.zone{8,end}=-1;
        set(gz,'uicontextmenu',handles.cmenu_export_crack);
    end
    if (handles.ana==1)||(handles.ana==3)
        oxfem=handles.fem_model.phantom_nodes;
        xfem=menu('Crack model','FEM','X-FEM');
        handles.fem_model.phantom_nodes=xfem-1;
        if ~(oxfem==handles.fem_model.phantom_nodes)
            handles=remesh(handles);
        end
    end
    handles=set_param(handles);
    set(handles.message_text,'String','Creating crack.....done')
    
else
    set(handles.message_text,'String','Crack definition.....cancelled')
end
set(gcf,'WindowButtonDownFcn',@double_clic);
set(gcf,'WindowButtonMotionFcn',@show_position);


% --------------------------------------------------------------------
function cmenu_crack_tip_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_tip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
set(handles.cmenu_crack_tip_txt,'Label',num2str(real(handles.fem_model.zone{7,iz})));

% --------------------------------------------------------------------
function cmenu_crack_tip_txt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_tip_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='tip_mesh_size';
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.active_zone=iz;
handles.scroll_value=handles.fem_model.zone{7,iz};
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust crack tip mesh size.....0 to disable attraction')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);

% --------------------------------------------------------------------
function cmenu_crack_face_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_face (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
set(handles.cmenu_crack_face_txt,'Label',num2str(real(handles.fem_model.zone{6,iz})));

% --------------------------------------------------------------------
function cmenu_crack_face_txt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_face_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='face_mesh_size';
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.active_zone=iz;
handles.scroll_value=handles.fem_model.zone{6,iz};
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust face mesh size.....0 to disable attraction')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_crack_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.fem_model.zone(:,iz)=[];
delete(gco)
switch handles.ana
    case {1,3}
        handles=remesh(handles);
end

% --------------------------------------------------------------------
function cmenu_crack_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
if handles.fem_model.zone{9,iz}
    set(handles.cmenu_crack_ends,'Checked','on');
else
    set(handles.cmenu_crack_ends,'Checked','off');
end
if handles.param.detect
    set(handles.cmenu_crack_detection,'Checked','on');
else
    set(handles.cmenu_crack_detection,'Checked','off');
end
% --------------------------------------------------------------------
function cmenu_crack_sif_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_sif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_mesh_sif_wil_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_sif_wil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_crack_sif_size_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_sif_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
rc=handles.fem_model.zone{8,iz};
set(handles.cmenu_crack_sif_size_txt,'Label',num2str(rc(1)));


% --------------------------------------------------------------------
function cmenu_crack_sif_size_txt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_sif_size_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='sif_zone_size';
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.active_zone=iz;
rc=handles.fem_model.zone{8,iz};
handles.scroll_value=rc(1);
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust extraction zone size.....0 to disable extraction')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_crack_ends_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_ends (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.fem_model.zone{9,iz}=~handles.fem_model.zone{9,iz};
if ~(handles.ana==2)
    handles=remesh(handles);
end

handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_export_crack_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_crack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
zone=handles.fem_model.zone(:,iz);
set(handles.cmenu_crack_export_remove,'Visible','off')
zsize=zone{8};
if abs(zsize(1))>0&&~iscell(handles.uvisu)
    if zsize(1)<0
        set(handles.cmenu_crack_export_remove,'Visible','on')
    end
end

% --------------------------------------------------------------------
function cmenu_crack_export_plot_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_export_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
plot_crack_data(handles,iz,0,1)
% --------------------------------------------------------------------
function cmenu_crack_export_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_export_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.fem_model.zone(:,iz)=[];
try delete(handles.fgage(iz));catch , end
try delete(100+handles.fgage(iz));catch , end
try delete(200+handles.fgage(iz));catch , end
try delete(300+handles.fgage(iz));catch , end
try delete(4000+handles.fgage(iz));catch , end
try delete(400+handles.fgage(iz));catch , end
try delete(500+handles.fgage(iz));catch , end
try delete(600+handles.fgage(iz));catch , end
delete(gco)
guidata(hObject,handles)

% --------------------------------------------------------------------
function cmenu_crack_export_export_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_export_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
plot_crack_data(handles,iz,1,1)


function plot_crack_data(handles,id,export,dok)
if nargin<4,dok=1;end
if nargin<3,export=0;end
if handles.preview
    handles.showplot=1;
    handles.fgage(id)=id;
    guidata(handles.figure1,handles);
    zone=handles.fem_model.zone(:,id);
    if zone{8}(1)<0||handles.fem_model.phantom_nodes==1
        %        display('NOT CODED YET')
    else
        indc=zone{10};
        nodes=indc{1};
        Nnodes=handles.mvisu.Nnodes;
        iim= handles.animation.iim;
        xon=handles.mvisu.xo(nodes(:,1));
        yon=handles.mvisu.yo(nodes(:,1));
        s=[0;cumsum(abs(diff(xon+1i*yon)))];
        t=-gradient(xon+1i*yon);
        t=t./abs(t);
        if (~iscell(handles.uvisu))
            if iim
                if handles.stereo
                    U=handles.uxyz(:,iim);
                    U=[reshape(U,prod(Nnodes),3)];
                else
                    if iscell(handles.uvisu)
                        U=handles.uvisu{iim};
                    else
                        U=handles.uvisu(:,iim);
                    end
                    U=reshape(U,prod(Nnodes),2);
                end
                nc=size(U,2);
                U=U(nodes(:),:);
                U=reshape(U(:),size(nodes,1),2*nc);
                if nc==2
                    todu=[1,0;-1,0;0,1;0,-1];
                else
                    todu=[1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1];
                end
                
                U=U*todu;
                Uc=U(:,1)+1i*U(:,2);
                U(:,1)=real((Uc.').*(t'));
                U(:,2)=real((Uc.').*((t*exp(1i*pi/2))'));
                
            else
                U=zeros(length(s),3);
            end
            if any(diff(s)>3*zone{7})
                ids=find(diff(s)>3*zone{7});
                for ii=1:length(ids)
                    t(ids(ii))=-diff(xon(ids(ii)+(-1:0))+1i*yon(ids(ii)+(-1:0)));
                    t(ids(ii)+1)=-diff(xon(ids(ii)+(0:1))+1i*yon(ids(ii)+(0:1)));
                    s=[s;0.5*(s(ids(ii)+1)+s(ids(ii)))];
                    U=[U;repmat(NaN,1,size(U,2))];
                    t=[t;NaN];
                    
                end
                t=t./abs(t);
                [s,ids]=sort(s);
                U=U(ids,:);
                t=t(ids);
            end
            fgage=figure(400+id);
            if ~isempty(findobj(gca,'DisplayName','dU_x'))
                set(findobj(gca,'DisplayName','dU_t'),'Ydata',U(:,1));
                set(findobj(gca,'DisplayName','dU_t'),'Xdata',s);
                set(findobj(gca,'DisplayName','dU_n'),'Ydata',U(:,2));
                set(findobj(gca,'DisplayName','dU_n'),'Xdata',s);
                if ~isempty(findobj(gca,'DisplayName','dU_z'))
                    set(findobj(gca,'DisplayName','dU_z'),'Ydata',U(:,3));
                    set(findobj(gca,'DisplayName','dU_z'),'Xdata',s);
                end
                set(fgage,'Name',sprintf('Crack #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
            else
                delete(gca)
                axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                    'FontName','Times');
                set(fgage,'Name',sprintf('Crack #%d: step %d/%d',id,handles.animation.iim,handles.animation.nbstep));
                set(fgage,'NumberTitle','off');
                box('on');
                hold('all');
                plot(s,U(:,1),'DisplayName','dU_t','Parent',axes1,'LineWidth',2);
                plot(s,U(:,2),'DisplayName','dU_n','Parent',axes1,'LineWidth',2);
                if handles.stereo
                    plot(s,Uxyz(:,3),'DisplayName','dU_z','Parent',axes1,'LineWidth',2);
                    ylabel('Displacement jump [m]','FontSize',20,'FontName','Times');
                else
                    ylabel('Displacement jump [pixel]','FontSize',20,'FontName','Times');
                end
                xlabel('s [pixel]','FontSize',20,'FontName','Times');
                leg=legend(axes1,'show','Location','NorthWest');
                set(leg,'Box','off')
            end
            
            if export
                if handles.stereo
                    U=handles.uxyz;
                else
                    U=handles.uvisu;
                end
                Ux=U(nodes(:,1),:)-U(nodes(:,2),:);
                Uy=U(nodes(:,1)+prod(Nnodes),:)-U(nodes(:,2)+prod(Nnodes),:);
                if handles.stereo
                    Uz=U(nodes(:,1)+2*prod(Nnodes),:)-U(nodes(:,2)+2*prod(Nnodes),:);
                end
                s(isnan(t))=[];
                [pp,filename,ext]=fileparts(handles.param.result_file);
                filexp=sprintf('%s-crack-%02d-du.csv',filename,id);
                if exist(filexp,'file')
                    set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                    [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save gage data...',filexp);
                end
                
                fid=fopen(filexp,'w');
                fprintf(fid,'Result file;%s\n',handles.param.result_file);
                if ~strcmp(handles.param.analysis,'mechanics')
                    fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
                end
                fprintf(fid,'Crack;%d\n',id);
                if handles.stereo
                    if isfield(handles.param,'deformed_image')
                        fprintf(fid,'Filename;;;');
                        for iim=1:size(Uy,2)
                            if iscell(handles.param.deformed_image)
                                fprintf(fid,'"%s";;;',handles.param.deformed_image{1,iim+handles.stereo});
                            else
                                fprintf(fid,'"%s";;;',handles.param.deformed_image);
                            end
                        end
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,'Step;;;');
                    for iim=1:size(Uy,2)
                        fprintf(fid,'%d;;;',iim);
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'x [pixel];y [pixel];s [pixel];Ux[m];Uy[m];Uz[m]\n');
                    for ip=1:length(s)
                        fprintf(fid,'%.3e;%.3e;%.3e;',xon(ip),yon(ip),s(ip));
                        for iim=1:size(Uy,2)
                            fprintf(fid,'%.3e;%.3e;%.3e;',Ux(ip,iim),Uy(ip,iim),Uz(ip,iim));
                        end
                        fprintf(fid,'\n');
                    end
                    
                else
                    if isfield(handles.param,'deformed_image')
                        fprintf(fid,'Filename;;;');
                        for iim=1:size(Uy,2)
                            if iscell(handles.param.deformed_image)
                                fprintf(fid,'"%s";;',handles.param.deformed_image{1,iim+handles.stereo});
                            else
                                fprintf(fid,'"%s";;',handles.param.deformed_image);
                            end
                        end
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,'Step;;;');
                    for iim=1:size(Uy,2)
                        fprintf(fid,'%d;;',iim);
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'x [pixel];y [pixel];s [pixel];Ux[pixel];Uy[pixel];\n');
                    for ip=1:length(s)
                        fprintf(fid,'%.3e;%.3e;%.3e;',xon(ip),yon(ip),s(ip));
                        for iim=1:size(Uy,2)
                            fprintf(fid,'%.3e;%.3e;',Ux(ip,iim),Uy(ip,iim));
                        end
                        fprintf(fid,'\n');
                    end
                end
                fclose(fid);
                
                
                param=handles.param;
                model=handles.fem_model;
                save(sprintf('%s-crack-%02d-du.res',filename,id),'Ux','Uy','zone','param','model','-v7.3')
                if handles.stereo
                    save(sprintf('%s-crack-%02d-du.res',filename,id),'Uz','-append')
                end
                set(handles.message_text,'String','Exporting crack data.....done')
                
                
            end
        end
    end
    zsize=zone{8};
    if abs(zsize(1))>0
        if zsize(1)<0
            if isempty(zone{5})
                
                if ~isfield(handles.fem_model,'material_parameters')
                    handles.fem_model.material_parameters.young=2e11;
                    handles.fem_model.material_parameters.nu=0.3;
                    handles.fem_model.material_parameters.mu=0.5*2e11/(1+.3);
                    handles.fem_model.material_parameters.kappa=(3-4*nu);
                end
                if ~isfield(handles.fem_model.material_parameters,'young')
                    handles.fem_model.material_parameters.young=2e11;
                    handles.fem_model.material_parameters.nu=0.3;
                end
                if ~isfield(handles.fem_model.material_parameters,'kappa')
                    handles.fem_model.material_parameters.mu=...
                        0.5*handles.fem_model.material_parameters.young/(1+handles.fem_model.material_parameters.nu);
                    handles.fem_model.material_parameters.kappa=3-4*handles.fem_model.material_parameters.nu;
                end
                cparam.mu=handles.fem_model.material_parameters.mu;
                cparam.kappa=handles.fem_model.material_parameters.kappa;
                cparam.pix2m=handles.param.pixel_size;
                cparam.steps=1:size(handles.uvisu,2);
                cparam.radius=round(10*handles.fem_model.mesh_size(1));
                cparam.mask_radius=round(1.5*handles.fem_model.mesh_size(1));
                cparam.mask_width=round(1.5*handles.fem_model.mesh_size(1));
                cparam.check=0;
                cparam.xyc=zone{2};
                go=1;
            else
                cparam=zone{5};
                if isfield(cparam,'young')&~isfield(cparam,'kappa')
                    cparam.mu=0.5*cparam.young/(1+cparam.nu);
                    cparam.kappa=3-4*cparam.nu;
                    cparam=rmfield(cparam,'young');
                    cparam=rmfield(cparam,'nu');
                    
                end
            end
            if ~iscell(handles.uvisu)
                if dok
                    thi=0:pi/50:2*pi;
                    gmax=plot(cparam.xyc(1,:)*[1;1i]+cparam.radius*exp(1i*thi),'r-','LineWidth',2);
                    gmin=plot(cparam.xyc(1,:)*[1;1i]+cparam.mask_radius*exp(1i*thi),'r-','LineWidth',2);
                    rmax=max(handles.sizeim);
                    iz=1;
                    for ri=50:100:rmax
                        gc{iz}=plot(cparam.xyc(1,:)*[1;1i]+ri*exp(1i*thi),'r--','LineWidth',1);
                        iz=iz+1;
                    end
                    answer=inputdlg({'mu','kappa','Steps','Rmax','Rmin','pixel_size','check'},sprintf('Extraction parameters'),1,...
                        {num2str(cparam.mu),num2str(cparam.kappa),sprintf('%d:%d',min(cparam.steps),max(cparam.steps)),num2str(cparam.radius),num2str(cparam.mask_radius),num2str(cparam.pix2m),num2str(cparam.check)});
                    if isempty(answer)
                        return
                    end
                    delete(gmax)
                    delete(gmin)
                    for iz=1:length(gc)
                        delete(gc{iz})
                    end
                    go=0;
                    if ~(cparam.mu==eval(answer{1}))
                        cparam.mu=eval(answer{1});
                        go=1;
                    end
                    if ~(cparam.kappa==eval(answer{2}))
                        cparam.kappa=eval(answer{2});
                        go=1;
                    end
                    steps=eval(answer{3});
                    if (~(min(cparam.steps)==min(steps)))||(~(max(cparam.steps)==max(steps)))
                        cparam.steps=steps;
                        go=1;
                    end
                    if ~(cparam.radius==eval(answer{4}))
                        cparam.radius=eval(answer{4});
                        go=1;
                    end
                    if ~(cparam.mask_radius==eval(answer{5}))
                        cparam.mask_radius=eval(answer{5});
                        cparam.mask_width=eval(answer{5});
                        go=1;
                    end
                    if ~(cparam.pix2m==eval(answer{6}))
                        cparam.pix2m=eval(answer{6});
                        go=1;
                    end
                    if ~(cparam.check==eval(answer{7}))
                        cparam.check=eval(answer{7});
                        go=1;
                    end
                    try
                        if any(~(cparam.xyc==zone{2}))
                            go=1;
                        end
                    catch
                        go=1;
                    end
                    cparam.xyc=zone{2};
                    zone{5}=cparam;
                    if isempty(zone{5})||go
                        run_williams_curved(id,zone,cparam,handles.param.result_file);
                    end
                    handles.fem_model.zone(:,id)=zone;
                    handles=set_param(handles);
                end
                load(sprintf('%s-crack-%02d-sif.res',strrep(handles.param.result_file,'.res',''),id),'-mat','da','xytips','K1','K2','T','B','param','model')
                images=cparam.steps;
            else
                load(handles.param.result_file,'-mat','da','xytips','K1','K2','T','B','param','model')
                handles.cracked=2;
                images=1:numel(K1);
                guidata(handles.figure1,handles);
            end
            detect=1;
            pix2m=cparam.pix2m;
        else
            load(handles.param.result_file,'-mat','Ks','tips','cracks','U')
            mu=0.5*handles.fem_model.material_parameters.young/(1+handles.fem_model.material_parameters.nu);
            Es=handles.fem_model.material_parameters.young/(1-(handles.fem_model.material_parameters.nu)^2);
            pix2m=handles.param.pixel_size;
            modes=[1,2];
            detect=handles.param.detect;
            if detect==1
                ind=-3:7;
            else
                ind=0:7;
            end
            xfem=handles.fem_model.phantom_nodes;
            if xfem==1
                ind=1:4;
            end
            scal=2*mu*sqrt(2*pi);
            nw=length(modes)*length(ind);
            scalamp=scal*(pix2m.^(1-(ind)*.5));
            T=[];B=[];
            if xfem&&~iscell(U)
                K1=sqrt(pix2m)*Ks(2*(cracks(id)-1)+1,1);
                K2=sqrt(pix2m)*Ks(2*(cracks(id)-1)+1,2);
                if size(Ks,2)>3
                    T=Ks(2*(cracks(id)-1)+1,6)/(4/sqrt(2*pi));
                    B=Ks(2*(cracks(id)-1)+1,7)/sqrt(pix2m);
                end
            else
                found=find(ind==1)+nw*(tips(cracks==id,1)-1);
                if xfem
                    K1=[sqrt(pix2m)*Ks(found,:)];
                    K2=[sqrt(pix2m)*Ks(found+1,:)];
                    T=Ks(found+5,:)/(4/sqrt(2*pi));
                    B=Ks(found+6,:)/sqrt(pix2m);
                else
                    K1=[scalamp(found)*Ks(found,:)];
                    K2=[scalamp(found)*Ks(found+length(ind),:)];
                end
            end
            found=find(ind==2);
            %        T=[scalamp(found)*Ks(found,:)];
            if isempty(T)
                T=[scalamp(found)*Ks(found,:)]*4/sqrt(2*pi);
            end
            found=find(ind==3);
            if isempty(B)
                B=[scalamp(found)*Ks(found,:)];
            end
            images=(1:length(K1));
            if detect==1
                
                load(handles.param.result_file,'-mat','da','ztips')
                path=zone{2}*[1;1i];
                s=[0;cumsum(abs(diff(path)))];
                xytips=interp1(s,path,ztips(tips(cracks==id,1),:),'linear','extrap');
                da=da(tips(cracks==id,1),:)*pix2m;
            end
            param=handles.param;
            model=handles.fem_model;
        end
        
        fgage=figure(100+id);
        if ~isempty(findobj(gca,'DisplayName','K_{I}'))
            set(findobj(gca,'DisplayName','K_{I}'),'Ydata',K1);
            set(findobj(gca,'DisplayName','K_{II}'),'Ydata',K2);
            set(findobj(gca,'DisplayName','K_{I}'),'Xdata',images);
            set(findobj(gca,'DisplayName','K_{II}'),'Xdata',images);
            set(fgage,'Name',sprintf('SIF Crack #%d',id));
        else
            delete(gca)
            axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                'FontName','Times');
            set(fgage,'Name',sprintf('SIF Crack #%d',id));
            set(fgage,'NumberTitle','off');
            box('on');
            hold('all');
            if numel(K1)==1
                plot(images,K1,'x','DisplayName','K_{I}','Parent',axes1,'LineWidth',2);
                plot(images,K2,'x','DisplayName','K_{II}','Parent',axes1,'LineWidth',2);
            else
                plot(images,K1,'DisplayName','K_{I}','Parent',axes1,'LineWidth',2);
                plot(images,K2,'DisplayName','K_{II}','Parent',axes1,'LineWidth',2);
            end
            ylabel('SIF','FontSize',20,'FontName','Times');
            xlabel('Step','FontSize',20,'FontName','Times');
            leg=legend(axes1,'show','Location','NorthWest');
            set(leg,'Box','off')
        end
        if ~isempty(T)
            fgage=figure(200+id);
            if ~isempty(findobj(gca,'DisplayName','T'))
                set(findobj(gca,'DisplayName','T'),'Ydata',T);
                set(findobj(gca,'DisplayName','T'),'Xdata',images);
                set(fgage,'Name',sprintf('T-stress Crack #%d',id));
            else
                delete(gca)
                axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                    'FontName','Times');
                set(fgage,'Name',sprintf('T-stress Crack #%d',id));
                set(fgage,'NumberTitle','off');
                box('on');
                hold('all');
                if numel(T)==1
                    plot(images,T,'x','DisplayName','T','Parent',axes1,'LineWidth',2);
                else
                    plot(images,T,'DisplayName','T','Parent',axes1,'LineWidth',2);
                end
                ylabel('T-stress','FontSize',20,'FontName','Times');
                xlabel('Step','FontSize',20,'FontName','Times');
            end
        end
        if ~isempty(B)
            fgage=figure(300+id);
            if ~isempty(findobj(gca,'DisplayName','B'))
                set(findobj(gca,'DisplayName','B'),'Ydata',B);
                set(findobj(gca,'DisplayName','B'),'Xdata',images);
                set(fgage,'Name',sprintf('B-stress Crack #%d',id));
            else
                delete(gca)
                axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                    'FontName','Times');
                set(fgage,'Name',sprintf('B-stress Crack #%d',id));
                set(fgage,'NumberTitle','off');
                box('on');
                hold('all');
                if numel(B)==1
                    plot(images,B,'x','DisplayName','B','Parent',axes1,'LineWidth',2);
                else
                    plot(images,B,'DisplayName','B','Parent',axes1,'LineWidth',2);
                end
                ylabel('B-stress','FontSize',20,'FontName','Times');
                xlabel('Step','FontSize',20,'FontName','Times');
            end
        end
        if detect==1
            
            fgage=figure(400+id);
            if ~isempty(findobj(gca,'DisplayName','a'))
                set(findobj(gca,'DisplayName','a'),'Ydata',da);
                set(findobj(gca,'DisplayName','a'),'Xdata',images);
                set(fgage,'Name',sprintf('Crack length Crack #%d',id));
            else
                delete(gca)
                axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                    'FontName','Times');
                set(fgage,'Name',sprintf('Crack length Crack #%d',id));
                set(fgage,'NumberTitle','off');
                box('on');
                hold('all');
                if numel(da)==1
                    plot(images,da,'x','DisplayName','a','Parent',axes1,'LineWidth',2);
                else
                    plot(images,da,'DisplayName','a','Parent',axes1,'LineWidth',2);
                end
                ylabel('Crack length [m]','FontSize',20,'FontName','Times');
                xlabel('Step','FontSize',20,'FontName','Times');
            end
            
            fgage=figure(4000+id);
            delete(gca)
            axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
                'FontName','Times');
            set(fgage,'Name',sprintf('Crack length v.s. KI Crack #%d',id));
            set(fgage,'NumberTitle','off');
            box('on');
            hold('all');
            if numel(da)==1
                plot(da,K1,'x','DisplayName','a','Parent',axes1,'LineWidth',2);
            else
                plot(da,K1,'DisplayName','a','Parent',axes1,'LineWidth',2);
            end
            xlabel('Crack length [m]','FontSize',20,'FontName','Times');
            ylabel('SIF','FontSize',20,'FontName','Times');
            
            
            
            
            
            
        end
        
        
        if export
            
            [pp,filename,ext]=fileparts(handles.param.result_file);
            filexp=sprintf('%s-crack-%02d-sif.csv',filename,id);
            if exist(filexp,'file')
                set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
                [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save crack data...',filexp);
            end
            
            fid=fopen(filexp,'w');
            fprintf(fid,'Result file;%s\n',handles.param.result_file);
            if ~strcmp(handles.param.analysis,'mechanics')
                fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
            end
            fprintf(fid,'Crack;%d\n',id);
            dic=1;
            if isfield(handles.param,'deformed_image')
                if isempty(handles.param.deformed_image)
                    dic=0;
                end
                if ~iscell(handles.param.deformed_image)
                    handles.param.deformed_image={handles.param.deformed_image};          
                end
            end
            if dic
                if detect==1
                    fprintf(fid,'Filename;Step;a [m];xtip [pixel];ytip [pixel];K1 [Pa.sqrt(m)];K2 [Pa.sqrt(m)];T [Pa];B [Pa/sqrt(m)]\n');
                    for iim=1:length(K1)
                        fprintf(fid,'"%s";%d;%.4e;%.4f;%.4f;%.3e;%.3e;%.3e;%.3e\n',handles.param.deformed_image{1,images(iim)+handles.stereo},images(iim),da(iim),real(xytips(iim)),imag(xytips(iim)),K1(iim),K2(iim),T(iim),B(iim));
                    end
                else
                    if ~isempty(T)
                        fprintf(fid,'Filename;Step;K1 [Pa.sqrt(m)];K2 [Pa.sqrt(m)];T [Pa];B [Pa/sqrt(m)]\n');
                        for iim=1:length(K1)
                            fprintf(fid,'"%s";%d;%.3e;%.3e;%.3e;%.3e\n',handles.param.deformed_image{1,images(iim)+handles.stereo},images(iim),K1(iim),K2(iim),T(iim),B(iim));
                        end
                    else
                        fprintf(fid,'Filename;Step;K1 [Pa.sqrt(m)];K2 [Pa.sqrt(m)]\n');
                        for iim=1:length(K1)
                            fprintf(fid,'"%s";%d;%.3e;%.3e\n',handles.param.deformed_image{1,images(iim)+handles.stereo},images(iim),K1(iim),K2(iim));
                        end
                        
                    end
                end
            else
                if ~isempty(T)
                    fprintf(fid,'Step;K1 [];K2 [];T [];B []\n');
                    for iim=1:length(K1)
                        fprintf(fid,'%d;%.3e;%.3e;%.3e;%.3e\n',images(iim),K1(iim),K2(iim),T(iim),B(iim));
                    end
                else
                    fprintf(fid,'Step;K1 [];K2 []\n');
                    for iim=1:length(K1)
                        fprintf(fid,'%d;%.3e;%.3e\n',images(iim),K1(iim),K2(iim));
                    end
                end
            end
            
            fclose(fid);
            if zsize(1)>=0
                save(sprintf('%s-crack-%02d-sif.res',filename,id),'images','K1','K2','T','B','zone','param','model','-v7.3')
                if detect==1
                    save(sprintf('%s-crack-%02d-sif.res',filename,id),'da','xytips','-append')
                end
            end
        end
        
    end
end

% --------------------------------------------------------------------
function cmenu_mesh_smoothing_median_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_smoothing_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.param.regularization_type='median';
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_mesh_smoothing_strain_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_smoothing_strain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.param.regularization_type='tiko';
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_mesh_smoothing_stress_Callback(~, eventdata, handles)
% hObject    handle to cmenu_mesh_smoothing_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.param.regularization_type='equilibrium_gap';
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_bcs_adjust_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_bcs_adjust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
xyp=handles.fem_model.zone{2,iz};
set(gco,'LineStyle','--','Marker','none','ButtonDownFcn',@zone_adjust)
set(gco,'Xdata',xyp(:,1),'Ydata',xyp(:,2))

% --------------------------------------------------------------------
function cmenu_bcs_remove_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_bcs_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
handles.fem_model.zone(:,iz)=[];
delete(gco)
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_bcs_load_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_bcs_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
loads=handles.fem_model.zone{7,iz};
answer=inputdlg({'Fx','Fy','Ux','Uy'},sprintf('Load for node set %d',iz),1,...
    {num2str(loads(1,2)),num2str(loads(2,2)),num2str(loads(1,1)),num2str(loads(2,1))});
if ~isempty(answer)
    loads(1,1)=eval(answer{3});
    loads(2,1)=eval(answer{4});
    loads(1,2)=eval(answer{1});
    loads(2,2)=eval(answer{2});
    handles.fem_model.zone{7,iz}=loads;
    handles=plot_mesh(handles);
    handles=set_param(handles);
end
% --------------------------------------------------------------------
function cmenu_bcs_zone_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_bcs_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function show_data_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to show_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showplot=~(handles.showplot);
if ~(handles.showplot)
    for id=1:length(handles.fgage)
        try delete(handles.fgage(id));catch, end
        try delete(100+handles.fgage(id));catch, end
        try delete(200+handles.fgage(id));catch, end
        try delete(300+handles.fgage(id));catch, end
        try delete(400+handles.fgage(id));catch, end
        try delete(4000+handles.fgage(id));catch, end
        try delete(500+handles.fgage(id));catch, end
        try delete(600+handles.fgage(id));catch, end
    end
    for id=1:length(handles.fbeam)
        try delete(handles.fbeam(id));catch, end
    end
else
    if (handles.ana==2)
        switch handles.fbasis
            case 'beam'
                cmenu_export_beam_plot_Callback(handles.figure1, [], handles)
                handles=guidata(handles.figure1);
            case 'vic'
                cmenu_export_vic_plot_Callback(handles.figure1, [], handles)
                handles=guidata(handles.figure1);
        end
    end
end
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_mesh_bcs_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_bcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message_text,'String','Click opposite corners to add a node set with boundary conditions.....')
handles=bcs_zone_set(handles);
set(handles.message_text,'String','Boundary conditions.....done')


% --------------------------------------------------------------------
function cmenu_mesh_material_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer=inputdlg({'E','nu'},sprintf('Elastic material parameters'),1,...
    {sprintf('%g',handles.fem_model.material_parameters.young),num2str(handles.fem_model.material_parameters.nu)});
if ~isempty(answer)
    handles.fem_model.material_parameters.nu=eval(answer{2});
    handles.fem_model.material_parameters.young=eval(answer{1});
    guidata(handles.figure1,handles);
end


% --------------------------------------------------------------------
function cmenu_animation_images_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles.param,'deformed_image')
    [fil0,path0,FilterIndex]=uigetfile({...
        '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2;*.raw','Image files'},'Open file(s)','MultiSelect','on');
    if FilterIndex
        fils=SortImageFiles(fil0);
        fil0=fils;
        iim=handles.animation.iim;
        handles.animation.nbstep=size(fil0,2);
        handles.animation.frames=1:size(fil0,2);
        handles.animation.iim=min(iim,handles.animation.nbstep);
        if numel(fil0)>1
            handles.param.deformed_image=fil0;
        else
            handles.param.deformed_image=fil0{1};
        end
        set(handles.message_text,'String','Loading images.....done')
        handles=set_param(handles);
        handles=display_frame(handles);
        
    end
else
    answer=inputdlg({'Sampling','Last frame'},sprintf('Video'),1,...
        {sprintf('%g',handles.param.video_sampling),num2str(handles.param.number_of_frames)});
    if ~isempty(answer)
        frq=eval(answer{1});
        handles.param.video_sampling=frq;
        nbf=eval(answer{2});
        nbf=min(nbf,handles.reader.NumberOfFrames);
        handles.param.number_of_frames=nbf;
        frames=2:handles.param.video_sampling:nbf;
        handles.animation.frames=frames;
        handles.animation.nbstep=length(frames);
        handles=set_param(handles);
        
    end
    
end

% --------------------------------------------------------------------
function cmenu_crack_propagation_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_propagation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer=inputdlg({'KIc','da'},sprintf('Crack propagation parameters'),1,...
    {sprintf('%g',handles.fem_model.material_parameters.KIc),num2str(handles.param.da)});
if ~isempty(answer)
    handles.fem_model.material_parameters.KIc=eval(answer{1});
    handles.param.da=eval(answer{2});
    guidata(handles.figure1,handles);
end


% --------------------------------------------------------------------
function cmenu_export_bcs_plot_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_bcs_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.fem_model.zone(3,:)),gco);
plot_bcs_data(handles,iz,0)


% --------------------------------------------------------------------
function cmenu_export_bcs_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_bcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function plot_bcs_data(handles,id,export)
if nargin<3,export=0;end
if handles.preview
    handles.showplot=1;
    handles.fgage(id)=id;
    guidata(handles.figure1,handles);
    zone=handles.fem_model.zone(:,id);
    fu=zone{7};
    fu=[zeros(4,1),fu];
    
    images=(1:length(fu))-1;
    fgage=figure(100+id);
    if isempty(findobj(gca,'DisplayName','F_x'))
        delete(gca)
        axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
            'FontName','Times');
        set(fgage,'Name',sprintf('Load #%d',id));
        set(fgage,'NumberTitle','off');
        box('on');
        hold('all');
        if any(sum(abs(fu(1:2,:)),1))
            plot(fu(1,:),fu(3,:),'DisplayName','F_x','Parent',axes1,'LineWidth',2);
            plot(fu(2,:),fu(4,:),'DisplayName','F_y','Parent',axes1,'LineWidth',2);
            ylabel('F','FontSize',20,'FontName','Times');
            xlabel('U','FontSize',20,'FontName','Times');
        else
            plot(fu(3,:),'DisplayName','F_x','Parent',axes1,'LineWidth',2);
            plot(fu(4,:),'DisplayName','F_y','Parent',axes1,'LineWidth',2);
            ylabel('F','FontSize',20,'FontName','Times');
            xlabel('Step','FontSize',20,'FontName','Times');
        end
        leg=legend(axes1,'show','Location','NorthWest');
        set(leg,'Box','off')
    end
end


% --------------------------------------------------------------------
function new_calibration_button_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to new_calibration_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fil0,path0,FilterIndex]=uigetfile({...
    '*.BMP;*.bmp;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.CR2;*.raw','Image files'...
    },'Image file(s)','MultiSelect', 'on');
filc=Ucalib(path0,fil0);


% --------------------------------------------------------------------
function cmenu_crack_detection_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_crack_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.param.detect=~(handles.param.detect);
iz=findg((handles.fem_model.zone(3,:)),gco);
rc=handles.fem_model.zone{8,iz};
if handles.param.detect==1
    handles.fem_model.zone{8,iz}=[rc(1),20];
else
    handles.fem_model.zone{8,iz}=rc(1);
end
if handles.fem_model.mesh_type==2
    handles=remesh(handles);
end
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_animation_images_vtk_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_images_vtk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExportImageToVTK(handles.param.result_file);
set(handles.message_text,'String','Exporting images to VTK.....done')


% --------------------------------------------------------------------
function cmenu_mesh_colorscale_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_mesh_colorscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cmenu_colorbar_user,'Checked','off')
set(handles.cmenu_colorbar_minmax,'Checked','off')
set(handles.cmenu_colorbar_auto,'Checked','off')
switch handles.field
    case 1
        comp=handles.ucomp;
    case 2
        comp=handles.ecomp;
    case {4,6}
        comp=0;
    case 3
        comp=handles.ercomp;
    case 5
        comp=handles.ccomp;
end
switch handles.scale_mode(handles.field,comp+1)
    case 1
        set(handles.cmenu_colorbar_auto,'Checked','on')
    case 2
        set(handles.cmenu_colorbar_user,'Checked','on')
    case 3
        set(handles.cmenu_colorbar_minmax,'Checked','on')
end


% --------------------------------------------------------------------
function cmenu_animation_field_vtk_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_field_vtk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.stereo
    postproVTK25D([handles.param.result_file,''],0);
else
    postproVTK([handles.param.result_file,''],0,0);
end

% --------------------------------------------------------------------
function cmenu_animation_export_field_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_animation_export_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[pp,filename,ext]=fileparts(handles.param.result_file);
handles=ComputeStrain(handles);
iim=handles.animation.iim;
filexp=sprintf('%s-data-%04d.csv',filename,iim);
if exist(filexp,'file')
    set(handles.message_text,'String',sprintf('%s already exists.....',filexp))
    [filexp,pp]=uiputfile({'*.csv','CSV ascii file (*.csv)'},'Save field data...',filexp);
end

fid=fopen(filexp,'w');
fprintf(fid,'Result file;%s\n',handles.param.result_file);
fprintf(fid,'Reference image;%s\n',handles.param.reference_image);
if isfield(handles.param,'deformed_image')
    if iscell(handles.param.deformed_image)
        if handles.stereo
            fprintf(fid,'Deformed images;%s;%s\n',handles.param.deformed_image{1,iim+handles.stereo},handles.param.deformed_image{2,iim+handles.stereo});
        else
            fprintf(fid,'Deformed image;%s\n',handles.param.deformed_image{iim});
        end
    else
        fprintf(fid,'Deformed image;%s\n',handles.param.deformed_image);
    end
end
fprintf(fid,'Step;%d\n',iim);
Nnodes=handles.mvisu.Nnodes;
if handles.stereo
    U=handles.uxyz;
    xo=handles.mvisu.xo;
    yo=handles.mvisu.yo;
    Xo=handles.mvisu.Xo;
    Yo=handles.mvisu.Yo;
    Zo=handles.mvisu.Zo;
    if handles.sonelt
        fprintf(fid,'x [pixel];y [pixel];X [m];Y [m];Z [m];UX[m];UY[m];UZ[m]\n');
        for ip=1:length(xo)
            fprintf(fid,'%.3f;%.3f;%.3e;%.3e;%.3e;',handles.param.roi(1)-1+xo(ip),handles.param.roi(3)-1+yo(ip),Xo(ip),Yo(ip),Zo(ip));
            fprintf(fid,'%.3e;%.3e;%.3e;',U(ip,iim),U(ip+prod(Nnodes),iim),U(ip+2*prod(Nnodes),iim));
            fprintf(fid,'\n');
        end
    else
        fprintf(fid,'x [pixel];y [pixel];X [m];Y [m];Z [m];UX[m];UY[m];UZ[m];dUXdX[];dUXdY[];dUYdX[];dUYdY[]\n');
        for ip=1:length(xo)
            fprintf(fid,'%.3f;%.3f;%.3e;%.3e;%.3e;',handles.param.roi(1)-1+xo(ip),handles.param.roi(3)-1+yo(ip),Xo(ip),Yo(ip),Zo(ip));
            fprintf(fid,'%.3e;%.3e;%.3e;',U(ip,iim),U(ip+prod(Nnodes),iim),U(ip+2*prod(Nnodes),iim));
            fprintf(fid,'%.3e;%.3e;%.3e;%.3e;',handles.evisu.xx(ip),handles.evisu.xy(ip),handles.evisu.yx(ip),handles.evisu.yy(ip));
            fprintf(fid,'\n');
        end
    end
else
    U=handles.uvisu;
    xo=handles.mvisu.xo;
    yo=handles.mvisu.yo;
    if handles.param.thermo==1
        fprintf(fid,'x [pixel];y [pixel];Ux [pixel];Uy [pixel];T[DL]');
    else
        fprintf(fid,'x [pixel];y [pixel];Ux [pixel];Uy [pixel]');
    end
    if ~handles.sonelt
        fprintf(fid,';dUxdx[];dUxdY[];dUydx[];dUydy[]');
    end
    fprintf(fid,'\n');
    
    for ip=1:length(xo)
        fprintf(fid,'%.3f;%.3f;',handles.param.roi(1)-1+xo(ip),handles.param.roi(3)-1+yo(ip));
        fprintf(fid,'%.3f;%.3f;',U(ip,iim),U(ip+prod(Nnodes),iim));
        if handles.param.thermo==1
            fprintf(fid,'%.3e;',U(ip+2*prod(Nnodes),iim));
        end
        if ~handles.sonelt
            fprintf(fid,'%.3e;%.3e;%.3e;%.3e;',handles.evisu.xx(ip),handles.evisu.xy(ip),handles.evisu.yx(ip),handles.evisu.yy(ip));
        end
        fprintf(fid,'\n');
    end
    
end
fclose(fid);
set(handles.message_text,'String',sprintf('CSV export to %s.....done',filexp))

function vic_circle_set(handles,attractor)
if nargin<2,attractor=0;end
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Click and drag to define a circular zone.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    
    handles.initial_point=point1;
    guidata(handles.figure1,handles);
    xyp=[point1;point1];
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_circle);
    
    waitforbuttonpress
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    if strcmp(get(handles.figure1,'selectiontype'),'normal')
        point2=getposition(handles);
        point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
        xyp=0.5*[(point1+point2),norm(point2-point1)];
        nface=round(2*pi*xyp(2)/5);
        Zcp=xyp(1:2)*[1;1i]+xyp(3)*exp(1i*[(0:2*pi/nface:2*pi)';2*pi]);
        Zcp=[real(Zcp),imag(Zcp)];
        set(gz,'Xdata',Zcp(:,1));
        set(gz,'Ydata',Zcp(:,2));
    else
        delete(gz);
        return;
        set(handles.message_text,'String','Zone definition.....cancelled')
    end
    gz.Tag=num2str(rand(1));
    handles.vic_model.zone{2,end+1}=Zcp;
    handles.vic_model.zone{3,end}=gz;
    handles.vic_model.zone{4,end}=3;
    handles.vic_model.zone{5,end}=xyp;
    handles.vic_model.zone{6,end}=20;
    handles.vic_model.zone{7,end}=2;
    handles.vic_model.zone{8,end}=10;
    handles.vic_model.zone{9,end}=0;
    if ~attractor
        handles.vic_model.zone{1,end}=1;
    else
        handles.vic_model.zone{1,end}=-1;
    end
    
    
    set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
    set(gz,'Color',[1,0.5,0]);
    set(gz,'uicontextmenu',handles.cmenu_vic);
    
    handles=set_param(handles);
    set(handles.message_text,'String','Creating circular zone.....done')
    
else
    set(handles.message_text,'String','Zone definition.....cancelled')
end
set(gcf,'WindowButtonDownFcn',@double_clic);
set(gcf,'WindowButtonMotionFcn',@show_position);
% --------------------------------------------------------------------
function cmenu_param_vic_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_vic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function cmenu_param_vic_line_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_vic_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vic_line_set(handles)

function vic_line_set(handles,attractor)
if nargin<2,attractor=0;end
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
%% pour Rana
set(handles.message_text,'String','Start defining the shape by picking its first point.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    xyp=[point1;point1];
    point2=-10000000;
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(handles.message_text,'String','.....and now pick up points on its path, right click to end.....')
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_poly);
    
    while abs((point2-xyp(1,:))*[1;1i])>0.01*min(handles.sizeim(1:2))
        waitforbuttonpress
        if strcmp(get(handles.figure1,'selectiontype'),'normal')
            point2=getposition(handles);
            point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
            xyp=[xyp;point2;point2];
            set(gz,'Xdata',xyp(:,1));
            set(gz,'Ydata',xyp(:,2));
        else
            if size(xyp,1)<3
                delete(gz);
                set(gcf,'WindowButtonDownFcn',@double_clic);
                set(gcf,'WindowButtonMotionFcn',@show_position);
                return;
                set(handles.message_text,'String','Zone definition.....cancelled')
            else
                break;
            end
        end
    end
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    xyp=xyp(1:2:end,:);
    set(gz,'Xdata',xyp(:,1));
    set(gz,'Ydata',xyp(:,2));
    gz.Tag=num2str(rand(1));
    handles.vic_model.zone{2,end+1}=xyp;
    handles.vic_model.zone{3,end}=gz;
    handles.vic_model.zone{4,end}=4;
    handles.vic_model.zone{6,end}=20;
    handles.vic_model.zone{7,end}=2;
    handles.vic_model.zone{8,end}=10;
    handles.vic_model.zone{9,end}=0;
    if ~attractor
        handles.vic_model.zone{1,end}=1;
    else
        handles.vic_model.zone{1,end}=-1;
    end
    
    set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
    set(gz,'Color',[1,0.5,0]);
    set(gz,'uicontextmenu',handles.cmenu_vic);
    handles=set_param(handles);
    set(handles.message_text,'String','Creating line.....done')
    
else
    set(handles.message_text,'String','Line definition.....cancelled')
end
set(gcf,'WindowButtonDownFcn',@double_clic);
set(gcf,'WindowButtonMotionFcn',@show_position);
%%
% set(handles.message_text,'String','Click both ends of the line.....')
% waitforbuttonpress
% if strcmp(get(handles.figure1,'selectiontype'),'normal')
%     point1=getposition(handles);
%     point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
%     xyp=[point1;point1];
%     point2=-10000000;
%     gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
%     set(gcf,'WindowButtonDownFcn','');
%     set(gcf,'WindowButtonMotionFcn',@follow_mouse_poly);
%
%     waitforbuttonpress
%     if strcmp(get(handles.figure1,'selectiontype'),'normal')
%         point2=getposition(handles);
%         point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
%         xyp=[xyp;point2;point2];
%         set(gz,'Xdata',xyp(:,1));
%         set(gz,'Ydata',xyp(:,2));
%     else
%         delete(gz);
%         set(gcf,'WindowButtonDownFcn',@double_clic);
%         set(gcf,'WindowButtonMotionFcn',@show_position);
%         return;
%         set(handles.message_text,'String','Zone definition.....cancelled')
%     end
%     set(gcf,'WindowButtonDownFcn',@double_clic);
%     set(gcf,'WindowButtonMotionFcn',@show_position);
%     xyp=xyp(1:2:end,:);
%     gz.Tag=num2str(rand(1));
%     handles.vic_model.zone{2,end+1}=xyp;
%     handles.vic_model.zone{3,end}=gz;
%     handles.vic_model.zone{4,end}=4;
%     handles.vic_model.zone{6,end}=20;
%     handles.vic_model.zone{7,end}=2;
%     handles.vic_model.zone{8,end}=10;
%     handles.vic_model.zone{9,end}=0;
%     if ~attractor
%         handles.vic_model.zone{1,end}=1;
%     else
%         handles.vic_model.zone{1,end}=-1;
%     end
%
%     set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
%     set(gz,'Color',[1,0.5,0]);
%     set(gz,'uicontextmenu',handles.cmenu_vic);
%     handles=set_param(handles);
%     set(handles.message_text,'String','Creating line.....done')
%
% else
%     set(handles.message_text,'String','Line definition.....cancelled')
% end
% set(gcf,'WindowButtonDownFcn',@double_clic);
% set(gcf,'WindowButtonMotionFcn',@show_position);


% --------------------------------------------------------------------
function cmenu_param_vic_circle_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_vic_circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vic_circle_set(handles)


% --------------------------------------------------------------------
function cmenu_vic_nbelt_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_nbelt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_vic_degree_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_vic_type_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_vic_type_solid_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_type_solid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.vic_model.zone(3,:)),gco);
handles.vic_model.zone{1,iz}=-1;
guidata(handles.figure1,handles)


% --------------------------------------------------------------------
function cmenu_vic_type_wire_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_type_wire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.vic_model.zone(3,:)),gco);
handles.vic_model.zone{1,iz}=1;
guidata(handles.figure1,handles)


% --------------------------------------------------------------------
function cmenu_vic_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.vic_model.zone(3,:)),gco);
set_vic_param(handles,iz);

% --------------------------------------------------------------------
function cmenu_vic_deg_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_deg_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='vic_deg';
iz=findg((handles.vic_model.zone(3,:)),gco);
handles.active_zone=iz;
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust the degree of the B-spline.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_vic_nbelt_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_nbelt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='vic_nbelt';
iz=findg((handles.vic_model.zone(3,:)),gco);
handles.active_zone=iz;
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust the number of elements.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_vic_delete_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iz=findg((handles.vic_model.zone(3,:)),gco);
handles.vic_model.zone(:,iz)=[];
delete(gco)
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_vic_ep_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_ep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_vic_ep_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_ep_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='vic_ep';
iz=findg((handles.vic_model.zone(3,:)),gco);
handles.active_zone=iz;
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust the transition length.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_field_vref_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_vref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.vref=~handles.vref;
if handles.showplot==1
    for id=1:length(handles.fgage)
        try delete(handles.fgage(id));catch, end
        try delete(100+handles.fgage(id));catch, end
        try delete(200+handles.fgage(id));catch, end
        try delete(300+handles.fgage(id));catch, end
        try delete(400+handles.fgage(id));catch, end
        try delete(4000+handles.fgage(id));catch, end
        try delete(500+handles.fgage(id));catch, end
        try delete(600+handles.fgage(id));catch, end
    end
    cmenu_export_vic_plot_Callback(handles.figure1, [], handles)
    handles=guidata(handles.figure1);
end
handles=plot_mesh(handles);
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_vic_t_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_vic_t_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_vic_t_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.scroll_adjusted='vic_t';
iz=findg((handles.vic_model.zone(3,:)),gco);
handles.active_zone=iz;
guidata(hObject,handles);
set(handles.message_text,'String','Scroll mouse wheel to adjust the wire thickness.....')
set(handles.figure1,'WindowButtonDownFcn',@stop_adjust);
set(handles.figure1,'WindowScrollWheelFcn',@scroll_adjust);


% --------------------------------------------------------------------
function cmenu_param_vic_point_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_param_vic_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vic_zone_set(handles)


function vic_zone_set(handles,attractor)
if nargin<2,attractor=0;end
set(0,'CurrentFigure',handles.figure1);
handles=reference_frame(handles);
set(handles.message_text,'String','Click opposite corners of a rectangular zone.....')
waitforbuttonpress
if strcmp(get(handles.figure1,'selectiontype'),'normal')
    point1=getposition(handles);
    point1=max([1,1],min(handles.sizeim(1:2),point1(1,1:2)));
    
    handles.initial_point=point1;
    guidata(handles.figure1,handles);
    xyp=[point1;point1];
    gz=plot(xyp(:,1),xyp(:,2),'Color',[0,0.5,0],'LineStyle','--','LineWidth',2);
    set(gcf,'WindowButtonDownFcn','');
    set(gcf,'WindowButtonMotionFcn',@follow_mouse_rect);
    
    waitforbuttonpress
    set(gcf,'WindowButtonDownFcn',@double_clic);
    set(gcf,'WindowButtonMotionFcn',@show_position);
    
    point2=getposition(handles);
    point2=max([1,1],min(handles.sizeim(1:2),point2(1,1:2)));
    ly=abs((point1(1)-point2(1)));
    lx=abs((point1(2)-point2(2)));
    if strcmp(get(handles.figure1,'selectiontype'),'normal')&&(~(max(lx,ly)<1))
        roi=[max(1,min(handles.sizeim(1),sort(round([point1(1),point2(1)])))),max(1,min(handles.sizeim(2),sort(round([point1(2),point2(2)]))))];
        xyp=[roi([1,2,2,1,1])',roi([3,3,4,4,3])'];
        gz=plot(xyp(:,1),xyp(:,2),'Color',[1,0.5,0.],'LineStyle','-','LineWidth',2);
        gz.Tag=num2str(rand(1));
        handles.vic_model.zone{2,end+1}=xyp;
        handles.vic_model.zone{3,end}=gz;
        handles.vic_model.zone{4,end}=1;
        handles.vic_model.zone{6,end}=20;
        handles.vic_model.zone{7,end}=2;
        handles.vic_model.zone{8,end}=10;
        handles.vic_model.zone{9,end}=0;
        if ~attractor
            handles.vic_model.zone{1,end}=1;
        else
            handles.vic_model.zone{1,end}=-1;
        end
        
        set(gz,'LineStyle','-','ButtonDownFcn',@zone_adjust)
        set(gz,'Color',[1,0.5,0]);
        set(gz,'uicontextmenu',handles.cmenu_vic);
        handles=set_param(handles);
        set(handles.message_text,'String','Creating line.....done')
        
    else
        delete(gz);
        set(handles.message_text,'String','Zone definition.....cancelled')
    end
else
    set(handles.message_text,'String','Zone definition.....cancelled')
end


% --------------------------------------------------------------------
function cmenu_export_vic_plot_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_vic_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ids=findg((handles.vic_model.zone(3,:)),gco);
if isempty(ids)
    ids=1:size(handles.vic_model.zone,2);
end
for ii=1:length(ids)
    id=ids(ii);
    handles.fgage(id)=id;
    guidata(handles.figure1,handles);
    pix2m=handles.param.pixel_size;
    load(handles.param.result_file,'-mat','s','curv','fleche')
    
    fleche=fleche{id};
    curv=curv{id};
    s=s{id};
    dt=ceil(size(curv,2)/7);
    fgage=figure(id);
    delete(gca)
    axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    set(fgage,'Name',sprintf('Zone #%d: curvature',id));
    set(fgage,'NumberTitle','off');
    box('on');
    hold('all');
    for iim=1:dt:size(curv,2)
        plot(s'*pix2m,curv(:,iim)/pix2m,'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
    end
    if pix2m==1
        ylabel('Curvature [1/pixel]','FontSize',20,'FontName','Times');
        xlabel('s [pixel]','FontSize',20,'FontName','Times');
    else
        ylabel('Curvature [1/m]','FontSize',20,'FontName','Times');
        xlabel('s [m]','FontSize',20,'FontName','Times');
    end
    
    leg=legend(axes1,'show','Location','NorthWest');
    set(leg,'Box','off')
    
    fgage=figure(id+100);
    delete(gca)
    axes1 = axes('Parent',fgage,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    set(fgage,'Name',sprintf('Zone #%d: displacement',id));
    set(fgage,'NumberTitle','off');
    box('on');
    hold('all');
    for iim=1:dt:size(curv,2)
        if handles.vref==0
            plot(s'*pix2m,(fleche(:,iim)-fleche(:,1))*pix2m,'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
        else
            plot(s'*pix2m,fleche(:,iim)*pix2m,'DisplayName',num2str((iim)),'Parent',axes1,'LineWidth',2);
        end
        
    end
    if pix2m==1
        ylabel('Displacement [pixel]','FontSize',20,'FontName','Times');
        xlabel('s [pixel]','FontSize',20,'FontName','Times');
    else
        ylabel('Displacement [1/m]','FontSize',20,'FontName','Times');
        xlabel('s [m]','FontSize',20,'FontName','Times');
    end
    
    leg=legend(axes1,'show','Location','NorthWest');
    set(leg,'Box','off')
end


% --------------------------------------------------------------------
function cmenu_export_vic_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_export_vic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function new_3d_analysis_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to new_3d_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
U3D


% --------------------------------------------------------------------
function cmenu_beam_kine_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_kine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_beam_kine_euler_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_kine_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.beam_type='euler';
handles=set_param(handles);


% --------------------------------------------------------------------
function cmenu_beam_kine_timosh_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_beam_kine_timosh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.beam_model.beam_type='timoshenko';
handles=set_param(handles);

% --------------------------------------------------------------------
function cmenu_field_data_shearstrain_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_field_data_shearstrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field=2;
handles.ecomp=8;
handles=plot_mesh(handles);
handles=set_param(handles);


