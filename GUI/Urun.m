function varargout = Urun(varargin)
% URUN MATLAB code for Urun.fig
%      URUN, by itself, creates a new URUN or raises the existing
%      singleton*.
%
%      H = URUN returns the handle to a new URUN or the handle to
%      the existing singleton*.
%
%      URUN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in URUN.M with the given input arguments.
%
%      URUN('Property','Value',...) creates a new URUN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Urun_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Urun_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Urun

% Last Modified by GUIDE v2.5 10-Jan-2014 14:44:15

% Begin initialization code - DO NOT EDIT
warning('off','all');
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Urun_OpeningFcn, ...
    'gui_OutputFcn',  @Urun_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Urun is made visible.
function Urun_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Urun (see VARARGIN)

% Choose default command line output for Urun
handles.output = hObject;
handles.jobs={};
handles.run=0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Urun wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Urun_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function open_job_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to open_job (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fils,paths]=uigetfile({'*.dat','Ufreckles input files';'*.ufr','Ufreckles input ascii files'},'Open data sets');
%if ~iscell(fils)
%   fils={fils};
%end
handles.jobs{1}=fils;
handles.jobs{2}=paths;
%  nbjobs=size(handles.jobs,1);
%  for id=1:length(fils)
%      handles.jobs{nbjobs+id,1}=fils{id};
%      handles.jobs{nbjobs+id,2}=paths;
%  end
guidata(hObject, handles);
set(handles.job_list,'String',handles.jobs{1});
set(handles.runing_text,'String',sprintf('No job running...'))


% --- Executes on selection change in job_list.

% --- Executes during object creation, after setting all properties.
function job_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to job_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function run_live(hrun,param,model)
global phix phiy
set(hrun,'String','Stop');
run=1;
nmod=0;
clear functions
param.onflight=1;
nscale=1;
model.nscale=1;
LoadParameters(param);
LoadParameters(model,nmod);
ReferenceImage(nmod);
LoadMask(nmod);
filres=param.result_file;
iscale=nscale;
switch model.basis
    case 'fem'
LoadMeshes(nmod);
LoadMat(nmod);
load(fullfile('TMP','0_mesh_0'),'Nnodes','Nelems','xo','yo','conn','elt','ng','rflag','rint');
tips=zeros(size(model.zone,2),2);
cracks=[];
Ks=[];
U=zeros(2*prod(Nnodes),1);
param.ulive=1;
save(filres,'U','Ks','tips','cracks','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
CreateBasisFunction(iscale,nmod);
ComputeGradFPhi(iscale,nmod);
CreateGradBasisFunction(iscale,nmod);
AssembleCorrelationOperator(iscale,nmod);
    case 'uni'
        gages=model.zone(2,:);
        assert(length(gages)==1)
param.ulive=1;ig=1;
    gage=gages{ig};
    roi=param.roi;
    sizeim=[roi(2)-roi(1),roi(4)-roi(3)]+1;
    xp=gage(1:4,1)+roi(1)-1;
    yp=gage(1:4,2)+roi(3)-1;
    xo=round(xp);
    yo=round(yp);
elt=4;
zo=1;
conn=1:4;
Nnodes=[length(xo),1,1];
Nelems=[length(elt),1,1];
        Smesh=[max(xo)-min(xo),max(yo)-min(yo)];
        Vmesh=[max(xo)-min(xo),max(yo)-min(yo),max(zo)-min(zo)];
ns=ones(1,3);
    rflag=1;
    rint=false;
    ng=0;
    selected=[];
      U1=zeros(6,1);
      U=zeros(2*prod(Nnodes),1);
      Up=cell(length(gages),1);
      for ig=1:length(gages)
          Up{ig}=zeros(9,1);
      end
    [Yi, Xi]=meshgrid(1:(sizeim(2)),1:(sizeim(1)));
    mask=inpolygon(Xi,Yi,xp([1:length(xp),1])-roi(1)+1,yp([1:length(xp),1])-roi(3)+1);
              save(filres,'U','Up','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
        mask=diag(sparse(mask(:)));
        save(fullfile('TMP',sprintf('%d_mask_%d',nmod,iscale-1)),'mask','-append');
CreateBasisFunction(iscale,nmod);
ComputeGradFPhi(iscale,nmod);
AssembleCorrelationOperator(iscale,nmod);
                    save(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'rflag','rint','xo','yo','zo','elt','conn','Nnodes','Nelems','Smesh','Vmesh','selected','ng','ns');

        try
            delete([strrep(filres,'.res',''),'-gage.csv']);
        end
        dlmwrite([strrep(filres,'.res',''),'-gage.csv'],zeros(1,4));
       
end
        switch model.basis
            case 'uni'
               indp=xo-roi(1)+1+(yo-roi(3)+1-1)*sizeim(1);
    phixo=phix(indp,:);
    phiyo=phiy(indp,:);
     
                
        end
[~,fil,ext]=fileparts(param.reference_image);
hmaster=findobj('Tag',param.umaster);
for ih=1:numel(hmaster)
    if ~isempty(strfind(hmaster(ih).Name,'UFreckles'))
        hmaster=hmaster(ih);
        break;
    end
end
handles=guidata(hmaster);
Umaster('reset_param',handles);
Umaster('loading_mat_file',hmaster,handles,param.result_file,[],1);
nim=1;
imdefs={param.deformed_image,param.deformed_image};
load(fullfile('TMP','params'),'param');
param.deformed_image=imdefs;
while run
    
    pause(0.1);
    listing=dir(sprintf('*%s',ext'));
    go=size(listing,1)>(nim);
    if go
        
        nim=nim+1;
        listing=listing(end);
        param.deformed_image{nim}=listing.name;
        save(fullfile('TMP','params'),'param');
        switch model.basis
            case 'fem'
        U=[U,zeros(2*prod(Nnodes),1)];
        [U]=SolvePreview(U,nscale,nmod,nim);
            case 'uni'
        U1=[U1,zeros(6,1)];
        [U1]=SolvePreview(U1,nscale,nmod,nim);
            Ux=phixo*U1;
            Uy=phiyo*U1;  
            exx=U1(4,:)/mean(sizeim);
                eyy=U1(5,:)/mean(sizeim);
                exy=U1(6,:)/mean(sizeim);

U=[Ux;Uy];
Up{1}=[U1;exx;eyy;exy];
        
        dlmwrite([strrep(filres,'.res',''),'-gage.csv'],[nim,exx(end),eyy(end),exy(end)],'-append');
        
        end
        handles=guidata(hmaster);
        handles.uvisu=U;
        handles.param.deformed_image=param.deformed_image;
        handles.animation.nbstep=size(U,2);
        handles.animation.frames=1:size(U,2);
        handles.animation.iim=handles.animation.nbstep;
        handles=Umaster('display_frame',handles);
        guidata(hmaster,handles);
        save(filres,'U','param','-append');
        switch model.basis
            case 'uni'
        save(filres,'Up','-append');  
        handles.uni_model.Up=Up;
        guidata(hmaster,handles);
        Umaster('plot_gage_data',handles,handles.uni_model.zone{2,1},1);
        end
    end
    rhandles=guidata(hrun);
    run=rhandles.run;
end
delete(fullfile('VTK',['camr*','-error.vtk']));
set(hrun,'String','Run');





% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%for id=1:size(handles.jobs,1)
if handles.run==0
    handles.run=1;
    guidata(handles.figure1,handles);
    filname=handles.jobs{1};
    path=handles.jobs{2};
    cd(path);
    [pp,jobname,ext]=fileparts(filname);
    set(handles.runing_text,'String',sprintf('Running %s...',jobname))
    pause(0.1);
    if strcmp(ext,'.ufr')
        [param,model]=readINPFile(filname);
    else
        load(filname,'-mat','param','model');
    end
    if ~isfield(param,'da'), param.da=0;end
    if ~isfield(param,'ulive'), param.ulive=0;end
    if ~isfield(param,'thermo'), param.thermo=0;end
    if ~isfield(param,'psample'), param.psample=1;end
    if param.thermo==1
        if iscell(param.ir_calibration_data)
            pini=cd;
            
            T=UTcalib(param.ir_calibration_data{1},param.ir_calibration_data{2});
            cd(pini);
            param.ir_calibration_data=T;
            save(filname,'param','-append')
        end
        
    end
    switch model.basis
        case 'fem'
            if ~isfield(model,'phantom_nodes'), model.phantom_nodes=0;end
            switch param.analysis
                case 'correlation'
                    if param.ulive==1
                        run_live(handles.run_button,param,model);
                    else
                        if isfield(param,'detect')
                            if param.detect
                                run_crack_detect_job(param,model);
                            else
                                if isfield(model,'degree')
                                    run_nurbs_job(param,model);
                                else
                                    run_fem_job(param,model);
                                    
                                end
                            end
                        else
                            if isfield(model,'degree')
                                run_nurbs_job(param,model);
                            else
                                run_fem_job(param,model);
                            end
                        end
                        handles.run=0;
                        guidata(handles.figure1,handles);
                    end
                    
                case 'mechanics'
                    if param.da>0
                        run_crack_propa_job(param,model);
                    else
                        if model.phantom_nodes==0
                            run_fea_job(param,model);
                        else
                            clear run_fdfea_job
                            run_fdfea_job(param,model);
                        end
                    end
                    handles.run=0;
                    guidata(handles.figure1,handles);
            end
        case 'uni'
                    if param.ulive==1
                        run_live(handles.run_button,param,model);
                    else
            run_uni_job(param,model);
            handles.run=0;
            guidata(handles.figure1,handles);
                    end
        case 'vic'
            run_vic_job(param,model);
            handles.run=0;
            guidata(handles.figure1,handles);
        case {'beam','beam-nurbs'}
            run_beam_job(param,model);
            handles.run=0;
            guidata(handles.figure1,handles);
            
    end
    set(handles.runing_text,'String',sprintf('End running %s !',jobname))
else
    handles.run=0;
    guidata(handles.figure1,handles);
end
%set(handles.job_list,'String','');


%end
