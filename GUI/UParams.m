function varargout = UParams(varargin)
% UPARAMS MATLAB code for UParams.fig
%      UPARAMS, by itself, creates a new UPARAMS or raises the existing
%      singleton*.
%
%      H = UPARAMS returns the handle to a new UPARAMS or the handle to
%      the existing singleton*.
%
%      UPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UPARAMS.M with the given input arguments.
%
%      UPARAMS('Property','Value',...) creates a new UPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UParams_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UParams

% Last Modified by GUIDE v2.5 01-Mar-2014 22:17:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UParams_OpeningFcn, ...
    'gui_OutputFcn',  @UParams_OutputFcn, ...
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


% --- Executes just before UParams is made visible.
function UParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UParams (see VARARGIN)

% Choose default command line output for UParams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(0,'CurrentFigure',handles.figure1)


param=varargin{1};
model=varargin{2};

pbuff{1}='Analysis type';
switch param.analysis
    case 'correlation'
        vbuff{1}=sprintf('DIC ');
    case 'mechanics'
        vbuff{1}=sprintf('FEM ');
end

pbuff{end+1}=sprintf('Algorithm');
if param.restart
    vbuff{end+1}=sprintf('independent ');
else
    if isfield(param,'time_step')
        vbuff{end+1}=sprintf('space-time');
        pbuff{end+1}='';
        vbuff{end+1}=sprintf('dt=%.1f',param.time_step);
    else
        pbuff{end+1}='Prediction';
        if param.do_pgd_prediction
            vbuff{end+1}='sequential ';
            vbuff{end+1}='automatic ';
        else
            vbuff{end+1}='sequential ';
            vbuff{end+1}='explicit ';
        end
    end
end
pbuff{end+1}='Physical pixel size';
vbuff{end+1}=sprintf('%.2e ',param.pixel_size);
pbuff{end+1}='Maximum iteration';
vbuff{end+1}=sprintf('%d ',param.iter_max);
pbuff{end+1}='Convergence criterion';
vbuff{end+1}=sprintf('%.2e ',param.convergance_limit);
pbuff{end+1}='GL Normalization';
if param.normalize_grey_level
    vbuff{end+1}='element-wise';
else
    vbuff{end+1}='global ';
end
pbuff{end+1}='Coarse graining level';
vbuff{end+1}=sprintf('%d ',model.nscale);
pbuff{end+1}='Pixel skip';
vbuff{end+1}=sprintf('%d ',param.psample);

pbuff{end+1}='';
vbuff{end+1}='';


pbuf{1}='Basis function';
switch model.basis
    case 'fem'
        if isfield(model,'degree')
        vbuf{1}=sprintf('NURBS');
                pbuf{end+1}='Mesh size';
                vbuf{end+1}=sprintf('[%.1f,%.1f] ',model.mesh_size);
                pbuf{end+1}='Degree';
                vbuf{end+1}=sprintf('[%d,%d] ',model.degree);
        else
        vbuf{1}=sprintf('FEM ');
        switch model.mesh_type
            case 1
                pbuf{end+1}='Mesh type';
                vbuf{end+1}=sprintf('structured ');
                pbuf{end+1}='Element type';
                vbuf{end+1}=sprintf('%d ',model.element_type);
                pbuf{end+1}='Mesh size';
                vbuf{end+1}=sprintf('[%.1f,%.1f] ',model.mesh_size);
            case 2
                pbuf{end+1}='Mesh type';
                vbuf{end+1}=sprintf('unstructured ');
                pbuf{end+1}='Element type';
                vbuf{end+1}=sprintf('%d ',3);
                pbuf{end+1}='Mesh size';
                vbuf{end+1}=sprintf('[%.1f,%.1f] ',model.mesh_size);
            case 3
                pbuf{end+1}='Mesh type';
                vbuf{end+1}=sprintf('external ');
                pbuf{end+1}='Mesh file';
                vbuf{end+1}=sprintf('%s ',model.mesh_file);
                
        end
        end
        if param.regularization_parameter
            pbuf{end+1}='Smoothing';
            switch param.regularization_type
                case 'median'
                    vbuf{end+1}=sprintf('median ');
            pbuf{end+1}='Number of neighboor';
            vbuf{end+1}=sprintf('%d ',param.regularization_parameter);
                case 'tiko'
                    vbuf{end+1}=sprintf('strain ');
            pbuf{end+1}='Cut-off wave length';
            vbuf{end+1}=sprintf('%d ',param.regularization_parameter);
                case 'equilibrium_gap'
                    vbuf{end+1}=sprintf('stress ');
            pbuf{end+1}='Cut-off wave length';
            vbuf{end+1}=sprintf('%d ',param.regularization_parameter);
            end
        else
            pbuf{end+1}='Smoothing';
            vbuf{end+1}=sprintf('off ');
        end
        
    case 'uni'
        vbuf{1}=sprintf('STRAIN-GAGE ');
    case 'beam'
        vbuf{1}=sprintf('BEAM ');
        pbuf{end+1}='Beam type';
        vbuf{end+1}=sprintf('4PB ');
        pbuf{end+1}='Additional axial strain';
        if model.exx
            vbuf{end+1}=sprintf('on ');
        else
            vbuf{end+1}=sprintf('off ');
        end
    case 'beam-nurbs'
        vbuf{1}=sprintf('BEAM ');
        pbuf{end+1}='Beam type';
        vbuf{end+1}=sprintf('NURBS ');
        pbuf{end+1}='Degree';
        vbuf{end+1}=sprintf('%d ',model.degree);
        pbuf{end+1}='Number of elements';
        vbuf{end+1}=sprintf('%d ',model.nb_element);
        pbuf{end+1}='Beam type';
        switch model.beam_type
            case 'euler'
        vbuf{end+1}=sprintf('Euler-Bernoulli ');
            case 'timoshenko'
         vbuf{end+1}=sprintf('Timoshenko ');
       end
        
        pbuf{end+1}='Additional axial strain';
        if model.exx
            vbuf{end+1}=sprintf('on ');
        else
            vbuf{end+1}=sprintf('off ');
        end
        
        
end

set(handles.params_text,'String',pbuff);
set(handles.params_values_text,'String',vbuff);

set(handles.model_text,'String',pbuf);
set(handles.model_values_text,'String',vbuf);





% UIWAIT makes UParams wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UParams_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
