function varargout = Uaba(varargin)
% UABA MATLAB code for Uaba.fig
%      UABA, by itself, creates a new UABA or raises the existing
%      singleton*.
%
%      H = UABA returns the handle to a new UABA or the handle to
%      the existing singleton*.
%
%      UABA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UABA.M with the given input arguments.
%
%      UABA('Property','Value',...) creates a new UABA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Uaba_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Uaba_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Uaba

% Last Modified by GUIDE v2.5 27-Feb-2014 19:59:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Uaba_OpeningFcn, ...
                   'gui_OutputFcn',  @Uaba_OutputFcn, ...
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


% --- Executes just before Uaba is made visible.
function Uaba_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Uaba (see VARARGIN)

% Choose default command line output for Uaba
handles.output = hObject;

% Update handles structure
handles.model=varargin{2};
handles.param=varargin{1};
  set(handles.nlayers_edit,'Enable','off');
  set(handles.nlayers_edit,'Enable','on');
  set(handles.thickness_edit,'Enable','off');
  set(handles.thickness_edit,'Enable','on');


if isfield(handles.param,'calibration_parameters')
    dflag=1;set(handles.thickness_text,'String','Thickness in m');
elseif(length(handles.param.roi)==6)
    dflag=1;set(handles.thickness_text,'String','Thickness in voxel');
else
    dflag=0;set(handles.thickness_text,'String','Thickness in pixel');
end
handles.dflag=dflag;
handles.model.rmesh=0;
if dflag
    set(handles.elt_type,'String',{'Element type';'3D solid';'3D shell';'3D membrane'});
handles.doextrusion=1;
    set(handles.extrusion_method,'String',{'Extrusion method','3D constant thickness','3D using Topo','Shell'})
    set(handles.extrusion_method,'Enable','on');
set(handles.extrusion_params,'Visible','on');
else
     set(handles.extrusion_method,'String',{'3D',})
    set(handles.extrusion_method,'Enable','off');
   handles.model.extrusion_parameters.type='3d';
   handles.model.extrusion_parameters.thickness=10;
     set(handles.elt_type,'String',{'Element type';'2D plane strain';'2D plane stress';'3D solid';'3D shell';'3D membrane'});
handles.doextrusion=0;
set(handles.extrusion_params,'Visible','off');
end
handles.model.cpflag=false;
handles.model.shell=false;
handles.model.menbrane=false;

        rint=true;
        if isfield(handles.model,'reduced_integration')
            rint=handles.model.reduced_integration;
        else
            handles.model.reduced_integration=rint;
        end
        set(handles.rint,'Value',double(rint));
        nlflag=false;
        if isfield(handles.model,'nlgeom')
            nlflag=handles.model.nlgeom;
        else
            handles.model.nlgeom=nlflag;
        end
        set(handles.nlgeom,'Value',double(nlflag));
 guidata(hObject, handles);

% UIWAIT makes Uaba wait for user response (see UIRESUME)
uiwait(handles.figure1);
%handles=guidata(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Uaba_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure

if~handles.doextrusion
   handles.model=rmfield(handles.model,'extrusion_parameters');
else
if ~handles.dflag&&strcmp(handles.model.extrusion_parameters.type,'3d')
handles.model.extrusion_parameters.thickness=0.5*handles.model.extrusion_parameters.thickness;
handles.model.extrusion_parameters.nlayers=ceil(0.5*handles.model.extrusion_parameters.nlayers);
end
end
varargout{1}=handles.model;
varargout{2}=handles.doextrusion;

delete(hObject);


% --- Executes on selection change in elt_type.
function elt_type_Callback(hObject, eventdata, handles)
% hObject    handle to elt_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns elt_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elt_type
handles.model.cpflag=false;
handles.model.shell=false;
handles.model.menbrane=false;
handles.doextrusion=0;
set(handles.nlayers_edit,'Enable','on');
if handles.dflag
    switch get(hObject,'Value')
        case 3
handles.model.shell=true;
handles.model.extrusion_parameters.nlayers=1;
set(handles.nlayers_edit,'String','1');
set(handles.nlayers_edit,'Enable','off');
nlayers_edit_Callback(handles.nlayers_edit,[],handles);
        case 4
handles.model.menbrane=false;
handles.model.extrusion_parameters.nlayers=1;
set(handles.nlayers_edit,'String','1');
set(handles.nlayers_edit,'Enable','off');
nlayers_edit_Callback(handles.nlayers_edit,[],handles);
    end
            handles.doextrusion=1;
else
    
     switch get(hObject,'Value')
         case 2
handles.model.cpflag=true;
             
        case 5
handles.model.shell=true;
handles.model.extrusion_parameters.nlayers=1;
set(handles.nlayers_edit,'String','1');
set(handles.nlayers_edit,'Enable','off');
nlayers_edit_Callback(handles.nlayers_edit,[],handles);

        case 6
handles.model.menbrane=false;
set(handles.nlayers_edit,'String','1');
set(handles.nlayers_edit,'Enable','off');
nlayers_edit_Callback(handles.nlayers_edit,[],handles);
    end
   if get(hObject,'Value')>3, handles.doextrusion=1;end
end
if handles.doextrusion
    set(handles.extrusion_params,'Visible','on');
else
    set(handles.extrusion_params,'Visible','off');
end
guidata(hObject, handles);
 



% --- Executes during object creation, after setting all properties.
function elt_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elt_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in rmesh.
function rmesh_Callback(hObject, eventdata, handles)
% hObject    handle to rmesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns rmesh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rmesh
handles.model.rmesh=get(hObject,'Value')-1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function rmesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rmesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nlgeom.
function nlgeom_Callback(hObject, eventdata, handles)
% hObject    handle to nlgeom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nlgeom
handles.model.nlgeom=get(hObject,'Value')>0; 
guidata(hObject, handles);

% --- Executes on button press in rint.
function rint_Callback(hObject, eventdata, handles)
% hObject    handle to rint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rint
handles.model.reduced_integration=get(hObject,'Value')>0; 
guidata(hObject, handles);


% --- Executes on selection change in extrusion_method.
function extrusion_method_Callback(hObject, eventdata, handles)
% hObject    handle to extrusion_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns extrusion_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from extrusion_method
set(handles.nlayers_edit,'Enable','on');
switch get(hObject,'Value')
    case 1
        handles.model.extrusion_parameters='';
    case 2
        handles.model.extrusion_parameters='solid';
    case 3
        handles.model.extrusion_parameters='solid-from-topo';
    case 4
        handles.model.extrusion_parameters='shell';
  set(handles.nlayers_edit,'Enable','off');
  set(handles.nlayers_edit,'String','1');
      
end        
  guidata(hObject, handles);
      

% --- Executes during object creation, after setting all properties.
function extrusion_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extrusion_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thickness_edit_Callback(hObject, eventdata, handles)
% hObject    handle to thickness_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thickness_edit as text
%        str2double(get(hObject,'String')) returns contents of thickness_edit as a double
handles.model.extrusion_parameters.thickness=eval(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function thickness_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thickness_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nlayers_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nlayers_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nlayers_edit as text
%        str2double(get(hObject,'String')) returns contents of nlayers_edit as a double
handles.model.extrusion_parameters.nlayers=round(eval(get(hObject,'String')))+1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nlayers_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nlayers_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1)


% --- Executes on button press in finish.
function finish_Callback(hObject, eventdata, handles)
% hObject    handle to finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)
