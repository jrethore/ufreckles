function varargout = USolver(varargin)
% USOLVER MATLAB code for USolver.fig
%      USOLVER, by itself, creates a new USOLVER or raises the existing
%      singleton*.
%
%      H = USOLVER returns the handle to a new USOLVER or the handle to
%      the existing singleton*.
%
%      USOLVER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in USOLVER.M with the given input arguments.
%
%      USOLVER('Property','Value',...) creates a new USOLVER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before USolver_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to USolver_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help USolver

% Last Modified by GUIDE v2.5 19-Jan-2015 13:34:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @USolver_OpeningFcn, ...
                   'gui_OutputFcn',  @USolver_OutputFcn, ...
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


% --- Executes just before USolver is made visible.
function USolver_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to USolver (see VARARGIN)

% Choose default command line output for USolver
handles.output = hObject;
handles.umaster=varargin{1};
%handles.data=guidata(varargin{1});

% Update handles structure
guidata(hObject, handles);




% UIWAIT makes USolver wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = USolver_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function pix2m_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pix2m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pix2m_edit as text
%        str2double(get(hObject,'String')) returns contents of pix2m_edit as a double
Umaster('update_param',handles.umaster);

% --- Executes during object creation, after setting all properties.
function pix2m_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix2m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nscale_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nscale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nscale_edit as text
%        str2double(get(hObject,'String')) returns contents of nscale_edit as a double
Umaster('update_param',handles.umaster);


% --- Executes during object creation, after setting all properties.
function nscale_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nscale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function convergence_edit_Callback(hObject, eventdata, handles)
% hObject    handle to convergence_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of convergence_edit as text
%        str2double(get(hObject,'String')) returns contents of convergence_edit as a double
Umaster('update_param',handles.umaster);

% --- Executes during object creation, after setting all properties.
function convergence_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to convergence_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxiter_edit_Callback(hObject, eventdata, handles)
% hObject    handle to maxiter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxiter_edit as text
%        str2double(get(hObject,'String')) returns contents of maxiter_edit as a double
Umaster('update_param',handles.umaster);


% --- Executes during object creation, after setting all properties.
function maxiter_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxiter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in restart_edit.
function restart_edit_Callback(hObject, eventdata, handles)
% hObject    handle to restart_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restart_edit
Umaster('update_param',handles.umaster);

% --- Executes on selection change in do_pgd_edit.
function do_pgd_edit_Callback(hObject, eventdata, handles)
% hObject    handle to do_pgd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns do_pgd_edit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from do_pgd_edit
Umaster('update_param',handles.umaster);


% --- Executes during object creation, after setting all properties.
function do_pgd_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to do_pgd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normalize_grey_levels_edit.
function normalize_grey_levels_edit_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_grey_levels_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_grey_levels_edit
Umaster('update_param',handles.umaster);



function psample_edit_Callback(hObject, eventdata, handles)
% hObject    handle to psample_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of psample_edit as text
%        str2double(get(hObject,'String')) returns contents of psample_edit as a double
Umaster('update_param',handles.umaster);


% --- Executes during object creation, after setting all properties.
function psample_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psample_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
