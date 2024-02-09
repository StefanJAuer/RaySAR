function varargout = RaySAR(varargin)
% RAYSAR M-file for RaySAR.fig
%      RAYSAR, by itself, creates a new RAYSAR or raises the existing
%      singleton*.
%
%      H = RAYSAR returns the handle to a new RAYSAR or the handle to
%      the existing singleton*.
%
%      RAYSAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAYSAR.M with the given input arguments.
%
%      RAYSAR('Property','Value',...) creates a new RAYSAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RaySAR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RaySAR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RaySAR

% Last Modified by GUIDE v2.5 20-Sep-2011 16:58:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RaySAR_OpeningFcn, ...
                   'gui_OutputFcn',  @RaySAR_OutputFcn, ...
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


% --- Executes just before RaySAR is made visible.
function RaySAR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RaySAR (see VARARGIN)

% Choose default command line output for RaySAR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RaySAR wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global startdir;

startdir=pwd;


% --- Outputs from this function are returned to the command line.
function varargout = RaySAR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Simulation_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Simulation_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Simulation_Data;

% --------------------------------------------------------------------
function Maps_Callback(hObject, eventdata, handles)
% hObject    handle to Maps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Maps;

% --------------------------------------------------------------------
function Elevation_Callback(hObject, eventdata, handles)
% hObject    handle to Elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Elevation;

% --------------------------------------------------------------------
function Focused_Pts_Callback(hObject, eventdata, handles)
% hObject    handle to Focused_Pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Focused_Pts;

% --------------------------------------------------------------------
function Intersection_Callback(hObject, eventdata, handles)
% hObject    handle to Intersection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Intersection;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Buttonname = questdlg('Exit RaySAR?', 'Exit','Yes');

TF = strcmp(Buttonname,'Yes');

if TF == 1
    
    figure_tags = findall(0,'type','figure');
    delete(figure_tags);
    delete(hObject);
    clear global;
end


% --- Executes on button press in Link_Model.
function Link_Model_Callback(hObject, eventdata, handles)
% hObject    handle to Link_Model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Model_path;

Model_path = uigetdir;

if Model_path == 0
    
    set(handles.path_model,'string','');
    
else
    
    set(handles.path_model,'string',Model_path);
    
end


% --- Executes on button press in Link_Output.
function Link_Output_Callback(hObject, eventdata, handles)
% hObject    handle to Link_Output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Output_path Model_path;
  
Output_path = uigetdir(Model_path);

if Output_path == 0
    
    set(handles.path_output,'string','');
    
else
    
    set(handles.path_output,'string',Output_path);
    
end


% --- Executes on button press in Link_SAR_Image.
function Link_SAR_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Link_SAR_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SAR_image_path Model_path;

SAR_image_path = uigetdir(Model_path);
    
if SAR_image_path == 0
    
    set(handles.path_SAR_image,'string','');
    
else
    
    set(handles.path_SAR_image,'string',SAR_image_path);
    
end
