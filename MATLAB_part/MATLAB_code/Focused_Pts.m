function varargout = Focused_Pts(varargin)
% FOCUSED_PTS M-file for Focused_Pts.fig
%      FOCUSED_PTS, by itself, creates a new FOCUSED_PTS or raises the existing
%      singleton*.
%
%      H = FOCUSED_PTS returns the handle to a new FOCUSED_PTS or the handle to
%      the existing singleton*.
%
%      FOCUSED_PTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOCUSED_PTS.M with the given input arguments.
%
%      FOCUSED_PTS('Property','Value',...) creates a new FOCUSED_PTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Focused_Pts_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Focused_Pts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Focused_Pts

% Last Modified by GUIDE v2.5 22-Sep-2010 12:59:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Focused_Pts_OpeningFcn, ...
                   'gui_OutputFcn',  @Focused_Pts_OutputFcn, ...
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


% --- Executes just before Focused_Pts is made visible.
function Focused_Pts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Focused_Pts (see VARARGIN)

% Choose default command line output for Focused_Pts
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Focused_Pts wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Focused_Pts_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function cube_size_Callback(hObject, eventdata, handles)
% hObject    handle to cube_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cube_size as text
%        str2double(get(hObject,'String')) returns contents of cube_size as a double


% --- Executes during object creation, after setting all properties.
function cube_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cube_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in extract_focused_pts.
function extract_focused_pts_Callback(hObject, eventdata, handles)
% hObject    handle to extract_focused_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Az bounce_rem;

if size(Az,1) == 0
    errordlg('Please load input data first.','Reflection Contributions');
else

    % Get size of cubes
    cube_size = str2double(get(handles.cube_size,'string'));

    % Bounce level of interest
    bofi = str2double(get(handles.bounce_of_interest,'string'));
    
    % Get name of figure
    height_ax = str2double(get(handles.height_axis,'String'));

    % Specular phase centers only?
    flag_spec = get(handles.Specular_Phase,'Value');

    % Get string containing transformation information
    sensor_POV = get(handles.Sensor_Geo,'string');

    % Get string containing transformation information
    trans_POV = get(handles.transformation_POV,'string');

    if isnan(cube_size) == 1 || cube_size <= 0

        msgbox('Please specify the cube size [meter].','Size of cubes','error');
    else

        if isnan(bofi) || bofi < 1 || bofi > bounce_rem

            msgbox('Please specify the bounce level of interest appropriately.','Bounce level','error');
        else
                            
            if height_ax ~= 2 && height_ax ~= 3
                msgbox('Please specify the height axis.','Height axis','error');
            else

                if isempty(sensor_POV) == 1
                    msgbox('Please specify the imaging geometry (sensor position, center of footprint).','Imaging geometry','error');
                else

                    % Write focused points to 3D model file
                    Write_Focused_Pts(cube_size,trans_POV,sensor_POV,flag_spec,bofi,height_ax);
                end
            end
        end
    end

end


% --- Executes on button press in help_transformation.
function help_transformation_Callback(hObject, eventdata, handles)
% hObject    handle to help_transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox({'Required only if the object model has been translated, rotated or scaled in POV-Ray (e.g. for defining the aspect angle with respect to the object).' ' ' 'Then: insert the transformation commands used in POV Ray Editor for transforming the 3D model scene in POV Ray geometry.' ' ' 'Syntax: x<a b c>' ' ' 'x: transformation method ("rotate", "translate", "scale").' 'a, b, c: rotation angles [degrees], scaling factors [dimensionless], shift values for translation [m]' ' ' 'Caution: the order of commands has to be the identical to the order used in POV Ray!' ' ' 'Separate the transformation commands by semicolon ";".'},'Transformation Parameters','help');


function transformation_POV_Callback(hObject, eventdata, handles)
% hObject    handle to transformation_POV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transformation_POV as text
%        str2double(get(hObject,'String')) returns contents of transformation_POV as a double


% --- Executes during object creation, after setting all properties.
function transformation_POV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transformation_POV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sensor_Geo_Callback(hObject, eventdata, handles)
% hObject    handle to Sensor_Geo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sensor_Geo as text
%        str2double(get(hObject,'String')) returns contents of Sensor_Geo as a double


% --- Executes during object creation, after setting all properties.
function Sensor_Geo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sensor_Geo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_sensor_geo.
function help_sensor_geo_Callback(hObject, eventdata, handles)
% hObject    handle to help_sensor_geo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox({'Insert the commands used in POV Ray editor for defining imaging geometry of virtual SAR sensor.' ' ' 'Syntax: x<a b c>' ' ' 'x: type of position ("location", "look_at"). Entry "location" corresponds to the sensor position, entry "look_at" defines the center of the footprint.' 'a, b, c: coordinates in POV Ray coordinate system [m]' ' ' 'Separate the transformation commands by semicolon ";".'},'Imaging Geometry','help');


% --- Executes on button press in Specular_Phase.
function Specular_Phase_Callback(hObject, eventdata, handles)
% hObject    handle to Specular_Phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Specular_Phase



function bounce_of_interest_Callback(hObject, eventdata, handles)
% hObject    handle to bounce_of_interest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bounce_of_interest as text
%        str2double(get(hObject,'String')) returns contents of bounce_of_interest as a double


% --- Executes during object creation, after setting all properties.
function bounce_of_interest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bounce_of_interest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Select_Fig.
function Select_Fig_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select figure
[str_fig,path] = uigetfile('*.fig','Select Figure');

if str_fig == 0
    
    set(handles.Fig_Model,'string','');
    
else

    % Entry in input mask
    set(handles.Fig_Model,'String',str_fig);
    
end



function height_axis_Callback(hObject, eventdata, handles)
% hObject    handle to height_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height_axis as text
%        str2double(get(hObject,'String')) returns contents of height_axis as a double


% --- Executes during object creation, after setting all properties.
function height_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_height_axis.
function help_height_axis_Callback(hObject, eventdata, handles)
% hObject    handle to help_height_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Specify the height axis in the coordinate system of POV-Ray object model: 2: y-axis, 3: z-axis','Height axis.','help');
