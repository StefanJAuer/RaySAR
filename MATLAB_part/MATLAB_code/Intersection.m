function varargout = Intersection(varargin)
% INTERSECTION M-file for Intersection.fig
%      INTERSECTION, by itself, creates a new INTERSECTION or raises the existing
%      singleton*.
%
%      H = INTERSECTION returns the handle to a new INTERSECTION or the handle to
%      the existing singleton*.
%
%      INTERSECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERSECTION.M with the given input arguments.
%
%      INTERSECTION('Property','Value',...) creates a new INTERSECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Intersection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Intersection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Intersection

% Last Modified by GUIDE v2.5 16-Jun-2011 19:03:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Intersection_OpeningFcn, ...
                   'gui_OutputFcn',  @Intersection_OutputFcn, ...
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


% --- Executes just before Intersection is made visible.
function Intersection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Intersection (see VARARGIN)

% Choose default command line output for Intersection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Intersection wait for user response (see UIRESUME)
% uiwait(handles.Intersection_Pts);


% --- Outputs from this function are returned to the command line.
function varargout = Intersection_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in extract_inters.
function extract_inters_Callback(hObject, eventdata, handles)
% hObject    handle to extract_inters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Index_All All_reflections az_min ra_min a_rem r_rem;
    
% Check if reflectivity map exists
if isempty(All_reflections) == 1
    
    errordlg('The reflectivity map has to be simulated first.','Lack of input data');
    
else  
    
    try
       
        arg1 = 1; % number of samples fixed to 1

        % Get pixel spacing in azimuth
        a_res = a_rem;

        % Get pixel spacing in slant/ground range
        r_res = r_rem;

        % Number of image pixels 
        num_pix_a = size(All_reflections,2); % get number of azimuth pixels
        num_pix_r = size(All_reflections,1); % get number of range pixels

        % Get size of cubes
        cube_size = str2double(get(handles.cube_size,'string'));

        % Get string containing transformation information
        trans_POV = get(handles.transformation_POV,'string');

        if isnan(cube_size) == 1 || cube_size <= 0

            msgbox('Please specify the cube size [meter].','Size of cubes','warn');

        else

                % Activate appropriate figure
                figure(Index_All);

                % Select Point
                ginput_Intersect(arg1,az_min,ra_min,a_res,r_res,num_pix_a,num_pix_r,cube_size,trans_POV);

        end
    
    catch
        msgbox('Pixel selection not possible. Has the simulated SAR image (= result of 2D map simulation) been closed?','Image not available');
    end
end



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


% --- Executes on button press in help_transformation.
function help_transformation_Callback(hObject, eventdata, handles)
% hObject    handle to help_transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox({'Required only if the object model has been translated, rotated or scaled in POV-Ray.' ' ' 'Then: insert the transformation commands used in POV Ray Editor for transforming the 3D model scene into the POV Ray geometry.' ' ' 'Syntax: x<a b c>' ' ' 'x: transformation method ("rotate", "translate", "scale").' 'a, b, c: rotation angles [degrees], scaling factors [dimensionless], shift values for the translation [m]' ' ' 'Caution: the order of commands has to be the identical to the order used in POV Ray!' ' ' 'Separate the transformation commands by semicolon ";".'},'Transformation Parameters','help');


% --- Executes when user attempts to close Intersection_Pts.
function Intersection_Pts_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Intersection_Pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(hObject);
