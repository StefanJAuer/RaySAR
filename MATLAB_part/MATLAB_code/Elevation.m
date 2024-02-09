function varargout = Elevation(varargin)
%ELEVATION M-file for Elevation.fig
%      ELEVATION, by itself, creates a new ELEVATION or raises the existing
%      singleton*.
%
%      H = ELEVATION returns the handle to a new ELEVATION or the handle to
%      the existing singleton*.
%
%      ELEVATION('Property','Value',...) creates a new ELEVATION using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Elevation_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ELEVATION('CALLBACK') and ELEVATION('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ELEVATION.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Elevation

% Last Modified by GUIDE v2.5 28-Jul-2010 13:43:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Elevation_OpeningFcn, ...
                   'gui_OutputFcn',  @Elevation_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before Elevation is made visible.
function Elevation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Elevation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Elevation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Elevation_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function elsamp_Callback(hObject, eventdata, handles)
% hObject    handle to elsamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elsamp as text
%        str2double(get(hObject,'String')) returns contents of elsamp as a double


% --- Executes during object creation, after setting all properties.
function elsamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elsamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function el_min_Callback(hObject, eventdata, handles)
% hObject    handle to el_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of el_min as text
%        str2double(get(hObject,'String')) returns contents of el_min as a double


% --- Executes during object creation, after setting all properties.
function el_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to el_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function el_max_Callback(hObject, eventdata, handles)
% hObject    handle to el_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of el_max as text
%        str2double(get(hObject,'String')) returns contents of el_max as a double


% --- Executes during object creation, after setting all properties.
function el_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to el_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Histogram_Elevation.
function Histogram_Elevation_Callback(hObject, eventdata, handles)
% hObject    handle to Histogram_Elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------

% global variables

global El;

% El: elevation coordinates of -''- [m]

%--------------------------------------------------------------------------

if size(El,1) == 0
    
    msgbox('No elevation data loaded!','Loading','warn');

else
    
    % Get resolution in elevation
    e_samp = str2num(get(handles.elsamp,'string'));
    
    % Get limits of histogram
    low_lim = round(min(El));
    high_lim = round(max(El));
    
    % Intervals for Histogram (along x-axis)
    el_interv = low_lim-2 : e_samp : high_lim+2;
       
    if isempty(e_samp) == 1 || e_samp <= 0
        uiwait(msgbox('Please specify the histogram intervals by entering the resolution in elevation!','Histogram steps','warn'));
        
    else
    
        figure;
        hist(El,el_interv); % Histogram
        title('Histogram of elevation values');
        xlabel('Elevation [m]'); ylabel('Number of entries');
        set(gcf,'Name','Histogram in Elevation','Numbertitle','off')
        xlim([el_interv(1),el_interv(length(el_interv))]);   
        grid;
    end
end


% --- Executes on button press in sl_azimuth.
function sl_azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to sl_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sl_azimuth


% --- Executes on button press in sl_range.
function sl_range_Callback(hObject, eventdata, handles)
% hObject    handle to sl_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sl_range


% --- Executes on button press in sl_elevation.
function sl_elevation_Callback(hObject, eventdata, handles)
% hObject    handle to sl_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sl_elevation


% --- Executes on button press in elevation_analysis.
function elevation_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to elevation_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------

% global variables

global All_reflections Index_All check_az check_ra check_el check_geom az_min ra_min a_rem r_rem;

% All_reflections: reflectivity map in a global matrix to be used in other parts of the program
% Index_All: handle to figure displaying reflectivity map (needed for elevation profiles)
% check_az: display slice along slant range axis at fixed azimuth interval? no --> value 0, yes --> value 1
% check_ra: display slice along azimuth axis at fixed slant range interval? no --> value 0, yes --> value 1
% check_el: display slice along elevation axis? no --> value 0, yes --> value 1
% check_geom: way of displaying height information; 
%             value 0 --> elevation with respect to center of synthetic aperture in elevation direction [m]
%             value 1 --> height over ground [m]
% az_min: minimum in azimuth saved in a global variable to be used for tomographic analysis [m]
% ra_min: minimum in range saved in a global variable to be for tomographic analysis [m]
% a_rem: variable storing resolution in azimuth within reflectivity map [m]
% r_rem: variable storing resolution in slant-range/ground-range within reflectivity map [m]

%--------------------------------------------------------------------------

% Check if reflectivity map exists
if isempty(All_reflections) == 1
    
    errordlg('A reflectivity map has to be simulated first.','Lack of input data');
    
else
    
    % Angle of incidence (dummy value)
    angle = -1;
           
    % Check which plots have to be displayed
    check_az = get(handles.sl_azimuth,'Value'); % azimuth
    check_ra = get(handles.sl_range,'Value'); % range
    check_el = get(handles.sl_elevation,'Value'); % elevation
        
    if get(handles.tag_elevation,'Value') == 1
        check_geom = 0;
    end
    
    if get(handles.tag_elevation,'Value') == 0
        
        check_geom = 1;
        
        if ~isempty(get(handles.angle,'String'))
            % Angle of incidence
            angle = str2double(get(handles.angle,'String'))*(pi/180); % [rad]
        end
    end
        
    %----------------------------------------------------------------------
        
    if isempty(get(handles.elsamp,'String'))
        
        errordlg('Please specify the resolution in elevation direction.','Missing of data');
        
    else
    
        if check_geom == 1 && angle == -1
            errordlg('Please specify the angle of incidence.','Missing of data');
        else

            % 1.) Provide necessary data

            % A.) Intervals to be used for azimuth and range slices:

            arg1 = 1; % number of samples fixed to 1
            e_res = -1; % dummy value for resolution in elevation direction

            %----------------------------------------------
            % B.) Get image properties

            % a) Get pixel spacing in azimuth
            a_res = a_rem;

            % b) Get pixel spacing in slant/ground range
            r_res = r_rem;

            % c) Get spacing in elevation

            % if elevation slice is to be displayed    
            if check_el == 1
                e_res = str2double(get(handles.elsamp,'String'));
            end

            % d) Get limits in elevation
            lim_e_min = get(handles.el_min,'String');
            lim_e_max = get(handles.el_max,'String');

            %----------------------------------------------
            % C.) Number of image pixels 
            num_pix_a = size(All_reflections,2); % get number of azimuth pixels
            num_pix_r = size(All_reflections,1); % get number of range pixels

            %-------------------------------------------------------------------------
            % 2.) Start selecting points
            
            try
            
                % Activate appropriate figure
                A = figure(Index_All);
            
                % Select Point
                ginput_POV(arg1,az_min,ra_min,a_res,r_res,e_res,lim_e_min,lim_e_max,num_pix_a,num_pix_r,angle);
            
            catch
                msgbox('Pixel selection not possible. Has the simulated SAR image (= result of 2D map simulation) been closed?','Image not available');
            end
         end
    end
end


% --- Executes on button press in tag_elevation.
function tag_elevation_Callback(hObject, eventdata, handles)
% hObject    handle to tag_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_elevation

set(handles.tag_ground,'Value',0.0);
set(handles.tag_elevation,'Value',1.0);
set(handles.angle,'Enable','Off');


% --- Executes on button press in tag_ground.
function tag_ground_Callback(hObject, eventdata, handles)
% hObject    handle to tag_ground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_ground

set(handles.tag_elevation,'Value',0.0);
set(handles.tag_ground,'Value',1.0);
set(handles.angle,'Enable','On');


function angle_Callback(hObject, eventdata, handles)
% hObject    handle to angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle as text
%        str2double(get(hObject,'String')) returns contents of angle as a double


% --- Executes during object creation, after setting all properties.
function angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in scatter_hist.
function scatter_hist_Callback(hObject, eventdata, handles)
% hObject    handle to scatter_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------

global All_reflections;

% Check if reflectivity map exists
if isempty(All_reflections) == 1
    
    errordlg('A reflectivity map has to be simulated first.','Lack of input data');
    
else
        
    % Get resolution in elevation direction
    e_samp = str2double(get(handles.elsamp,'String'));
 
    % Check if input data avaiable
    if isnan(e_samp) || e_samp <= 0
        errordlg('Please specify the resolution in elevation direction.','Lack of input data');
    else
        
        % Create Histogram
        Gen_Histo(e_samp);
    end
end


% --- Executes on button press in Elevation_Map.
function Elevation_Map_Callback(hObject, eventdata, handles)
% hObject    handle to Elevation_Map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global All_reflections;

% Check if reflectivity map exists
if isempty(All_reflections) == 1
    
    errordlg('A reflectivity map has to be simulated first.','Lack of input data');
    
else
             
    % Create Histogram
    Gen_Height_Map;
end


% --- Executes on button press in res_help.
function res_help_Callback(hObject, eventdata, handles)
% hObject    handle to res_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox({'"Elevation sampling" is valid when providing profiles for selected pixels (extent of resolution cell in elevation).' ' ' '"Elevation resolution" is valid when providing scatterer histograms (distance between two spatially separable scatterers).'},'Application','help')
