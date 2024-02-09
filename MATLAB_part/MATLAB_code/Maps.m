function varargout = Maps(varargin)

%MAPS M-file for Maps.fig
%      MAPS, by itself, creates a new MAPS or raises the existing
%      singleton*.
%
%      H = MAPS returns the handle to a new MAPS or the handle to
%      the existing singleton*.
%
%      MAPS('Property','Value',...) creates a new MAPS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Maps_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MAPS('CALLBACK') and MAPS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MAPS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Maps

% Last Modified by GUIDE v2.5 16-Jun-2011 18:47:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Maps_OpeningFcn, ...
                   'gui_OutputFcn',  @Maps_OutputFcn, ...
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


% --- Executes just before Maps is made visible.
function Maps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

global Az str_fig startdir;

cd(startdir);

% Choose default command line output for Maps
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Maps wait for user response (see UIRESUME)
% uiwait(handles.Mask_Single_Image);

if isempty(Az) == 0

    set(handles.Data_Loaded,'String',[str_fig ' loaded.']);
    
    % Get number of loaded point samples
    if size(Az,1) > 10^6
        num_points_temp = (round((100*size(Az,1))/10^6))/100; % two digits after the comma to be rounded
        num_points_temp = strcat(num2str(num_points_temp),' Million');
        num_points = strvcat('Number of loaded points:',num_points_temp);
    else
        num_points = strvcat('Number of loaded points:',num2str(size(Az,1)));
    end

    set(handles.No_of_points,'String',num_points);
end

set(handles.Min_Azimuth,'String','Min');
set(handles.Max_Azimuth,'String','Max');
set(handles.Min_Range,'String','Min');
set(handles.Max_Range,'String','Max');
set(handles.dB_Min,'String','Min');
set(handles.dB_Max,'String','Max');
set(handles.clipping_thresh,'String','Max');
set(handles.Bounce_Level,'String','Max');


% --- Outputs from this function are returned to the command line.
function varargout = Maps_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function azimuth_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to azimuth_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of azimuth_spacing as text
%        str2double(get(hObject,'String')) returns contents of azimuth_spacing as a double


% --- Executes during object creation, after setting all properties.
function azimuth_spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to azimuth_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function range_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to range_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of range_spacing as text
%        str2double(get(hObject,'String')) returns contents of range_spacing as a double


% --- Executes during object creation, after setting all properties.
function range_spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to range_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function incidence_angle_Callback(hObject, eventdata, handles)
% hObject    handle to incidence_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of incidence_angle as text
%        str2double(get(hObject,'String')) returns contents of incidence_angle as a double


% --- Executes during object creation, after setting all properties.
function incidence_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to incidence_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Slant_Range.
function Slant_Range_Callback(hObject, eventdata, handles)
% hObject    handle to Slant_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Slant_Range

set(handles.Ground_Range,'Value',0.0);
set(handles.incidence_angle,'Enable','Off');
set(handles.incidence_angle,'BackgroundColor',[0.925 0.914 0.847]);

%--------------------------------------------------------------------------

% global variables

global r_geom;

% r_geom: image geometry in range direction; 0: slant range, 1: ground range

%--------------------------------------------------------------------------

r_geom = 0;

% --- Executes on button press in Ground_Range.
function Ground_Range_Callback(hObject, eventdata, handles)
% hObject    handle to Ground_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ground_Range

set(handles.Slant_Range,'Value',0.0);
set(handles.incidence_angle,'Enable','On');
set(handles.incidence_angle,'BackgroundColor',[1 1 1]);

%--------------------------------------------------------------------------

% global variables

global r_geom;

% r_geom: image geometry in range direction; 0: slant range, 1: ground range

%--------------------------------------------------------------------------

r_geom = 1;

% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------

% global variables

global Az r_geom range_dir ang a_rem r_rem;

% Az: azimuth coordinates of reflection contributions [m]
% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% range_dir: range direction within simulated images; value 0: bottom up; value 1: top down
% ang: incidence angle of simulated signal [rad]
% a_rem, r_rem: global variables storing pixel spacing in azimuth and range

%--------------------------------------------------------------------------

if size(Az,1) == 0
    errordlg('Please load ray tracing data first (reflection contributions).','Reflection Contributions');
else
    %------------------------------------------
    % Range Geometry

    % Assume ground range
    r_geom = 1;

    % Check if slant range has to be used
    if get(handles.Slant_Range,'Value') == 1
        r_geom = 0;
    end
    %------------------------------------------

    % Range Direction: Top Down or Bottom Up?

    % Default values
    range_dir = 0; % Assume bottom up
    go.specular_flag = 0; % Assume: no specular analysis required
    go.sinc_sim = 0; % Assume: sinc not to be displayed
    go.a_res = -1; % resolution in azimuth
    go.r_res = -1; % resolution in range
    
    % Clipping or dB
    if get(handles.dB_range,'Value') == 1
        go.dBnodB = 1;
    else
        go.dBnodB = 0;
    end

    % Check whether "top down" has to be used
    if get(handles.Top_Down,'Value') == 1
        range_dir = 1;
    end
    
    % Check whether specular analysis has to be conducted
    if get(handles.specular,'Value') == 1
        go.specular_flag = 1;
    end
    
    if get(handles.impulse_resp,'Value') == 1
        go.sinc_sim = 1;
    end
        
    %------------------------------------------

    % Close old figures
    F = findobj('type','figure');

    for i = 1:length(F)
        if isempty(get(F(i),'Tag')) == 1
           close(F(i));
        end
    end

    % Default angle
    ang = -1;

    % Test flag for error inquiry    
    test_flag = 0;

    %----------------------------------------------------------------------
    % Detecting errors
    %----------------------------------------------------------------------

    % Range geometry

    % Requirements for ground range
    if r_geom == 1

        if isempty(get(handles.incidence_angle,'String')) == 1
            errordlg('Please enter value for incidence angle or change into slant range geometry.','Incidence Angle','help');
            test_flag = test_flag-1;
        end

        if isempty(get(handles.incidence_angle,'String')) == 0

            if str2double(get(handles.incidence_angle,'String')) <= 0
                errordlg('Incidence angle has to be bigger than zero.','Incidence Angle');
                test_flag = test_flag-1;
            end

            if str2double(get(handles.incidence_angle,'String')) > 90
                errordlg('Incidence angle has to be smaller than 90 degrees','Incidence Angle');
                test_flag = test_flag-1;
            end
            
            if isnan(str2double(get(handles.incidence_angle,'String'))) == 1
                errordlg('Incidence angle has to be a number.','Incidence Angle');
                test_flag = test_flag-1;
            end
        end
    end

    %----------------------------------------------------------------------

    % Azimuth spacing
    if isempty(get(handles.azimuth_spacing,'String')) == 1
        errordlg('Please enter value for azimuth spacing.','Azimuth Spacing','help');
        test_flag = test_flag-1;
    end

    if isempty(get(handles.azimuth_spacing,'String')) == 0
                        
        if isnan(str2double(get(handles.azimuth_spacing,'String'))) == 1
            errordlg('Azimuth spacing has to be a number.','Azimuth Spacing');
            test_flag = test_flag-1;
        end

        if str2double(get(handles.azimuth_spacing,'String')) <= 0
            errordlg('Azimuth spacing has to be bigger than zero.','Azimuth Spacing'); 
            test_flag = test_flag-1;
        end
    end

    %----------------------------------------------------------------------
    
    % Range spacing
    if isempty(get(handles.range_spacing,'String')) == 1
        errordlg('Please enter value for slant/ground range spacing.','Range Spacing','help');
        test_flag = test_flag-1;
    end

    if isempty(get(handles.range_spacing,'String')) == 0
                
        if isnan(str2double(get(handles.range_spacing,'String'))) == 1
            errordlg('Range spacing has to be a number.','Range Spacing');
            test_flag = test_flag-1;
        end

        if str2double(get(handles.range_spacing,'String')) <= 0
            errordlg('Range spacing has to be bigger than zero.','Range Spacing');
            test_flag = test_flag-1;
        end
    end
    
    %----------------------------------------------------------------------
    
    % Azimuth interval
    if isempty(get(handles.Min_Azimuth,'String')) == 1
         errordlg('Please enter minimum coordinate in azimuth.','Azimuth Interval','help');
         test_flag = test_flag-1;
    end
    
    if isempty(get(handles.Min_Azimuth,'String')) == 0
    
        if isnan(str2double(get(handles.Min_Azimuth,'String'))) == 1 && strcmp(get(handles.Min_Azimuth,'String'),'Min') == 0
            errordlg('Minimum azimuth coordinate has to be a number or Min.','Azimuth Interval');
            test_flag = test_flag-1;
        end
    end
    
    if isempty(get(handles.Max_Azimuth,'String')) == 1
        errordlg('Please enter maximum coordinate in azimuth.','Azimuth Interval','help');
        test_flag = test_flag-1;
    end
    
    if isempty(get(handles.Max_Azimuth,'String')) == 0
    
        if isnan(str2double(get(handles.Max_Azimuth,'String'))) == 1 && strcmp(get(handles.Max_Azimuth,'String'),'Max') == 0
            errordlg('Maximum azimuth coordinate has to be a number or Max.','Azimuth Interval'); 
            test_flag = test_flag-1;
        end
    end
    
    %----------------------------------------------------------------------
    
    % Range interval
    if isempty(get(handles.Min_Range,'String')) == 1
        errordlg('Please enter minimum coordinate in range/ground range.','Azimuth Interval','help');
        test_flag = test_flag-1;
    end
    
    if isempty(get(handles.Min_Range,'String')) == 0
    
        if isnan(str2double(get(handles.Min_Range,'String'))) == 1 && strcmp(get(handles.Min_Range,'String'),'Min') == 0
            errordlg('Minimum range coordinate has to be a number or Min.','Range Interval');
            test_flag = test_flag-1;
        end
    end
    
    if isempty(get(handles.Max_Range,'String')) == 1
        errordlg('Please enter maximum coordinate in range/ground range.','Azimuth Interval','help');
        test_flag = test_flag-1;
    end
    
    if isempty(get(handles.Max_Range,'String')) == 0
    
        if isnan(str2double(get(handles.Max_Range,'String'))) == 1 && strcmp(get(handles.Max_Range,'String'),'Max') == 0
            errordlg('Maximum range coordinate has to be a number or Max.','Range Interval');
            test_flag = test_flag-1;
        end
    end
    
    %----------------------------------------------------------------------
    
    % Bounce level
    if isempty(get(handles.Bounce_Level,'String')) == 1
        errordlg('Please enter the number of bounces to be displayed (maximum: 5).','Number of Bounces','help');
        test_flag = test_flag-1;
    end

    if isempty(get(handles.Bounce_Level,'String')) == 0
        
        if isnan(str2double(get(handles.Bounce_Level,'String'))) == 1 && strcmp(get(handles.Bounce_Level,'String'),'Max') == 0
            errordlg('The number of bounces has to be a number or Max.','Number of Bounces');
            test_flag = test_flag-1;
        end

        if str2double(get(handles.Bounce_Level,'String')) <= 0
            errordlg('The number of bounces has to be higher than zero.','Number of Bounces');
            test_flag = test_flag-1;
        end

        if str2double(get(handles.Bounce_Level,'String')) > 5
            errordlg('The number of bounces has to be smaller than 5.','Number of Bounces');
            test_flag = test_flag-1;
        end
    end
    
    %----------------------------------------------------------------------
    
    % dB vs. clipping
    if go.dBnodB == 1
        
        % dB
        if isempty(get(handles.dB_Min,'String')) == 1
            errordlg('Please enter value for minimum dB level ("Min" or number).','dB limit','help');
            test_flag = test_flag-1;
        end
        
        if isempty(get(handles.dB_Min,'String')) == 0
    
            if isnan(str2double(get(handles.dB_Min,'String'))) == 1 && strcmp(get(handles.dB_Min,'String'),'Min') == 0
                errordlg('Minimum dB level has to be a number or Min.','dB limit');
                test_flag = test_flag-1;
            end
        end

        if isempty(get(handles.dB_Max,'String')) == 1
            errordlg('Please enter value for maximum dB level ("Max" or number).','dB limit','help');
            test_flag = test_flag-1;
        end
        
        if isempty(get(handles.dB_Max,'String')) == 0
    
            if isnan(str2double(get(handles.dB_Max,'String'))) == 1 && strcmp(get(handles.dB_Max,'String'),'Max') == 0
                errordlg('Maximum dB level has to be a number or Max.','dB limit');
                test_flag = test_flag-1;
            end
        end
    
    else
        
        % clipping
        if isempty(get(handles.clipping_thresh,'String')) == 1
            errordlg('Please enter value for clipping threshold ("Max" or number).','Clipping threshold','help');
            test_flag = test_flag-1;
        end
        
        if isempty(get(handles.clipping_thresh,'String')) == 0
    
            if isnan(str2double(get(handles.clipping_thresh,'String'))) == 1 && strcmp(get(handles.clipping_thresh,'String'),'Max') == 0
                errordlg('The clipping threshold has to be a number or Max.','Clipping threshold');
                test_flag = test_flag-1;
            end
        end
    end
    
    if get(handles.add_noise,'Value') == 1
        
        % phase noise
        if isempty(get(handles.phase_angle,'String')) == 1
            errordlg('Please enter value for phase noise (0-90 degrees).','Phase noise','help');
            test_flag = test_flag-1;
        end

        if isempty(get(handles.phase_angle,'String')) == 0

            if isnan(str2double(get(handles.phase_angle,'String'))) == 1
                errordlg('The phase noise value has to be a number.','Phase noise');
                test_flag = test_flag-1;
            end
        end 
    
    end
            
    %----------------------------------------------------------------------
    
    % Initiate simulation process
    
    if test_flag == 0
        
        % Case: ground range --> otherwise: ang = -1 (default)
        if r_geom == 1

            % Get incidence angle
            ang = str2double(get(handles.incidence_angle,'string'));
            ang = ang*(pi/180); % [rad]

        end

        % Get pixel spacing in azimuth
        go.a_pix = str2double(get(handles.azimuth_spacing,'string'));

        % Get pixel spacing in slant/ground range
        go.r_pix = str2double(get(handles.range_spacing,'string'));

        % Store pixel spacing in global variables
        a_rem = go.a_pix;
        r_rem = go.r_pix;   

        % Get number of bounces
        go.bounce_level = get(handles.Bounce_Level,'string');

        if go.dBnodB == 1 
            % Power range (dB)
            go.db_min = get(handles.dB_Min,'string');
            go.db_max = get(handles.dB_Max,'string');
        else
            % Clipping threshold
            go.clip = get(handles.clipping_thresh,'string');
        end

        % minimum azimuth to be displayed
        go.a_min = get(handles.Min_Azimuth,'string');

        % maximum azimuth to be displayed
        go.a_max = get(handles.Max_Azimuth,'string');

        % minimum range to be displayed
        go.r_min = get(handles.Min_Range,'string');

        % maximum range to be displayed
        go.r_max = get(handles.Max_Range,'string');

        % get flag for smoothing process
        go.filt = get(handles.binomial,'Value');

        % get flag for coherent summation of signal contributions
        go.coh = get(handles.coherent,'Value');

        % get flag for adding noise
        go.noise = get(handles.add_noise,'Value');

        if go.noise == 1

            % Get maximum phase noise
            go.phase_noise = str2double(get(handles.phase_angle,'string'))*(pi/180);

        end

        % get flag for storing image matrix
        go.store = get(handles.store_image,'Value');

        % Check parameters of system impulse response
        if go.sinc_sim == 1

            % Get width of impulse response function
            go.a_res = str2double(get(handles.res_azimuth,'string'));
            go.r_res = str2double(get(handles.res_range,'string'));

            % Check whether resolution values are positive
            if go.a_res < 0 || go.r_res < 0

                errordlg('Width of impulse response has to be bigger than zero.','Resolution Parameters');

            % Check whether resolution entries are available
            elseif isnan(go.a_res) || isnan(go.r_res)

                 errordlg('Please choose resolution parameters appropriately.','Resolution Parameters');

            else
                % Create images
                Gen_Refl_Map(go);
            end

        else                            
            % Create images 
            Gen_Refl_Map(go);
        end
    end
end

%----------------------------------------------------------------------


function Bounce_Level_Callback(hObject, eventdata, handles)
% hObject    handle to Bounce_Level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bounce_Level as text
%        str2double(get(hObject,'String')) returns contents of Bounce_Level as a double


% --- Executes during object creation, after setting all properties.
function Bounce_Level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bounce_Level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Maximum.
function Maximum_Callback(hObject, eventdata, handles)
% hObject    handle to Maximum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------

% global variables

global Az Tr_L;

% Az: azimuth coordinates of reflection contributions [m]
% Tr_L: bounce level of -''- [value between 1 and 5]

%--------------------------------------------------------------------------

% Error, if no reflection contributions are available
if size(Az,1) == 0
    uiwait(errordlg('Please load ray tracing data first (reflection contributions).','Reflection Contributions'));
end

% If reflection contributions available
if size(Az,1) > 1
    Max_B = max(Tr_L);
    uiwait(msgbox(['Maximum number of bounces: ' num2str(Max_B)],'Maximum','help'));
end


% --- Executes on button press in contributions.
function contributions_Callback(hObject, eventdata, handles)
% hObject    handle to contributions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if exist('Az','var') == 1
    clear Az Ra El Intens Tr_L Sp;
end

%--------------------------------------------------------------------------

% global variables

global Az Ra El Intens Tr_L Sp S_X S_Y S_Z str_fig Model_path;

% Az: azimuth coordinates of reflection contributions [m]
% Ra: range coordinates of -''- [m]
% El: elevation coordinates of -''- [m]
% Intens: intensity of -''- [value between 0 and 1]
% Tr_L: bounce level of -''- [value between 1 and 5]
% Sp: flag for specular reflection; value 0: diffuse; value 1: specular contribution
% S_X: vector including x-coordinates of intersection points detected at objects
% S_Y: vector including y-coordinates of intersection points -''-
% S_Z: vector including z-coordinates of intersection points -''-
% str_fig: string containing name of input data

%--------------------------------------------------------------------------

% Close old figures
F = findobj('type','figure');

for i = 1:length(F)
    if isempty(get(F(i),'Tag')) == 1
       close(F(i));
    end
end

% Select figure
[str_fig,path] = uigetfile([Model_path,'/*.txt'],'Select ray tracing data');

if str_fig == 0
    
    set(handles.Data_Loaded,'string','');
    
else

    % Entry in input mask
    set(handles.Data_Loaded,'String',str_fig);

    M_box = msgbox('Please wait. Signal reflection data are loaded ...','Loading','help');

    % set default value for flag
    flag_load = 0;

    % check whether azimuth vector already has entries
    if size(Az,1) > 2

        % set flag to 1 --> entries in Maps are not changed during
        % loading process
        flag_load = 1;
    end
    
    % Full path for loading
    str_fig_full = [Model_path,'/',str_fig];

    % load
    R = load(str_fig_full);
    Ind_keep = 1:1:size(R,1);
         
    % Remove irrelevant data samples
    % --> i.e. intersection points without signal contribution 
    if size(R,2) > 6
        Ind_keep = Remove_irrelevant_data(R);
    end
   
    % Set Minimas and Maximas
    if  flag_load == 0;

        set(handles.Min_Azimuth,'String','Min');
        set(handles.Max_Azimuth,'String','Max');
        set(handles.Min_Range,'String','Min');
        set(handles.Max_Range,'String','Max');
        set(handles.dB_Min,'String','Min');
        set(handles.dB_Max,'String','Max');
        set(handles.clipping_thresh,'String','Max');
        set(handles.Bounce_Level,'String','Max');

    end

    % Store input data in vectors
    Az = R(Ind_keep,1); % azimuth coordinate
    Ra = R(Ind_keep,2); % slant range coordinate
    El = R(Ind_keep,3);
    Intens = R(Ind_keep,4); % intensities
    Tr_L = R(Ind_keep,5); % trace level
    Sp = R(Ind_keep,6); % flag for specular reflection

    % Check whether intersection points have been loaded
    if size(R,2) > 6

        % Capture coordinates of intersection points in vectors
        S_X = R(Ind_keep,7); % x-coordinates
        S_Y = R(Ind_keep,8); % y
        S_Z = R(Ind_keep,9); % z
    end
    
    close(M_box);
    uiwait(msgbox('Loading finished!','Loading','help'));
    set(handles.Data_Loaded,'String',[str_fig ' loaded.']);
    
    clear R

    % Get number of loaded point samples
    if size(Az,1) > 10^6
        num_points_temp = (round((100*size(Az,1))/10^6))/100; % two digits after the comma to be rounded
        num_points_temp = strcat(num2str(num_points_temp),' Million');
        num_points = strvcat('Number of loaded points:',num_points_temp);
    else
        num_points = strvcat('Number of loaded points:',num2str(size(Az,1)));
    end

    set(handles.No_of_points,'String',num_points);
end


function smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothing as text
%        str2double(get(hObject,'String')) returns contents of smoothing as a double


% --- Executes during object creation, after setting all properties.
function smoothing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Top_Down.
function Top_Down_Callback(hObject, eventdata, handles)
% hObject    handle to Top_Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Top_Down

%--------------------------------------------------------------------------

% global variables

global range_dir;

% range_dir: range direction within simulated images; value 0: bottom up; value 1: top down

%--------------------------------------------------------------------------

range_dir = get(handles.Top_Down,'Value');



function Max_Range_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Max_Range as text
%        str2double(get(hObject,'String')) returns contents of Max_Range as a double


% --- Executes during object creation, after setting all properties.
function Max_Range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Histogram_Range.
function Histogram_Range_Callback(hObject, eventdata, handles)
% hObject    handle to Histogram_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------

% global variables

global Ra;

% Ra: range coordinates of reflection contributions [m]

%--------------------------------------------------------------------------

if size(Ra,1) == 0
    errordlg('Please load ray tracing data first (signal contributions).','Loading');
end

% If Range values available
if size(Ra,1) > 0
    
    %------------------------------------------
    % Range Geometry

    % Assume ground range
    r_geom = 1;

    % Check if slant range has to be used
    if get(handles.Slant_Range,'Value') == 1
        r_geom = 0;
    end    
    %------------------------------------------

    % Case: ground range 
    %if get(handles.Ground_Range,'Value') == 1
    if r_geom == 1
        
        if isempty(get(handles.range_spacing,'string'))
            uiwait(msgbox('Please specify histogram intervals by entering range spacing.','Histogram steps','warn'));
        else
            if isempty(get(handles.incidence_angle,'string'))
                uiwait(msgbox('Please specify incidence angle or change into slant range geometry.','Angle','warn'));
            else
                
                % Get incidence angle
                ang = str2double(get(handles.incidence_angle,'string'))*(pi/180);
                
                % Get resolution in range
                r_res = str2double(get(handles.range_spacing,'string'));
                
                % Get limits of histogram
                low_lim = round(min(Ra))*(1/sin(ang));
                high_lim = round(max(Ra))*(1/sin(ang));
            
                % Intervals for Histogram (along x-axis)
                r_interv = low_lim-2 : r_res : high_lim+2;

                figure;
                hist(Ra.*(1/sin(ang)),r_interv); % Histogram
                title('Histogram of range values - ground range');
                xlabel('Range [m]'); ylabel('Number of entries');
                set(gcf,'Name','Histogram in Ground Range','Numbertitle','off');
                grid;
            end
        end
    end

    % Case: slant range 
    %if get(handles.Ground_Range,'Value') == 0
    if r_geom == 0
  
        % Get resolution in range
        r_res = str2double(get(handles.range_spacing,'string'));
        
        % Get limits of histogram
        low_lim = round(min(Ra));
        high_lim = round(max(Ra));
       
        % Intervals for Histogram (along x-axis)
        r_interv = low_lim-3 : r_res : high_lim+3;

        if isempty(r_res)
            uiwait(msgbox('Please specify histogram intervals by entering range spacing!','Histogram steps','warn'));
        else
            figure;
            hist(Ra,r_interv); % Histogram
            title('Histogram of range values - slant range');
            xlabel('Range [m]'); ylabel('Number of entries');
            set(gcf,'Name','Histogram in Slant Range','Numbertitle','off');
            grid;
        end
    end
end


% --- Executes when user attempts to close Mask_Single_Image.
function Mask_Single_Image_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Mask_Single_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(hObject);


function Range_Max_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Max_Range as text
%        str2double(get(hObject,'String')) returns contents of Max_Range as a double


% --- Executes during object creation, after setting all properties.
function Range_Max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Min_Azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to Min_Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Min_Azimuth as text
%        str2double(get(hObject,'String')) returns contents of Min_Azimuth as a double


% --- Executes during object creation, after setting all properties.
function Min_Azimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Histogramm_Azimuth.
function Histogramm_Azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to Histogramm_Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------

% global variables

global Az;

% Az: azimuth coordinates of reflection contributions [m]

%--------------------------------------------------------------------------

if size(Az,1) == 0
    errordlg('Please load ray tracing data first (signal contributions).','Loading');
end

if size(Az,1) > 0
    
    % Get resolution in slant/ground range
    a_res = str2double(get(handles.azimuth_spacing,'string'));
    
    if isnan(a_res) == 1
        
        % Warning due to missing parameter
         uiwait(msgbox('Please specify histogram interval by entering azimuth spacing','Histogram steps','warn'));
        
    else

        % Get limits of histogram
        low_lim = round(min(Az));
        high_lim = round(max(Az));

        % Intervals for Histogram (along x-axis)
        a_interv = low_lim-2 : a_res : high_lim+2;

        if ~isempty(a_res)
            figure;
            hist(Az,a_interv); % Histogram
            title('Histogram of azimuth values');
            xlabel('Azimuth [m]'); ylabel('Number of entries'); 
            set(gcf,'Name','Histogram in Azimuth','Numbertitle','off');
            grid;
        end
    end
end

function Min_Range_Callback(hObject, eventdata, handles)
% hObject    handle to Min_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Min_Range as text
%        str2double(get(hObject,'String')) returns contents of Min_Range as a double


% --- Executes during object creation, after setting all properties.
function Min_Range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Max_Azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Max_Azimuth as text
%        str2double(get(hObject,'String')) returns contents of Max_Azimuth as a double


% --- Executes during object creation, after setting all properties.
function Max_Azimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dB_Min_Callback(hObject, eventdata, handles)
% hObject    handle to dB_Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dB_Min as text
%        str2double(get(hObject,'String')) returns contents of dB_Min as a double


% --- Executes during object creation, after setting all properties.
function dB_Min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dB_Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dB_Max_Callback(hObject, eventdata, handles)
% hObject    handle to dB_Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dB_Max as text
%        str2double(get(hObject,'String')) returns contents of dB_Max as a double


% --- Executes during object creation, after setting all properties.
function dB_Max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dB_Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in binomial.
function binomial_Callback(hObject, eventdata, handles)
% hObject    handle to binomial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of binomial


% --- Executes on button press in specular.
function specular_Callback(hObject, eventdata, handles)
% hObject    handle to specular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of specular


% --- Executes on button press in SAR_Image.
function SAR_Image_Callback(hObject, eventdata, handles)
% hObject    handle to SAR_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in speckle.
function speckle_Callback(hObject, eventdata, handles)
% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speckle


% --- Executes on button press in impulse_resp.
function impulse_resp_Callback(hObject, eventdata, handles)
% hObject    handle to impulse_resp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of impulse_resp

if get(handles.impulse_resp,'Value') == 0
    
    set(handles.res_azimuth,'Enable','Off');
    set(handles.res_azimuth,'BackgroundColor',[0.925 0.914 0.847]);
    set(handles.res_range,'Enable','Off');
    set(handles.res_range,'BackgroundColor',[0.925 0.914 0.847]);
    
else

    set(handles.res_azimuth,'Enable','On');
    set(handles.res_azimuth,'BackgroundColor',[1 1 1]);
    set(handles.res_range,'Enable','On');
    set(handles.res_range,'BackgroundColor',[1 1 1]);

end


function res_azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to res_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_azimuth as text
%        str2double(get(hObject,'String')) returns contents of res_azimuth as a double


% --- Executes during object creation, after setting all properties.
function res_azimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function res_range_Callback(hObject, eventdata, handles)
% hObject    handle to res_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_range as text
%        str2double(get(hObject,'String')) returns contents of res_range as a double


% --- Executes during object creation, after setting all properties.
function res_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dB_range.
function dB_range_Callback(hObject, eventdata, handles)
% hObject    handle to dB_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dB_range

set(handles.clipping,'Value',0);
set(handles.clipping_thresh,'Enable','Off');
set(handles.clipping_thresh,'BackgroundColor',[0.925 0.914 0.847]);
set(handles.dB_Min,'Enable','On');
set(handles.dB_Min,'BackgroundColor',[1 1 1]);
set(handles.dB_Max,'Enable','On');
set(handles.dB_Max,'BackgroundColor',[1 1 1]);



% --- Executes on button press in clipping.
function clipping_Callback(hObject, eventdata, handles)
% hObject    handle to clipping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clipping

set(handles.dB_range,'Value',0);    
set(handles.dB_Min,'Enable','Off');
set(handles.dB_Min,'BackgroundColor',[0.925 0.914 0.847]);
set(handles.dB_Max,'Enable','Off');
set(handles.dB_Max,'BackgroundColor',[0.925 0.914 0.847]);
set(handles.clipping_thresh,'Enable','On');
set(handles.clipping_thresh,'BackgroundColor',[1 1 1]);


function clipping_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to clipping_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clipping_thresh as text
%        str2double(get(hObject,'String')) returns contents of clipping_thresh as a double


% --- Executes during object creation, after setting all properties.
function clipping_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clipping_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in coherent.
function coherent_Callback(hObject, eventdata, handles)
% hObject    handle to coherent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coherent

if get(handles.coherent,'Value') == 0
    
    set(handles.add_noise,'Enable','Off');
    set(handles.add_noise,'BackgroundColor',[0.925 0.914 0.847]);
    
else

    set(handles.add_noise,'Enable','On');
    set(handles.add_noise,'BackgroundColor',[1 1 1]);

end


% --- Executes on button press in add_noise.
function add_noise_Callback(hObject, eventdata, handles)
% hObject    handle to add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_noise

if get(handles.add_noise,'Value') == 0
    
    set(handles.phase_angle,'Enable','Off');
    set(handles.phase_angle,'BackgroundColor',[0.925 0.914 0.847]);
    
else

    set(handles.phase_angle,'Enable','On');
    set(handles.phase_angle,'BackgroundColor',[1 1 1]);

end


function phase_angle_Callback(hObject, eventdata, handles)
% hObject    handle to phase_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phase_angle as text
%        str2double(get(hObject,'String')) returns contents of phase_angle as a double


% --- Executes during object creation, after setting all properties.
function phase_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in store_image.
function store_image_Callback(hObject, eventdata, handles)
% hObject    handle to store_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of store_image
