function varargout = Simulation_Data(varargin)
% SIMULATION_DATA M-file for Simulation_Data.fig
%      SIMULATION_DATA, by itself, creates a new SIMULATION_DATA or raises the existing
%      singleton*.
%
%      H = SIMULATION_DATA returns the handle to a new SIMULATION_DATA or the handle to
%      the existing singleton*.
%
%      SIMULATION_DATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATION_DATA.M with the given input arguments.
%
%      SIMULATION_DATA('Property','Value',...) creates a new SIMULATION_DATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Simulation_Data_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Simulation_Data_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Simulation_Data

% Last Modified by GUIDE v2.5 23-Nov-2011 12:45:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Simulation_Data_OpeningFcn, ...
                   'gui_OutputFcn',  @Simulation_Data_OutputFcn, ...
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


% --- Executes just before Simulation_Data is made visible.
function Simulation_Data_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Simulation_Data (see VARARGIN)

% Choose default command line output for Simulation_Data
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Simulation_Data wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.scene_center_x,'string','0');
set(handles.scene_center_y,'string','0');
set(handles.scene_center_z,'string','0');


% --- Outputs from this function are returned to the command line.
function varargout = Simulation_Data_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Sensor_Parameters.
function Sensor_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to Sensor_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SAR_image_path;
    
% Get height direction
height_dir = str2double(get(handles.height_axis,'string'));
    
if height_dir == 1 || height_dir == 2

    % Select figure
    [str_par,path] = uigetfile([SAR_image_path,'/*.xml'],'Select image meta file.');
        
    if str_par == 0

        % No input file selected
        set(handles.sensor_par,'string','');

    else
        
        % Angle of incidence to be interpolated?
        interpolate_flag = get(handles.Angle_Interp,'value');
        
        % Entry in input mask
        set(handles.sensor_par,'String',str_par);
        
        % Get scene center from input mask
        scene_center = zeros(1,3);
        scene_center(1) = str2double(get(handles.scene_center_x,'string'));
        scene_center(2) = str2double(get(handles.scene_center_y,'string'));
        scene_center(3) = str2double(get(handles.scene_center_z,'string'));
        
        if isnan(scene_center(1)) == 1 || isnan(scene_center(2)) == 1 || isnan(scene_center(3)) == 1
            errordlg('Scene center coordinates have to be numbers.','Scene center');  
        else

            if sum(isnan(scene_center)) == 3
                scene_center = zeros(1,3);
            end

            % Path for loading image meta file
            meta_file_path = [path,'/',str_par];

            % Scan TerraSAR-X image meta file
            TerraSAR_to_POV(meta_file_path,height_dir,scene_center,interpolate_flag);
        end
    end
else
        
    msgbox('Please specify the height axis.','help');
end


% --- Executes on button press in manual_changes.
function manual_changes_Callback(hObject, eventdata, handles)
% hObject    handle to manual_changes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global startdir;

% Get system architecture
get_architecture = computer('arch');

% Case: windows
if strcmp(get_architecture,'win32') == 1 || strcmp(get_architecture, 'win64')

    %startdir=pwd;
    cd('../POV_Ray/POV_reg');
    system('pvengine.exe');
    cd(startdir);
end
    
% Case: linux    
if strcmp(get_architecture,'glnx86') == 1 || strcmp(get_architecture, 'glnxa64')
    
    msgbox('Note that editing in the POV-Ray editor is not supported by the Linux version of POV-Ray. Hence, please edit the POV-Ray file in a Linux editor.','Editing','help');
    
end
    

% --- Executes on button press in Render_it.
function Render_it_Callback(hObject, eventdata, handles)
% hObject    handle to Render_it (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get system architecture
get_architecture = computer('arch');

% Case: windows
if strcmp(get_architecture,'win32') == 1 || strcmp(get_architecture, 'win64')

    startdir=pwd;
    cd('../POV_Ray/RaySAR');
    system('pvengine.exe');
    cd(startdir);
end

% Case: linux    
if strcmp(get_architecture,'glnx86') == 1 || strcmp(get_architecture, 'glnxa64')
    
    msgbox('Note that no POV-Ray editor can be opened for the Linux version of POV-Ray. Hence, please run the POV-Ray file in the terminal (example command: povray scenario_name.pov +W1000 +H1000).','Run Process','help');
    
end

    
% --- Executes on button press in Triangulate_Cloud.
function Triangulate_Cloud_Callback(hObject, eventdata, handles)
% hObject    handle to Triangulate_Cloud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Model_path SAR_image_path;

format long;

[str_cloud,path] = uigetfile([Model_path,'/*.txt'],'Point cloud (delimiter in coordinate row: , ; or space)');

if str_cloud == 0
    
    set(handles.cloud_name,'string','');
    
else
    
    M_box = msgbox('Point cloud model is being created. Depending on the number of points, this may take a moment :-).','Loading','help');
    
    % Entry in input mask
    set(handles.cloud_name,'String',str_cloud);
    
    % Load point cloud
    path_name = [path,str_cloud];
    point_cloud = load(path_name);
            
    x = point_cloud(:,1);
    y = point_cloud(:,2);
    z = point_cloud(:,3);
    clear point cloud
        
    %--------------------------------------------------------------

    % Delaunay triangulation
    tri = delaunay(x,y);
    cloud_size = size(tri,1); % no. of triangles
    
    %--------------------------------------------------------------
    
    if isempty(Model_path) == 0
    
        % Store triangles in POV-Ray file    
        fid = fopen([Model_path,'/Point_cloud_scene.pov'],'w');

        %----------------------------------------
        % Add header if available

        if exist([SAR_image_path,'/POV_Ray_Header.pov'],'file') == 2

            fid2 = fopen([SAR_image_path,'/POV_Ray_Header.pov'],'r'); 

            while 1

                tline = fgetl(fid2);
                if ~ischar(tline), break, end        

                % Store string in LiDAR POV-Ray model file
                fprintf(fid,'%s \n',tline);
            end

            fprintf(fid,'%s \n',' ');
        end
        %----------------------------------------

        fprintf(fid,'%s \n','//Insert triangles of DEM');
        fprintf(fid,'%s \n \n','#declare DEM = union{');

        for k = 1 : 1: cloud_size

            % Point indices
            p_1 = tri(k,1);
            p_2 = tri(k,2);
            p_3 = tri(k,3);

            % Store triangle in POV-format
            fprintf(fid,'%s \n','triangle {');
            fprintf(fid,'%s %s %s %s %s %s %s \n','<', num2str(x(p_1)),',',num2str(z(p_1)),',',num2str(y(p_1)),'>');
            fprintf(fid,'%s %s %s %s %s %s %s \n','<', num2str(x(p_2)),',',num2str(z(p_2)),',',num2str(y(p_2)),'>');
            fprintf(fid,'%s %s %s %s %s %s %s \n','<', num2str(x(p_3)),',',num2str(z(p_3)),',',num2str(y(p_3)),'>');
            fprintf(fid,'%s \n \n','}');
        end

        fprintf(fid,'%s \n \n','}');
        fprintf(fid,'%s \n','object{DEM');
        fprintf(fid,'%s \n','  pigment {color <1,1,1>}');
        fprintf(fid,'%s \n','  finish {reflection{0.3} ambient 0 diffuse 0.3}');
        fprintf(fid,'%s \n','}');

        close(M_box);      

        if exist([SAR_image_path,'/POV_Ray_Header.pov'],'file') == 2

             fclose(fid);
             fclose(fid2);
             msgbox('POV-Ray DEM has been stored in the model folder (signal source and sensor added).','DEM stored to folder...');
        else

             fclose(fid);
             msgbox('POV-Ray DEM has been stored in the model folder (note: signal source and sensor definition still missing). For enabling the POV-Ray header generation, the TerraSAR-X image meta file has to be scanned first.','DEM stored to folder...');   
        end
        
    else
        msgbox('Path to model destination folder does not exist, yet. Please specify path to model folder in the GUI RaySAR.','POV-Ray header');
    end
end


function scene_center_z_Callback(hObject, eventdata, handles)
% hObject    handle to scene_center_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scene_center_z as text
%        str2double(get(hObject,'String')) returns contents of scene_center_z as a double


% --- Executes during object creation, after setting all properties.
function scene_center_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scene_center_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scene_center_y_Callback(hObject, eventdata, handles)
% hObject    handle to scene_center_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scene_center_y as text
%        str2double(get(hObject,'String')) returns contents of scene_center_y as a double


% --- Executes during object creation, after setting all properties.
function scene_center_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scene_center_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scene_center_x_Callback(hObject, eventdata, handles)
% hObject    handle to scene_center_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scene_center_x as text
%        str2double(get(hObject,'String')) returns contents of scene_center_x as a double


% --- Executes during object creation, after setting all properties.
function scene_center_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scene_center_x (see GCBO)
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


% --- Executes on button press in Open_Header_File.
function Open_Header_File_Callback(hObject, eventdata, handles)
% hObject    handle to Open_Header_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open editor visualizing POV-Ray header

global SAR_image_path;

if exist([SAR_image_path,'/POV_Ray_Header.pov'],'file') == 2
    
    % Get system architecture
    get_architecture = computer('arch');
    
    % Case: windows
    if strcmp(get_architecture,'win32') == 1 || strcmp(get_architecture, 'win64')
    
        % Open windows editor
        startdir=pwd;
        cd(SAR_image_path);
        winopen('POV_Ray_Header.pov'); % Open file with tool assigned in Windows
        cd(startdir);
    end
    
    % Case: Linux
    if strcmp(get_architecture,'glnx86') == 1 || strcmp(get_architecture, 'glnxa64')

        % Open file in linux
        startdir=pwd;
        cd(SAR_image_path);
        system('gedit POV_Ray_Header.pov');
        cd(startdir);
    end

else
    
    msgbox('A POV-Ray header file is not available. Please scan the TerraSAR-X image meta file first.','help');
    
end


% --- Executes on button press in Angle_Interp.
function Angle_Interp_Callback(hObject, eventdata, handles)
% hObject    handle to Angle_Interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Angle_Interp


% --- Executes on button press in Info_Triangulation.
function Info_Triangulation_Callback(hObject, eventdata, handles)
% hObject    handle to Info_Triangulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Optional for tool "Orbit file to POV-Ray", generating a POV-Ray header based on a TerraSAR-X image meta file (xml). Make sure that the SAR image for the pixel selection is defined as follows: azimuth (time) increases from right to left; range (time) increases bottom up.','Angle Interpolation');


% --- Executes on button press in Open_Body.
function Open_Body_Callback(hObject, eventdata, handles)
% hObject    handle to Open_Body (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Model_path;

[str_body,path] = uigetfile([Model_path,'/*.pov'],'POV-Ray file body');

if str_body ~= 0
    
    % Get system architecture
    get_architecture = computer('arch');

    % Case: windows
    if strcmp(get_architecture,'win32') == 1 || strcmp(get_architecture, 'win64')
        
        % Open file in windows
        startdir=pwd;
        cd(Model_path);
        winopen(str_body); % Open file with tool assigned in Windows
        cd(startdir);
    
    end
    
    % Case: Linux
    if strcmp(get_architecture,'glnx86') == 1 || strcmp(get_architecture, 'glnxa64')

       % Open file in linux
        startdir=pwd;
        cd(Model_path);
        system(['gedit ' str_body]);
        cd(startdir); 
        
    end
        
end


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


% --- Executes on button press in angle_to_POV.
function angle_to_POV_Callback(hObject, eventdata, handles)
% hObject    handle to angle_to_POV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get height direction
height_dir = str2double(get(handles.height_axis,'string'));
    
if height_dir == 1 || height_dir == 2


    % Get scene center from input mask
    scene_center = zeros(1,3);
    scene_center(1) = str2double(get(handles.scene_center_x,'string'));
    scene_center(2) = str2double(get(handles.scene_center_y,'string'));
    scene_center(3) = str2double(get(handles.scene_center_z,'string'));
    
    if isnan(scene_center(1)) == 1 || isnan(scene_center(2)) == 1 || isnan(scene_center(3)) == 1
       errordlg('Scene center coordinates have to be numbers.','Scene center');  
    else

        if sum(isnan(scene_center)) == 3
            scene_center = zeros(1,3);
        end

        % Get local angle of incidence (entered manually in GUI)
        angle_local = str2double(get(handles.angle,'string'))*(pi/180); % [rad]

        if angle_local >= 0 && angle_local <= pi/2

            % Define POV-Ray header
            POV_header(angle_local,height_dir,scene_center);

        else
            msgbox('The parameter for the incidence angle has to be between 0 and 90 degrees.','error');
        end
    end
else
        
    msgbox('Please specify height axis first.','help');
end


function ra_sampling_Callback(hObject, eventdata, handles)
% hObject    handle to ra_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ra_sampling as text
%        str2double(get(hObject,'String')) returns contents of ra_sampling as a double


% --- Executes during object creation, after setting all properties.
function ra_sampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ra_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in el_sampling.
function el_sampling_Callback(hObject, eventdata, handles)
% hObject    handle to el_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get angle of incidence
inc_angle = str2double(get(handles.angle,'String'))*(pi/180);

% get range sampling
range_sampling = str2double(get(handles.ra_sampling,'String'));

if isnan(inc_angle) == 1
    
    msgbox('Please specify the angle of incidence in the input mask.','Incidence angle','error');
    
else
    
    if isnan(range_sampling) == 1
        
        msgbox('Please specify sampling in range.','Range sampling','error');
        
    else
        % calculate angle
        inv_angle = pi/2-inc_angle;

        % calculate elevation sampling
        elevation_sampling = tan(inv_angle)*range_sampling;

        % display value in window
        msgbox(['Minimum stepwidth along the elevation axis for ray tracing: ', num2str(elevation_sampling), ' m. Further alternatives are obtained by deviding the value by integers.'],'Sampling along elevation axis','none');
    end
end
