function Write_Focused_Pts(c_size,trans_POV,sensor_POV,flag_spec,bofi,height_ax)

% function file for writing the 3D distribution of signatures to a CAD file

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input:
% - c_size: size of cubes for representation [m]
% - trans_POV: string containing POV-Ray transformation/scaling information to be
%              compensated
% - sensor_POV: string containing information about the imaging geometry of
%              the virtual sensor (sensor position, center of footprint)
% - flag_spec: flag for displaying specular components; 0 --> display all phase
%              centers; 1 --> display specular components only
% - bofi: bounce of interest to be displayed [integer],
%              e.g. 3 corresponds to triple bounce
% - height_ax: integer indicating the height axis in POV-Ray object file: 2 -->
%              y-axis, 3 --> z-axis


% output: geometry file
% - geometry file representing the 3D distribution of signatures (Format:
%   Wavefront, .obj; note: add ending .obj to resulting file)

%--------------------------------------------------------------------------

format long;

% global parameters

global Az Ra El Tr_L Sp Output_path;

% Az: azimuth coordinates of reflection contributions [m]
% Ra: range coordinates of -''- [m]
% Tr_L: bounce level of -''- [value between 1 and 5]
% Sp: vector containing flags for specular direction (0 --> non specular, 1
%     --> specular) for each intensity contribution
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% 1.) Select data of interest

% Limit for removing identical points if phase centers are located within the same area in space
% --> motivation: reduce number of points to be displayed!

% Get indices of points of interest
if flag_spec == 0
    Ind_M = find(Tr_L == bofi); % Triple bounce
else
    Ind_M_temp1 = find(Tr_L == bofi); % Triple bounce
    Ind_M_temp2 = find(Sp == 1); % Specular
    Ind_M = intersect(Ind_M_temp1, Ind_M_temp2);
end

% Select data
Az_M = Az(Ind_M);
Ra_M = Ra(Ind_M);
El_M = El(Ind_M);
Tr_M = Tr_L(Ind_M);

%--------------------------------------------------------------------------

% Create matrix containing points of interest
M_ROI = zeros(length(Ind_M),4);
M_ROI(:,1) = Az_M;
M_ROI(:,2) = Ra_M;
M_ROI(:,3) = El_M;
M_ROI(:,4) = Tr_M;

% Round for millimeter
M_ROI = M_ROI.*1000;
M_ROI = round(M_ROI);
M_ROI = M_ROI./1000;

% Remove redundant information and capture number of occurrence

% at first: check number of entries per cell by resorting pixel list and determining
%           first and last appearance of each pixel coordinate in list of input points
% - sort ascending (idx_data)
[M_ROI_sorted, sort_idx] = sortrows(M_ROI,[1,2,3]);

% - check difference between unique(...,'first') and unique(...,'last')
[M_ROI_sorted_1, I_1st, J_1st]    = unique(M_ROI_sorted,'rows','first');
[M_ROI_sorted_2, I_last, J_last] = unique(M_ROI_sorted,'rows','last'); 
% - difference of 'I's + 1 gives hint to how many entries per "scatterer"
no_entries = I_last - I_1st + 1;

% Update matrix of signals
M_ROI = M_ROI_sorted_1;

clear M_ROI_sorted_1 M_ROI_sorted_2
%--------------------------------------------------------------------------

% 2.) Transformation of signal samples: Image coordinate system --> Model coordinate system

% set flag to 1
flag = 1;

% Interprete "sensor_POV" --> camera position, position of footprint center
[M_sens,flag]=Check_Sensor(sensor_POV,flag);

% Get coordinates
for i = 1:2
    
    % Sensor Position
    if M_sens(i,4) == 0
        Sens_Pos = M_sens(i,1:3);
    end
    
    % Center of 3D Site Scene (position sensor is looking at)
    if M_sens(i,4) == 1
        Sc_Center = M_sens(i,1:3);
    end
end

%--------------------------------------------------------------------------

% 2.1.) Refine nadir information 
% --> nadir pointing in positive or negative direction?

% Get distance between sensor and scene footprint center in nadir direction

if height_ax == 3

    dif_nad = Sens_Pos(3)-Sc_Center(3);

    if dif_nad < 0

        % Correct direction of nadir 
        nadir = -1;
    else
        nadir = 1; 
    end
end

if height_ax == 2

    dif_nad = Sens_Pos(2)-Sc_Center(2);

    if dif_nad < 0

        % Correct direction of nadir 
        nadir = -1;
    else
        nadir = 1; 
    end
end

%--------------------------------------------------------------------------

% 2.2.) Phase centers in world coordinate system
M_ROI_temp = phase_pos(M_ROI,Sens_Pos,Sc_Center,nadir,height_ax);
%--------------------------------------------------------------------------

% 3.) Remove rotations, translations and scaling performed in POV Ray Editor

% Get transformation parameters for compensating POV Ray transformations
if isempty(trans_POV) == 0
    [M_trans,flag] = Check_Trans(trans_POV,flag);
    
    % Conduct necessary transformations 
    for j = 1:size(M_trans,1)

        if M_trans(j,4) == 0
            
            % Rotation
            rot_ang = M_trans(j,1:3)*(pi/180)*(-1); % scaling factor
            M_ROI_temp = rotation_M(M_ROI_temp,rot_ang); 
        end

        if M_trans(j,4) == 1
                        
            % Translation
            transf = M_trans(j,1:3); % scaling factor
            M_ROI_temp = translate_M(M_ROI_temp,transf);  
        end

        if M_trans(j,4) == 2

            % Scaling
            sc_f = M_trans(j,1:3); % scaling factor
            M_ROI_temp = scale_M(M_ROI_temp,sc_f); 
        end
    end
end

M_ROI(:,1:3) = M_ROI_temp(:,1:3);
clear M_ROI_temp;

% Create string
if bofi == 1
   string_add = 'Single'; 
end

if bofi == 2
   string_add = 'Double';  
end

if bofi == 3
   string_add = 'Triple';  
end

if bofi == 4
   string_add = 'Fourfold';  
end

if bofi == 5
   string_add = 'Fivefold';  
end

%--------------------------------------------------------------------------

% Get new size
s_ROI = size(M_ROI,1);
%--------------------------------------------------------------------------

% 4.) Write signal samples to geometry file

if flag == 1

    if exist([Output_path,'/Maps/Focused_Pts'],'dir') ~= 7
       mkdir([Output_path,'/Maps/Focused_Pts']); 
    end

    if exist([Output_path,'/Maps/Focused_Pts/Focused_Pts.obj'],'file') == 2
       delete([Output_path,'/Maps/Focused_Pts/Focused_Pts.obj']); 
    end

    if exist([Output_path,'/Maps/Focused_Pts/Model_Color.mtl'],'file') == 2
        delete([Output_path,'/Maps/Focused_Pts/Model_Color.mtl']);
    end

    % Create new or recreate
    fid = fopen([Output_path,'/Maps/Focused_Pts/Cubes_' string_add '.obj'],'w');
    fid2 = fopen([Output_path,'/Maps/Focused_Pts/Model_Color.mtl'],'w');
    fid3 = fopen([Output_path,'/Maps/Focused_Pts/Phase_Centers_' string_add '.txt'],'w');
    
    % header
    fprintf(fid,'mtllib Model_Color.mtl\n');

    % Create material file
    fprintf(fid2,'newmtl Single_Bounce\n Ns 10\n Ka 0.000000 0.000000 1.000000\n Kd 0.000000 0.000000 1.000000\n Ks 0.000000 0.000000 1.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
    fprintf(fid2,'newmtl Double_Bounce\n Ns 10\n Ka 0.000000 1.000000 0.000000\n Kd 0.000000 1.000000 0.000000\n Ks 0.000000 1.000000 0.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
    fprintf(fid2,'newmtl Triple_Bounce\n Ns 10\n Ka 1.000000 0.000000 0.000000\n Kd 1.000000 0.000000 0.000000\n Ks 1.000000 0.000000 0.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
    fprintf(fid2,'newmtl Fourfold_Bounce\n Ns 10\n Ka 1.000000 0.000000 1.000000\n Kd 1.000000 0.000000 1.000000\n Ks 1.000000 0.000000 1.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
    fprintf(fid2,'newmtl Fivefold_Bounce\n Ns 10\n Ka 0.000000 1.000000 1.000000\n Kd 0.000000 1.000000 1.000000\n Ks 0.000000 1.000000 1.000000\n Ni 1.000000\n d 1.000000\n illum 2');
        
    for k = 1:s_ROI
        
        % Write point coordinates to file
        fprintf(fid3,'%f %s %f %s %f %s %d\n',M_ROI(k,1),';',M_ROI(k,3),';',M_ROI(k,2),';',no_entries(k));

        % limits on x-axis
        x_min = M_ROI(k,1)-c_size;
        x_max = M_ROI(k,1)+c_size;

        % limits on y-axis
        y_min = M_ROI(k,2)-c_size;
        y_max = M_ROI(k,2)+c_size;

        % limits on z-axis
        z_min = M_ROI(k,3)*(-1)-c_size; % -1 for adaptation to Wavefront coordinate system
        z_max = M_ROI(k,3)*(-1)+c_size; 

        m_ind = [1 2 6 5; 4 1 5 8; 4 3 7 8; 3 2 6 7; 4 3 2 1; 7 6 5 8]+(k-1)*8;

        % Writing process
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_min,z_max);
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_min,z_max);
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_max,z_max);
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_max,z_max);
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_min,z_min);
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_min,z_min);
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_max,z_min);
        fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_max,z_min);

        if M_ROI(k,4) == 1
            fprintf(fid,'usemtl Single_Bounce\n');
        end

        if M_ROI(k,4) == 2
            fprintf(fid,'usemtl Double_Bounce\n');
        end

        if M_ROI(k,4) == 3
            fprintf(fid,'usemtl Triple_Bounce\n');
        end

        if M_ROI(k,4) == 4
            fprintf(fid,'usemtl Fourfold_Bounce\n');
        end

        if M_ROI(k,4) == 5
            fprintf(fid,'usemtl Fivefold_Bounce\n');
        end

        fprintf(fid,'s off\n');
        fprintf(fid,'f %7.0d %7.0d %7.0d %7.0d\n',m_ind(1,1),m_ind(1,2),m_ind(1,3),m_ind(1,4));
        fprintf(fid,'f %7.0d %7.0d %7.0d %7.0d\n',m_ind(2,1),m_ind(2,2),m_ind(2,3),m_ind(2,4));
        fprintf(fid,'f %7.0d %7.0d %7.0d %7.0d\n',m_ind(3,1),m_ind(3,2),m_ind(3,3),m_ind(3,4));
        fprintf(fid,'f %7.0d %7.0d %7.0d %7.0d\n',m_ind(4,1),m_ind(4,2),m_ind(4,3),m_ind(4,4));
        fprintf(fid,'f %7.0d %7.0d %7.0d %7.0d\n',m_ind(5,1),m_ind(5,2),m_ind(5,3),m_ind(5,4));
        fprintf(fid,'f %7.0d %7.0d %7.0d %7.0d\n',m_ind(6,1),m_ind(6,2),m_ind(6,3),m_ind(6,4));    
    end

    fclose(fid);
    fclose(fid2);
    fclose(fid3);

    msgbox({'Export of 3D phase centers accomplished (see folder Focused_Pts including: .obj model file, .mtl material file, 3D phase center coordinates + number of phase center occurrence.'},'Export of 3D Phase Centers','Help');
end






