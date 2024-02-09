function[M_Pixel] = phase_pos(M_ROI,Sens_Pos,Sc_Center,nadir,height_ax)

% function file for calculating (model) world coordinates of simulated signatures

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input:
% - M_ROI: azimuth, range, elevation coordinates of signal samples [m]
% - Sens_Pos: world coordinates of sensor (1x3 vector) [m]
% - Sc_Center: world coordinates of the center of footprint (1x3 vector) [m]
% - nadir: variable indicating nadir direction, i.e. nadir = -1 --> 
%          height coordinate of sensor is lower than height coordinate of center of
%          footprint --> height axis of model defined in negative
%          direction (alternative: nadir = 1)
% - height_ax: number indicating height axis in POV-Ray object file: 2 -->
%              y-axis, 3 --> z-axis

% output:
% - M_Pixel: position of signal samples in the world coordinate system of the simulated object model [m]

%--------------------------------------------------------------------------

if height_ax == 2
        
    M_Pixel_temp = M_ROI(:,1:3);
    M_Pixel_temp(:,3) = M_Pixel_temp(:,3)*(-1); % left to right system

    % Define local coordinate system
    % origin: center of scene
    % --> new sensor coordinates!
    Sc_Center_temp = [0 0 0];
    Sens_Pos_temp = Sens_Pos - Sc_Center;

    % Line of sight vector
    Vec_Dir = Sc_Center_temp-Sens_Pos_temp; % direction vector
    abs_Vec = sqrt(Vec_Dir(1)^2+Vec_Dir(2)^2+Vec_Dir(3)^2); % absolute
    Vec_Dir = Vec_Dir./abs_Vec; % normalize to unit vector

    % nadir vector
    if nadir < 0 % z_sensor < z_footprint
        nad = [0 0 1]; % pointing in positive direction
    else % z_sensor > z_footprint
        nad = [0 0 -1]; % -''- positive direction
    end

    % incidence angle
    ang_nad = acos(Vec_Dir*nad');
    ang_vec2 = [-(pi/2-ang_nad) 0 0];   

    % change of aspect angle
    ang_nor = atan2(Vec_Dir(1),Vec_Dir(2));
    ang_vec = [0 0 -ang_nor];

    % Rotation according to look angle
    M_Pixel_temp = rotation_M(M_Pixel_temp,ang_vec2);

    % Rotation according to aspect angle
    M_Pixel_temp = rotation_M(M_Pixel_temp,ang_vec);
        
    % Coordinates in local system
    M_Pixel_temp(:,1) =  Sens_Pos_temp(1) - M_Pixel_temp(:,1);
    M_Pixel_temp(:,2) =  Sens_Pos_temp(2) + M_Pixel_temp(:,2);
    M_Pixel_temp(:,3) =  Sens_Pos_temp(3) - M_Pixel_temp(:,3);
            
    % World coordinates
    M_Pixel(:,1) = M_Pixel_temp(:,1) + Sc_Center(1);
    M_Pixel(:,2) = M_Pixel_temp(:,2) + Sc_Center(2);
    M_Pixel(:,3) = M_Pixel_temp(:,3) + Sc_Center(3);
end

if height_ax == 3
        
    M_Pixel_temp = M_ROI(:,1:3);
    M_Pixel_temp(:,3) = M_Pixel_temp(:,3)*(-1); % left to right system

    % Define local coordinate system
    % origin: center of scene
    % --> new sensor coordinates!
    Sc_Center_temp = [0 0 0];
    Sens_Pos_temp = Sens_Pos - Sc_Center;

    % Line of sight vector
    Vec_Dir = Sc_Center_temp-Sens_Pos_temp; % direction vector
    abs_Vec = sqrt(Vec_Dir(1)^2+Vec_Dir(2)^2+Vec_Dir(3)^2); % absolute
    Vec_Dir = Vec_Dir./abs_Vec; % normalize to unit vector

    % nadir vector
    if nadir < 0 % z_sensor < z_footprint
        nad = [0 0 1]; % pointing in positive direction
    else % z_sensor > z_footprint
        nad = [0 0 -1]; % -''- positive direction
    end
          
    % incidence angle
    ang_nad = acos(Vec_Dir*nad');
    ang_vec2 = [-(pi/2-ang_nad) 0 0];

    % change of aspect angle
    ang_nor = atan2(Vec_Dir(1),Vec_Dir(2));
    ang_vec = [0 0 -ang_nor];
        
    % Rotation according to look angle
    M_Pixel_temp = rotation_M(M_Pixel_temp,ang_vec2);

    % Rotation according to aspect angle
    M_Pixel_temp = rotation_M(M_Pixel_temp,ang_vec);
    
    % Coordinates in local system
    M_Pixel_temp(:,1) =  Sens_Pos_temp(1) - M_Pixel_temp(:,1);
    M_Pixel_temp(:,2) =  Sens_Pos_temp(2) + M_Pixel_temp(:,2);
    M_Pixel_temp(:,3) =  Sens_Pos_temp(3) - M_Pixel_temp(:,3);
        
    % World coordinates
    M_Pixel(:,1) = M_Pixel_temp(:,1) + Sc_Center(1);
    M_Pixel(:,2) = M_Pixel_temp(:,2) + Sc_Center(2);
    M_Pixel(:,3) = M_Pixel_temp(:,3) + Sc_Center(3); 
   
end