function POV_header(angle,height_dir,sc_center)

% function file for creating header of POV-Ray simulation file

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% Input: 
% - angle: angle of incidence [rad]
% - height_dir: direction of height axis in POV-Ray coordinate system
%   1: y-axis (positive)
%   -1: y-axis (negative)
%   2: z-axis (positive)
%   -2: z-axis (negative)
% - sc_center: position of scene center (x,y,z in POV-Ray coordinates) [m]

% Output: none (POV-Ray header saved to file 'POV_Ray_Header.pov')

%--------------------------------------------------------------------------

global SAR_image_path;

% Prepare parameters

% Default values
comp_horizontal = 500; % horizontal component of distance from object scene
plane_dimension = 100; % height and width of sensor plane

% Vertical component
comp_vertical = comp_horizontal/tan(angle);

% Prepare vectors
satellite_pos = zeros(1,3);
ground_normal = zeros(1,3);

%--------------------------------------------------------------------------

% Distinguish: height axis in POV-Ray

% Height axis points along y-axis in POV-Ray coordinate system
if height_dir == 1 || height_dir == -1
       
    % satellite height = positive value
    if height_dir == 1
    
        % Calculate satellite position
        satellite_pos(1) = sc_center(1);
        satellite_pos(2) = sc_center(2) + comp_vertical;
        satellite_pos(3) = sc_center(3) - comp_horizontal;
        
        % Ground normal vector
        ground_normal(2) = 1;
    end
        
    % satellite height = negative value    
    if height_dir == -1  
            
        % Calculate satellite position
        satellite_pos(1) = sc_center(1);
        satellite_pos(2) = sc_center(2) - comp_vertical;
        satellite_pos(3) = sc_center(3) - comp_horizontal;
        
        % Ground normal vector
        ground_normal(2) = -1;
    end    
end

% Height axis along z-axis
if height_dir == 2 || height_dir == -2
        
    % satellite height = positive value
    if height_dir == 2
    
        % Calculate satellite position
        satellite_pos(1) = sc_center(1);
        satellite_pos(2) = sc_center(2) - comp_horizontal;
        satellite_pos(3) = sc_center(3) + comp_vertical;
        
        % Ground normal vector
        ground_normal(3) = 1;
    end
        
    % satellite height = negative value    
    if height_dir == -2
        
        % Calculate satellite position
        satellite_pos(1) = sc_center(1);
        satellite_pos(2) = sc_center(2) - comp_horizontal;
        satellite_pos(3) = sc_center(3) - comp_vertical;
        
        % Ground normal vector
        ground_normal(3) = -1;
    end
end

%--------------------------------------------------------------------------

if isempty(SAR_image_path) == 0
    
    % Write header to file
    fid = fopen([SAR_image_path,'/POV_Ray_Header.pov'],'w');

    fprintf(fid,'%s \n','#include "colors.inc"');
    fprintf(fid,'%s \n \n','#include "finish.inc"');
    fprintf(fid,'%s \n','// Global RaySAR simulation settings');
    fprintf(fid,'%s \n','// --> SAR_Output_Data: get data for 3D simulation (3D signal position, intensity, bounce level, specular flag), [1: yes, 0: no]');
    fprintf(fid,'%s \n \n','// --> SAR_Intersection: get 3D coordinates of object-ray intersections, [1: yes, 0: no]');
    fprintf(fid,'%s \n \n','//global_settings{SAR_Output_Data 1 SAR_Intersection 1}'); 
    fprintf(fid,'%s %s%s \n','// Radar sensor (angle of incidence: ',num2str(angle*(180/pi)),' degrees)');
    fprintf(fid,'%s \n','#declare Cam = camera {');
    fprintf(fid,'%s \n','    orthographic');
    fprintf(fid,'%s %s %s %s %s %s %s \n','    location <',num2str(satellite_pos(1)),',',num2str(satellite_pos(2)),',',num2str(satellite_pos(3)),'>');
    fprintf(fid,'%s %s %s %s %s %s %s \n','    look_at <',num2str(sc_center(1)),',',num2str(sc_center(2)),',',num2str(sc_center(3)),'>');    
    fprintf(fid,'%s %s %s \n','    right ', num2str(plane_dimension), '*x  // width of sensor plane');
    fprintf(fid,'%s %s %s \n','    up ', num2str(plane_dimension), '*y  // height of sensor plane');
    fprintf(fid,'%s \n \n','}');
    fprintf(fid,'%s \n \n','camera{Cam} ');
    fprintf(fid,'%s \n','// Radar antenna');
    fprintf(fid,'%s \n','light_source {');
    fprintf(fid,'%s \n','    0*x  // position of signal source');
    fprintf(fid,'%s \n','    color rgb <1,1,1> // color of emitted signal');
    fprintf(fid,'%s \n','    parallel');
    fprintf(fid,'%s %s %s %s %s %s %s \n','    translate <',num2str(satellite_pos(1)),',',num2str(satellite_pos(2)),',',num2str(satellite_pos(3)),'> // shift to antenna position');
    fprintf(fid,'%s %s %s %s %s %s %s \n','    point_at <',num2str(sc_center(1)),',',num2str(sc_center(2)),',',num2str(sc_center(3)),'> // point for defining signal direction');
    fprintf(fid,'%s \n \n','}');
    fprintf(fid,'%s \n \n','//--------------------------------------------------');                                                                  
    fprintf(fid,'%s \n','// Insert ground plane');
    fprintf(fid,'%s \n','#declare ground = plane {');
    fprintf(fid,'%s %s %s %s %s %s %s \n','    <',num2str(ground_normal(1)),',',num2str(ground_normal(2)),',',num2str(ground_normal(3)),'> // unit surface normal');
    fprintf(fid,'%s \n','    0 // distance from the origin in the direction of the surface normal');    
    fprintf(fid,'%s \n','}');  

    fclose(fid);

    msgbox('POV-Ray header file stored in SAR image folder.','POV-Ray header');
    
else 
    msgbox('Path to destination folder does not exist, yet. Please specify path to SAR image folder in the GUI RaySAR.','POV-Ray header');
end

