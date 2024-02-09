function TerraSAR_to_POV(path,height_dir,scene_center,interpolate_flag)

% function file for extracting satellite orbit and image parameters from orbit file

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% Input: 
% - path: path to TerraSAR-X image meta file
% - height_dir: direction of height axis (1 --> positive y, -1 --> negative y, 2 --> positive z, -2 --> negative z)
% - scene_center: 3D coordinates of scene center (POV-Ray coordinate system) 
% - interpolate_flag: flag for triggering the interpolation of the local incidence angle [0: no, 1: yes]

% Output: 
% - parameters: none
% - files: simulation parameters stored in file 'Sensor_Parameters.txt'
% - function: initiate function POV_header

%--------------------------------------------------------------------------

format long;

global SAR_image_path;

% Open TerraSAR-X log-file
fid = fopen(path,'r');

% Delete old parameter file if necessary
if exist([SAR_image_path,'/Sensor_Parameters.txt'],'file') == 2
     delete([SAR_image_path,'/Sensor_Parameters.txt']); 
end

% Create new parameter file
fid2 = fopen([SAR_image_path,'/Sensor_Parameters.txt'],'w');

% Test variables
% restriction: one search per parameter
test_1 = 0; test_2 = 0; test_3 = 0; test_4 = 0; test_5 = 0;
test_7 = 0; test_8 = 0; test_9 = 0; test_10 = 0; test_11 = 0; test_12 = 0;
test_13 = 0; test_14 = 0; test_15 = 0; test_16 = 0; test_17 = 0; test_18 = 0;

% Search block
flag_search_block = 0;

%--------------------------------------------------------------------------

% Scan TerraSAR-X log-file for sensor parameters

while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
               
    % Search for angle of incidence
    if test_1 == 0 && flag_search_block == 0
        flag_inc_angle = strfind(tline,'<incidenceAngle>'); 
        if isempty(flag_inc_angle) == 0
            
            % Set flag value to 1
            test_1 = 1;
            flag_search_block = 1; % Block further incidence angle search for this round
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            inc_angle_center = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for spacing in east (for GEC) or range (for SSC)
    if test_2 == 0
        flag_row_spacing = strfind(tline,'<rowSpacing'); 
        if isempty(flag_row_spacing) == 0
            
            % Set flag value to 1
            test_2 = 1;
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            row_spacing = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for spacing in north (for GEC) and range (for SSC)
    if test_3 == 0
        flag_column_spacing = strfind(tline,'<columnSpacing'); 
        if isempty(flag_column_spacing) == 0
            
            % Set flag value to 1
            test_3 = 1;
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            column_spacing = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for azimuth resolution
    if test_4 == 0
        flag_az_resolution = strfind(tline,'<azimuthResolution>'); 
        if isempty(flag_az_resolution) == 0
            
            % Set flag value to 1
            test_4 = 1;
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            az_resolution = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for ground range resolution
    if test_5 == 0
        flag_grra_resolution = strfind(tline,'<groundRangeResolution>'); 
        if isempty(flag_grra_resolution) == 0
            
            % Set flag value to 1
            test_5 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            ground_range_resolution = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for orbit direction
    if test_7 == 0
        flag_orbit_direction = strfind(tline,'<orbitDirection>'); 
        if isempty(flag_orbit_direction) == 0
            
            % Set flag value to 1
            test_7 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            heading_mode = tline_entry;
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for signal polarization
    if test_8 == 0
        flag_polarization = strfind(tline,'<polLayer>'); 
        if isempty(flag_polarization) == 0
            
            % Set flag value to 1
            test_8 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            polarization = tline_entry;
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for heading angle
    if test_9 == 0
        flag_heading = strfind(tline,'<headingAngle>'); 
        if isempty(flag_heading) == 0
            
            % Set flag value to 1
            test_9 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            heading_angle = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for slant range resolution
    if test_10 == 0
        flag_slant_range_resolution = strfind(tline,'<slantRangeResolution>'); 
        if isempty(flag_slant_range_resolution) == 0
            
            % Set flag value to 1
            test_10 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            slant_range_resolution = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for slant range spacing
    if test_11 < 2
        
        flag_slant_range_spacing = strfind(tline,'<slantRange>');
        
        % Take second hit (first hit refers to Quicklook in XML file :-)
        if isempty(flag_slant_range_spacing) == 0
            
            % Set flag value to 1
            test_11 = test_11 + 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            slant_range_spacing = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for time of data take
    if test_12 == 0
        flag_time = strfind(tline,'<timeUTC>'); 
        if isempty(flag_time) == 0
            
            % Set flag value to 1
            test_12 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            time_of_data_take = tline_entry;
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for product type
    if test_13 == 0
        flag_product_type = strfind(tline,'<productType>'); 
        if isempty(flag_product_type) == 0
            
            % Set flag value to 1
            test_13 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            product_type = tline_entry;
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for sensor velocity
    if test_14 == 0
        flag_sensor_velocity = strfind(tline,'<velocity>'); 
        if isempty(flag_sensor_velocity) == 0
            
            % Set flag value to 1
            test_14 = 1;
            
             % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            sensor_velocity = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Get incidence angles at image corners
    if test_15 == 0 && flag_search_block == 0
        flag_inc_c1_angle = strfind(tline,'<incidenceAngle>'); 
        
        if isempty(flag_inc_c1_angle) == 0
            
            % Set flag value to 1
            test_15 = 1; 
            flag_search_block = 1; % Block further incidence angle search for this round
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            inc_angle_c1 = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    if test_16 == 0 && flag_search_block == 0
        flag_inc_c2_angle = strfind(tline,'<incidenceAngle>'); 
        
        if isempty(flag_inc_c2_angle) == 0
            
            % Set flag value to 1
            test_16 = 1;
            flag_search_block = 1; % Block further incidence angle search for this round
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            inc_angle_c2 = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    if test_17 == 0 && flag_search_block == 0
        flag_inc_c3_angle = strfind(tline,'<incidenceAngle>'); 
        
        if isempty(flag_inc_c3_angle) == 0
            
            % Set flag value to 1
            test_17 = 1;
            flag_search_block = 1; % Block further incidence angle search for this round
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            inc_angle_c3 = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end
    
    % Search for angle of incidence
    if test_18 == 0 && flag_search_block == 0
        flag_inc_c4_angle = strfind(tline,'<incidenceAngle>'); 
        
        if isempty(flag_inc_c4_angle) == 0
            
            % Set flag value to 1
            test_18 = 1;
            
            % Extract substring
            str_start = strfind(tline,'>');
            str_end = strfind(tline,'<');
            
            % Extract number from line
            tline_entry = tline(str_start(1)+1:str_end(2)-1);
            
            % Read sensor parameter
            inc_angle_c4 = str2double(tline_entry);
            
            clear str_start str_end tline_entry
        end
    end

    % Reset search block
    flag_search_block = 0;

end

%--------------------------------------------------------------------------

try

    % Interpolation of local incidence angle
    % --> local incidence angle based on marking the object location in the SAR image

    if interpolate_flag == 1

        % Get pixel coordinates from SAR image

        % Select figure
        [str_image,path_image] = uigetfile([SAR_image_path,'/*.*'],'Please select SAR image.');   
        string_load = [path_image,'/',str_image];

        % Read image
        SAR_Image = imread(string_load);
        size_image = size(SAR_Image); % get size
        
        % Define coordinates of image corners
        corner_top_left = [0,0];
        corner_top_right = [size_image(2),0];
        corner_bottom_left = [0,size_image(1)];
        corner_bottom_right = [size_image(2),size_image(1)];

        % Open image
        figure;
        imagesc(SAR_Image); colormap('gray'); title('SAR Image for pixel selection.');
        
        % Define dummies
        pix_x = -1; pix_y = -1;
              
        while pix_x < 0 || pix_x > size_image(2) || pix_y < 0 || pix_y > size_image(1)
            [pix_x,pix_y] = ginput(1);
        end 
                
        close(gcf); % Close figure after pixel selection
        
        % Calculate distances to image corners
        dist_corner_4 = sqrt((pix_y-corner_top_left(2))^2 +(pix_x-corner_top_left(1))^2);
        dist_corner_2 = sqrt((pix_y-corner_top_right(2))^2 +(pix_x-corner_top_right(1))^2);
        dist_corner_3 = sqrt((pix_y-corner_bottom_left(2))^2 +(pix_x-corner_bottom_left(1))^2);
        dist_corner_1 = sqrt((pix_y-corner_bottom_right(2))^2 +(pix_x-corner_bottom_right(1))^2);
        
        % Total distance for normalization
        dist_total = dist_corner_1 + dist_corner_2 + dist_corner_3 + dist_corner_4;
        
        % Weights based on distance to corners
        proportion_1 = dist_total/dist_corner_1;
        proportion_2 = dist_total/dist_corner_2;
        proportion_3 = dist_total/dist_corner_3;
        proportion_4 = dist_total/dist_corner_4;
        
        % Total of proportions
        proportion_sum = proportion_1 + proportion_2 + proportion_3 + proportion_4;
        
        % Normalized weights
        weight_angle_1 = proportion_1/proportion_sum;
        weight_angle_2 = proportion_2/proportion_sum;
        weight_angle_3 = proportion_3/proportion_sum;
        weight_angle_4 = proportion_4/proportion_sum;
        
        % Interpolate incidence angle (weighted mean)
        angle_local = (inc_angle_c1 * weight_angle_1 + inc_angle_c2 * weight_angle_2 + inc_angle_c3 * weight_angle_3 + inc_angle_c4 * weight_angle_4)*(pi/180);

    else

        % Take incidence angle from scene center
        angle_local = inc_angle_center*(pi/180); % [rad]
    end

    %--------------------------------------------------------------------------

    % Write sensor parameters to text file
    parameter = 'Product type: ';
    fprintf(fid2,'%s %s \n',parameter,product_type);
    parameter = 'Heading mode: ';
    fprintf(fid2,'%s %s \n',parameter,heading_mode);
    parameter = 'Polarization: ';
    fprintf(fid2,'%s %s \n',parameter,polarization);
    parameter = 'Time of data take: ';
    fprintf(fid2,'%s %s \n',parameter,time_of_data_take);

    if interpolate_flag == 1
        parameter = 'Angle of incidence (interpolated from corners) [deg]: ';
        fprintf(fid2,'%s %f \n',parameter,round(angle_local*(180/pi),4));
    else
        parameter = 'Angle of incidence (scene center) [deg]: ';
        fprintf(fid2,'%s %f \n',parameter,round(angle_local*(180/pi),4));
    end
    parameter = 'Heading angle [deg]: ';
    fprintf(fid2,'%s %f \n',parameter,heading_angle);
        
    if strcmp(product_type(1:3),'GEC') == 1
        parameter = 'Pixel spacing in east [m]: ';
        fprintf(fid2,'%s %f \n',parameter,row_spacing);
    end
    
    if strcmp(product_type(1:3),'SSC') == 1
        parameter = 'Pixel spacing in azimuth [m] (calculated from orbit parameters): ';
        fprintf(fid2,'%s %f \n',parameter,column_spacing*sensor_velocity);
    elseif strcmp(product_type(1:3),'GEC') == 1
        parameter = 'Pixel spacing in north [m]: ';
        fprintf(fid2,'%s %f \n',parameter,column_spacing);
    else        
        parameter = 'Pixel spacing for unknown product and format: ';
        fprintf(fid2,'%s %d \n',parameter,column_spacing);
    end
    
    parameter = 'Spacing in slant range [m]: ';
    fprintf(fid2,'%s %f \n',parameter,slant_range_spacing);
    parameter = 'Azimuth resolution [m]: ';
    fprintf(fid2,'%s %f \n',parameter,az_resolution);
    parameter = 'Slant range resolution [m]: ';
    fprintf(fid2,'%s %f \n',parameter,slant_range_resolution);
    parameter = 'Ground range resolution [m]: ';
    fprintf(fid2,'%s %f \n \n',parameter,ground_range_resolution);

    parameter = 'Please note that the definition of TerraSAR-X xml-files may be subject to changes. Then, parameter export has to be adapted accordingly.';
    fprintf(fid2,'%s \n',parameter);

catch
    
    msgbox('Problems occured while interpreting the xml-file. Please check whether the file fulfills the standard TerraSAR-X file format.','Problem occured!');
    
end

% Close text files
fclose(fid);
fclose(fid2);

%--------------------------------------------------------------------------

% Create POV-Ray header (light source + sensor plane)
POV_header(angle_local,height_dir,scene_center);

%--------------------------------------------------------------------------