function Gen_Height_Map

% function-file for creating an image whose values indicate the maximum elevation
% coordinate for each image pixel

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------


% local paramters: none

% Output:
% - Height map of scatterers within reflectivity map
%   --> maximum elevation value for each image pixel

%--------------------------------------------------------------------------

% global variables

global Az Ra El a_rem r_rem r_geom ang az_min az_max ra_min ra_max Output_path;

% Az: azimuth coordinates of reflection contributions [m]
% Ra: range coordinates of -''- [m]
% El: elevation coordinates of -''- [m]
% Tr_L: bounce level of -''- [value between 1 and 5]
% bounce_rem: maximum bounce level to be analyzed [value between 1 and 5]
% a_rem: variable storing resolution in azimuth within reflectivity map [m]
% r_rem: variable storing resolution in slant-range/ground-range within reflectivity map [m]
% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% ang: incidence angle of simulated signal [rad]
% az_min, az_max: minimum/maximum in azimuth saved in a global variable to be used for tomographic analysis [m]
% ra_min, ra_max: minimum/maximum in slant-range/ground-range saved in a global variable to be for tomographic analysis [m]
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% Copy range information for current function

if r_geom == 0
    % Slant range geometry
    Ra_g = Ra;
else
    % Change to ground range geometry
    Ra_g = Ra.*(1/sin(ang)); % ground range
end

% Number of lines
s = size(Az(:,1),1);

%--------------------------------------------------------------------------

% 1.) Set image properties
% --> according to reflectivity map

min_a = az_min; % minimum in azimuth
max_a = az_max; % maximum in azimuth
min_r = ra_min; % minimum in range
max_r = ra_max; % maximum in range
azimuth_pixel = a_rem; % resolution in azimuth
range_pixel = r_rem; % resolution in elevation

%--------------------------------------------------------------------------

% Intervals on axes
diff_a = max_a - min_a; % interval azimuth
diff_r = max_r - min_r; % interval slant range

num_a = ceil(diff_a/azimuth_pixel); % azimuth
num_r = ceil(diff_r/range_pixel); % range

%--------------------------------------------------------------------------

% 2.) Predefine size of image matrices

%--------------------------------------------------------------------------

% Matrix containing elevation data

Elevation_Map = zeros(num_r,num_a,1);
Elevation_Map(:,:) = NaN;

%--------------------------------------------------------------------------

% 3.) Fill image matrix with height data

%--------------------------------------------------------------------------

% Define steps in elevation for obtaining a reference system

% Loop over all contributions
for i = 1 : 1 : s
    
    % Analyze contributions
    
    %----------------------------------------------------------------------
        
    % If contribution is located within defined margins
     if Az(i,1) >= min_a && Az(i,1) <= max_a
         if Ra_g(i,1) >= min_r && Ra_g(i,1) <= max_r
                      
            %----------------------------------------------------------------------

            % B.)Pixel center
            % --> rounded coordinate in azimuth + 0.5
            % --> rounded coordinate in range + 0.5

            %e.g. for first pixel: row_pix = 0.5, column_pix = 0.5
           
            row_pix = floor((Ra_g(i,1)-min_r)/range_pixel)+0.5;
            column_pix = floor((Az(i,1)-min_a)/azimuth_pixel)+0.5;   
            
            %----------------------------------------------------------------------

            % D.) Add contribution to pixels 

            % Round up pixel coordinates for obtaining matrix indices
            column_pix = round(column_pix);
            row_pix = round(row_pix);

            if row_pix <= num_r && column_pix <= num_a
                
                if isnan(Elevation_Map(row_pix,column_pix)) == 1
         
                    % input of first entry
                    Elevation_Map(row_pix,column_pix) = El(i,1);
                else
                    
                    % Compare elevation values
                    if Elevation_Map(row_pix,column_pix) < El(i,1)
                        
                        % update entry due to higher elevation value
                        Elevation_Map(row_pix,column_pix) = El(i,1);
                    end
                end             
            end
        end
     end  
end

%--------------------------------------------------------------------------

% 4.) Preparations before displaying images

% Size of images to be displayed
r_p = 1:1:size(Elevation_Map,1); % rows --> range
a_p = 1:1:size(Elevation_Map,2); % columns --> azimuth

%--------------------------------------------------------------------------

% 5.) Plot Elevation Map

%--------------------------------------------------------------------------

if exist([Output_path,'/Maps/Figures/Elevation_Map'],'dir') ~= 7
   mkdir([Output_path,'/Maps/Figures/Elevation_Map']);
   mkdir([Output_path,'/Maps/Frames/Elevation_Map']);
end

% Name of figure to be displayed
H_Text = 'Elevation Map';

% Histogram containing all scatterers
Height_Map(Elevation_Map,a_p,r_p,H_Text);



    
    