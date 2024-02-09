function Gen_Histo(el_res)

% function-file for generating scatterer histogram (number of scatterers per pixel)

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------


% local paramters:

% Input: 
% - el_res: resolution in elevation [m] --> corresponds to interval in
%           elevation direction where two scatterers are considered as one

% Output:
% - function H_Map is iniciated

%--------------------------------------------------------------------------

% global variables

global Az Ra El Tr_L bounce_rem a_rem r_rem r_geom ang az_min az_max ra_min ra_max Output_path;

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

M_box = msgbox('A scatterer histogram is being generated. Please wait as this may take a moment.','Scatterer histogram');

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
bounce_level = bounce_rem; % maximum bounce level to be analyzed
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

% Separated matrices 

%-----------------------------------

% Histogram

Histogram = zeros(num_r,num_a); 

% cell array for capturing number of scatterers and elevation heights
Histogram_heights = cell(size(Histogram)); 
%-----------------------------------
 
% Single Bounce
 
Single = zeros(num_r,num_a); 
Single_heights = cell(size(Single));
                                 
%-----------------------------------

% Double Bounce

if bounce_level > 1
    Double = zeros(num_r,num_a); 
    Double_heights = cell(size(Double));
end

%-----------------------------------

% Triple Bounce

if bounce_level > 2
    Triple = zeros(num_r,num_a); 
    Triple_heights = cell(size(Triple));                                  
end

%-----------------------------------

% Fourfold Bounce

if bounce_level > 3
    Fourfold = zeros(num_r,num_a); 
    Fourfold_heights = cell(size(Fourfold));                                                                        
end

%-----------------------------------

% Fivefold Bounce

if bounce_level > 4
    Fivefold = zeros(num_r,num_a); 
    Fivefold_heights = cell(size(Fivefold));                                                                     
end

%--------------------------------------------------------------------------

% 3.) Fill image matrices with height data

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
                
                %----------------------------------------------------------
                
                % Single bounce
                if Tr_L(i,1) == 1
                
                    % First entry
                    if Single(row_pix,column_pix) == 0

                       % Increase number of scatterers by 1
                       Single(row_pix,column_pix) = Single(row_pix,column_pix)+1; 

                       % Remember elevation value (first entry)
                       Single_heights{row_pix,column_pix,1} = El(i,1);

                    % Further entries   
                    else                     

                       % index 
                       j = 1;
                       test_ind = 0;

                       % check: scatterer with new elevation height?                        
                       while j <= Single(row_pix,column_pix)

                           if abs(El(i,1)-Single_heights{row_pix,column_pix,j}) <= el_res 

                               % Change flag since elevation value already exists
                               test_ind = 1;

                           end

                           % increase index
                           j = j + 1;
                       end 

                       if test_ind == 0

                           % Increase number of scatterers by 1
                           Single(row_pix,column_pix) = Single(row_pix,column_pix)+1;

                           % Remember elevation value (further entries)
                           Single_heights{row_pix,column_pix,j} = El(i,1);
                       end                        
                    end
                end
                
                %----------------------------------------------------------
                
                % Double bounce
                if Tr_L(i,1) == 2
                
                    % First entry
                    if Double(row_pix,column_pix) == 0

                       % Increase number of scatterers by 1
                       Double(row_pix,column_pix) = Double(row_pix,column_pix)+1; 

                       % Remember elevation value (first entry)
                       Double_heights{row_pix,column_pix,1} = El(i,1);

                    % Further entries   
                    else                     

                       % index 
                       j = 1;
                       test_ind = 0;

                       % check: scatterer with new elevation height?                        
                       while j <= Double(row_pix,column_pix)

                           if abs(El(i,1)-Double_heights{row_pix,column_pix,j}) <= el_res 

                               % Change flag since elevation value already exists
                               test_ind = 1;

                           end

                           % increase index
                           j = j + 1;
                       end 

                       if test_ind == 0

                           % Increase number of scatterers by 1
                           Double(row_pix,column_pix) = Double(row_pix,column_pix)+1;

                           % Remember elevation value (further entries)
                           Double_heights{row_pix,column_pix,j} = El(i,1);
                       end                        
                    end 
                end
                
                %----------------------------------------------------------
                
                % Triple bounce
                if Tr_L(i,1) == 3
                
                    % First entry
                    if Triple(row_pix,column_pix) == 0

                       % Increase number of scatterers by 1
                       Triple(row_pix,column_pix) = Triple(row_pix,column_pix)+1; 

                       % Remember elevation value (first entry)
                       Triple_heights{row_pix,column_pix,1} = El(i,1);

                    % Further entries   
                    else                     

                       % index 
                       j = 1;
                       test_ind = 0;

                       % check: scatterer with new elevation height?                        
                       while j <= Triple(row_pix,column_pix)

                           if abs(El(i,1)-Triple_heights{row_pix,column_pix,j}) <= el_res 

                               % Change flag since elevation value already exists
                               test_ind = 1;

                           end

                           % increase index
                           j = j + 1;
                       end 

                       if test_ind == 0

                           % Increase number of scatterers by 1
                           Triple(row_pix,column_pix) = Triple(row_pix,column_pix)+1;

                           % Remember elevation value (further entries)
                           Triple_heights{row_pix,column_pix,j} = El(i,1);
                       end                        
                    end      
                end
                %----------------------------------------------------------
                
                % Fourfold bounce
                if Tr_L(i,1) == 4
                
                    % First entry
                    if Fourfold(row_pix,column_pix) == 0

                       % Increase number of scatterers by 1
                       Fourfold(row_pix,column_pix) = Fourfold(row_pix,column_pix)+1; 

                       % Remember elevation value (first entry)
                       Fourfold_heights{row_pix,column_pix,1} = El(i,1);

                    % Further entries   
                    else                     

                       % index 
                       j = 1;
                       test_ind = 0;

                       % check: scatterer with new elevation height?                        
                       while j <= Fourfold(row_pix,column_pix)

                           if abs(El(i,1)-Fourfold_heights{row_pix,column_pix,j}) <= el_res 

                               % Change flag since elevation value already exists
                               test_ind = 1;

                           end

                           % increase index
                           j = j + 1;
                       end 

                       if test_ind == 0

                           % Increase number of scatterers by 1
                           Fourfold(row_pix,column_pix) = Fourfold(row_pix,column_pix)+1;

                           % Remember elevation value (further entries)
                           Fourfold_heights{row_pix,column_pix,j} = El(i,1);
                       end                        
                    end   
                end
                %----------------------------------------------------------
                
                % Fivefold bounce
                if Tr_L(i,1) == 5
                
                    % First entry
                    if Fivefold(row_pix,column_pix) == 0

                       % Increase number of scatterers by 1
                       Fivefold(row_pix,column_pix) = Fivefold(row_pix,column_pix)+1; 

                       % Remember elevation value (first entry)
                       Fivefold_heights{row_pix,column_pix,1} = El(i,1);

                    % Further entries   
                    else                     

                       % index 
                       j = 1;
                       test_ind = 0;

                       % check: scatterer with new elevation height?                        
                       while j <= Fivefold(row_pix,column_pix)

                           if abs(El(i,1)-Fivefold_heights{row_pix,column_pix,j}) <= el_res 

                               % Change flag since elevation value already exists
                               test_ind = 1;

                           end

                           % increase index
                           j = j + 1;
                       end 

                       if test_ind == 0

                           % Increase number of scatterers by 1
                           Fivefold(row_pix,column_pix) = Fivefold(row_pix,column_pix)+1;

                           % Remember elevation value (further entries)
                           Fivefold_heights{row_pix,column_pix,j} = El(i,1);
                       end                        
                    end      
                end
                %----------------------------------------------------------
                
                % Histogram counting all scatterers
                
                % First entry
                if Histogram(row_pix,column_pix) == 0
                    
                   % Increase number of scatterers by 1
                   Histogram(row_pix,column_pix) = Histogram(row_pix,column_pix)+1; 
                   
                   % Remember elevation value (first entry)
                   Histogram_heights{row_pix,column_pix,1} = El(i,1);
                   
                % Further entries   
                else                     

                   % index 
                   j = 1;
                   test_ind = 0;

                   % check: scatterer with new elevation height?                        
                   while j <= Histogram(row_pix,column_pix)
                       
                       if abs(El(i,1)-Histogram_heights{row_pix,column_pix,j}) <= el_res 

                           % Change flag since elevation value already exists
                           test_ind = 1;

                       end

                       % increase index
                       j = j + 1;
                   end 

                   if test_ind == 0

                       % Increase number of scatterers by 1
                       Histogram(row_pix,column_pix) = Histogram(row_pix,column_pix)+1;
                                              
                       % Remember elevation value (further entries)
                       Histogram_heights{row_pix,column_pix,j} = El(i,1);
                   end                        
                end                 
            end
        end
     end  
end

%--------------------------------------------------------------------------

% Clear cell arrays

clear Single_heights Double_heights Triple_heights Fourfold_heights Fivefold_heights Histogram_heights

%--------------------------------------------------------------------------

% 4.) Preparations before displaying images

% Size of images to be displayed
r_p = 1:1:size(Histogram,1); % rows --> range
a_p = 1:1:size(Histogram,2); % columns --> azimuth

close(M_box);

%--------------------------------------------------------------------------

% 5.) Plot Histograms

%--------------------------------------------------------------------------

if exist([Output_path,'/Maps/Figures/Scatterer_Map'],'dir') ~= 7
    
   mkdir([Output_path,'/Maps/Figures/Scatterer_Map']); 
   mkdir([Output_path,'/Maps/Figures/Scatterer_Map/JPG']); 
   mkdir([Output_path,'/Maps/Figures/Scatterer_Map/FIG']); 
   mkdir([Output_path,'/Maps/Frames/Scatterer_Map']);
end


H_Text = 'Scatterers (Single Bounce)';

% Maximum for map containing all scatterers
max_Single = max(max(Single(:,:)));
% Histogram containing all scatterers
Histo_Map(Single,a_p,r_p,max_Single,H_Text);


if exist('Double','var') && max(max(Double(:,:))) > 0
    
    H_Text = 'Scatterers (Double Bounce)';
    
    % Maximum -''- for double bounce
    max_Double = max(max(Double(:,:))); 
    % Plot double bounce histogram
    Histo_Map(Double,a_p,r_p,max_Double,H_Text);
end

if exist('Triple','var') && max(max(Triple(:,:))) > 0
    
    H_Text = 'Scatterers (Triple Bounce)';
    
    % Maximum -''- for triple bounce
    max_Triple = max(max(Triple(:,:))); 
    % Plot triple bounce histogram
    Histo_Map(Triple,a_p,r_p,max_Triple,H_Text);
end

if exist('Fourfold','var') && max(max(Fourfold(:,:))) > 0
    
    H_Text = 'Scatterers (Fourfold Bounce)';
    
    % Maximum -''- for fourfold bounce
    max_Fourfold = max(max(Fourfold(:,:))); 
    % Plot fourfold bounce histogram
    Histo_Map(Fourfold,a_p,r_p,max_Fourfold,H_Text);
end

if exist('Fivefold','var') && max(max(Fivefold(:,:))) > 0
    
    H_Text = 'Scatterers (Fivefold Bounce)';
    
    % Maximum -''- for fivefold bounce
    max_Fivefold = max(max(Fivefold(:,:))); 
    % Plot fivefold bounce histogram
    Histo_Map(Fivefold,a_p,r_p,max_Fivefold,H_Text);
end

H_Text = 'Histogram (all scatterers)';

% Maximum for map containing all scatterers
scat_max = max(max(Histogram(:,:)));
% Histogram containing all scatterers
Histo_Map(Histogram,a_p,r_p,scat_max,H_Text);


    
    