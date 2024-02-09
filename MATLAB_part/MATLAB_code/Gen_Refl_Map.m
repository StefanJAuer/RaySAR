function Gen_Refl_Map(go)

% function-file for creating reflectivity maps, maps for specular analysis and intensity distribution maps

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters:

% Elements within struct "go": 
% - go.a_pix: azimuth pixel spacing of reflectivity map to be generated [m]
% - go.r_pix: range pixel spacing of reflectivity map to be generated [m]
% - go.a_min, go.a_max, go.r_min, go.r_max: image limits on azimuth and range axis [m]
% - go.db_min, go.db_max: power range to be displayed [dB in string format]
% - go.clip: threshold for clipping of image intensities [string]
% - go.filt: flag for applying binomial filter [no: value 0; yes: value 1]
% - go.bounce_level: maximum bounce level to be analyzed [value between 1 and 5]
% - go.specular_flag: flag for preparing specular maps [no --> value 0; yes --> value 1]
% - go.sinc_sim: flag for convolving the simulated image with the SAR 
%   system response [no --> value 0; yes --> value 1]
% - go.a_res,go.r_res: resolution in azimuth and range (3dB width of 2D impulse response)
% - go.coh: flag for coherent summation of signal contributions (0 = absolute values only; 1 = summation of complex values)
% - go.dBnodB: flag for deciding between dB amplitude scale and clipping [value 0: clipping, value 1: dB]
% - go.store: flag for storing image data to file [value 0: do not store pixel values, value 1: store image matrix to .mat file]

% Output:
% - geometrical distribution of scatterers in the azimuth-range plane
% - reflectivity map and separate layers for different reflection levels
%   (clipped or dB)
% - specular map: binary map marking pixels containing specular signal
%   contributions (angle tolerance: +-1 degrees)

%--------------------------------------------------------------------------

% global parameters

global Az Ra Intens Tr_L Sp ang r_geom range_dir All_reflections Index_All az_min az_max ra_min ra_max bounce_rem Output_path;

% Az: azimuth coordinates of reflection contributions [m]
% Ra: range coordinates of -''- [m]
% Intens: intensity of -''- [value between 0 and 1]
% Tr_L: bounce level of -''- [value between 1 and 5]
% Sp: flag for specular reflection; value 0: diffuse; value 1: specular contribution
% ang: incidence angle of simulated signal [rad]
% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% All_reflections: save reflectivity map in a global matrix to be used in other parts of the program
% Index_All: handle to figure displaying reflectivity map (to be used for tomographic analysis)
% az_min, az_max: minimum/maximum in azimuth saved in a global variable to be used for tomographic analysis [m]
% ra_min, ra_max: minimum/maximum in slant-range/ground-range saved in a global variable to be for tomographic analysis [m]
% bounce_rem: maximum bounce level to be used for tomographic analysis 
%             --> corresponds to maximum value used for creating the reflectivity map)
%             [value between 1 and 5]
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% Delete older result folder

if exist([Output_path,'/Maps'],'dir') == 7
    rmdir([Output_path,'/Maps'],'s');
end

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

% Signal wavelength (X-Band)
lambda = 0.031; % [m]

% Signal wavelength (Ku-Band)
%lambda = 0.008543; % [m]

%--------------------------------------------------------------------------

% 1.) Set image properties

%--------------------------------------------------------------------------

% Set boundaries of image

% minimum in azimuth
if strcmp(go.a_min,'Min') 
     go.a_min = round(min(Az(:,1))); 
else 
    go.a_min = str2double(go.a_min);
end

% maximum in azimuth
if strcmp(go.a_max,'Max')
     go.a_max = round(max(Az(:,1))); 
else
    go.a_max = str2double(go.a_max);
end

% minimum in range
if strcmp(go.r_min,'Min') 
     go.r_min = round(min(Ra_g(:,1))); 
else
    go.r_min = str2double(go.r_min);
end

% maximum in range
if strcmp(go.r_max,'Max') 
     go.r_max = round(max(Ra_g(:,1))); 
else
    go.r_max = str2double(go.r_max);
end

% Maximum bounce level to be displayed
if strcmp(go.bounce_level,'Max') 
    go.bounce_level = max(Tr_L(:,1)); 
else
    go.bounce_level = str2double(go.bounce_level);
end

%--------------------------------------------------------------------------

% Save image limits globally for tomographic analysis
az_min = go.a_min; % Minimum in Azimuth [m]
az_max = go.a_max; % Maximum in -''- [m]
ra_min = go.r_min; % Minimum in Slant-Range/Ground-Range [m]
ra_max = go.r_max; % Maximum in -''- [m]
bounce_rem = go.bounce_level; % Maximum bounce level to be analyzed [value between 1 and 5]

%--------------------------------------------------------------------------

% Intervals on axes
diff_a = go.a_max - go.a_min; % interval azimuth
diff_r = go.r_max - go.r_min; % interval slant range

num_a = ceil(diff_a/go.a_pix); % azimuth
num_r = ceil(diff_r/go.r_pix); % range

%--------------------------------------------------------------------------

% 2.) Predefine size of image matrices

%--------------------------------------------------------------------------

% Reflectivity Maps

Single = zeros(num_r,num_a); % For gathering single bounce contributions

if go.bounce_level > 1
    Double = zeros(num_r,num_a); % double ...
end

if go.bounce_level > 2
    Triple = zeros(num_r,num_a); % triple ...
end

if go.bounce_level > 3
    Fourfold = zeros(num_r,num_a); % fourfold ...
end

if go.bounce_level > 4
    Fivefold = zeros(num_r,num_a); % fivefold ...
end

%--------------------------------------------------------------------------

% Specular Maps

Specular_Single = zeros(num_r,num_a); % Only single bounce
Specular_All = zeros(num_r,num_a); % Single and multiple bounce

if go.bounce_level > 1
    Specular_Double = zeros(num_r,num_a); % Only double bounce
end

if go.bounce_level > 2
    Specular_Triple = zeros(num_r,num_a); % Only triple bounce
end

if go.bounce_level > 3
    Specular_Fourfold = zeros(num_r,num_a); % Only fourfold bounce
end

if go.bounce_level > 4
    Specular_Fivefold = zeros(num_r,num_a); % Only fivefold bounce
end

%--------------------------------------------------------------------------

% 3.) Fill values into image matrices

%--------------------------------------------------------------------------

% Loop over all contributions
for k = 1 : 1 : s
         
    % Assign signal contributions to image layers
    %----------------------------------------------------------------------
            
    % If contribution is located within defined margins
     if Az(k,1) >= go.a_min && Az(k,1) <= go.a_max
         if Ra_g(k,1) >= go.r_min && Ra_g(k,1) <= go.r_max
                                 
            %--------------------------------------------------------------
           
            % B.)Pixel center
            % --> rounded coordinate in azimuth + 0.5
            % --> rounded coordinate in range + 0.5

            %e.g. for first pixel: row_pix = 0.5, column_pix = 0.5
           
            row_pix = floor((Ra_g(k,1)-go.r_min)/go.r_pix)+0.5;
            column_pix = floor((Az(k,1)-go.a_min)/go.a_pix)+0.5;   
            
            %--------------------------------------------------------------

            % D.) Add contribution to pixels 

            % Round pixel coordinates for obtaining matrix indices
            column_pix = round(column_pix);
            row_pix = round(row_pix);
            
            if row_pix <= num_r && column_pix <= num_a
                                               
                % Preparation for coherent summation of signal contributions
                if go.coh == 1
                                        
                    % Calculate phase
                    Phi_temp = (-(4*pi)/lambda)*Ra(k,1); % ping pong
                    %Phi_temp = (-(2*pi)/lambda)*Ra(k,1); % no ping pong
                    Cycles = floor(Phi_temp/(2*pi));
                    Phi = Phi_temp - Cycles*2*pi;
                                        
                    % Complex signal contribution
                    Signal = Intens(k,1)*cos(Phi)+1i*Intens(k,1)*sin(Phi);
                                                        
                 % Preparation for uncoherent summation (absolute values only)
                else  
                    % Absolute
                    Signal = Intens(k,1);
                end
               
                % Single bounce
                if Tr_L(k,1) == 1
                    Single(row_pix, column_pix) = Single(row_pix, column_pix)+Signal; 

                    % Mark specular reflected contributions
                    if Sp(k,1) == 1
                        Specular_All(row_pix,column_pix) = 1;
                        Specular_Single(row_pix,column_pix) = 1; 
                    end
                end

                % Double bounce
                if Tr_L(k,1) == 2 && go.bounce_level > 1
                    Double(row_pix, column_pix) = Double(row_pix, column_pix)+Signal;

                    % Mark specular reflected contributions
                    if Sp(k,1) == 1
                        Specular_All(row_pix,column_pix) = 1;
                        Specular_Double(row_pix,column_pix) = 1; 
                    end
                end
                
                % Triple bounce
                if Tr_L(k,1) == 3 && go.bounce_level > 2 
                    
                    Triple(row_pix, column_pix) = Triple(row_pix, column_pix)+Signal;

                    % Mark specular reflected contributions
                    if Sp(k,1) == 1 
                        Specular_All(row_pix,column_pix) = 1;
                        Specular_Triple(row_pix,column_pix) = 1; 
                    end
                end

                % Fourfold bounce
                if Tr_L(k,1) == 4 && go.bounce_level > 3
                    Fourfold(row_pix, column_pix) = Fourfold(row_pix, column_pix)+Signal;

                   % Mark specular reflected contributions
                    if Sp(k,1) == 1     
                        Specular_All(row_pix,column_pix) = 1;
                        Specular_Fourfold(row_pix,column_pix) = 1; 
                    end
                end

                % Fivefold bounce
                if Tr_L(k,1) == 5 && go.bounce_level > 4
                   Fivefold(row_pix, column_pix) = Fivefold(row_pix, column_pix)+Signal;

                   % Mark specular reflected contributions
                    if Sp(k,1) == 1   
                        Specular_All(row_pix,column_pix) = 1;
                        Specular_Fivefold(row_pix,column_pix) = 1; 
                    end
                end
            end
         end
    end
end

%--------------------------------------------------------------------------

% 4.) Preparations before displaying images

%--------------------------------------------------------------------------

% Plots for displaying results

% Find scatterers within area of interest

% Azimuth
Gr_Az_1 = find(Az >= go.a_min);
Gr_Az_2 = find(Az <= go.a_max);
Gr_Az = intersect(Gr_Az_1,Gr_Az_2);
clear Gr_Az_1 Gr_Az_2

% Range
Gr_Ra_1 = find(Ra_g >= go.r_min);
Gr_Ra_2 = find(Ra_g <= go.r_max);
Gr_Ra = intersect(Gr_Ra_1,Gr_Ra_2);
clear Gr_Ra_1 Gr_Ra_2

% Intersection of Azimuth and Range
Gr_Az_Ra = intersect(Gr_Az,Gr_Ra);
clear Gr_Az Gr_Ra

% Groups: Bounce Level
A_temp = find(Tr_L == 1); % Indices of contributions having bounce level 1
B_temp = find(Tr_L == 2); % ...
C_temp = find(Tr_L == 3);
D_temp = find(Tr_L == 4);
E_temp = find(Tr_L == 5);

% Gather contributions inside area of interest
A = intersect(A_temp,Gr_Az_Ra);
B = intersect(B_temp,Gr_Az_Ra);
C = intersect(C_temp,Gr_Az_Ra);
D = intersect(D_temp,Gr_Az_Ra);
E = intersect(E_temp,Gr_Az_Ra);
clear A_temp B_temp C_temp D_temp E_temp

%--------------------------------------------------------------------------

if isempty(Gr_Az_Ra) == 1
    
    % Error
    msgbox('No data detected within selected area of interest. Please check limits of reflectivity map to be generated.','No data selected','error')
else

    % Summing up signals according to chosen bounce maximum
    % in case of imaginary numbers: maps are added coherently --> absolute
    % values are calculated afterwards

    All_reflections = Single;
        
    if go.bounce_level > 1 && size(B,1) > 0
        % Sum of signals
        All_reflections = All_reflections+Double;
    end

    if go.bounce_level > 2 && size(C,1) > 0
        % Sum of signals
        All_reflections = All_reflections+Triple;
    end

    if go.bounce_level > 3 && size(D,1) > 0
        % Sum of signals
        All_reflections = All_reflections+Fourfold;
    end

    if go.bounce_level > 4 && size(E,1) > 0
        % Sum of signals
        All_reflections = All_reflections+Fivefold;
    end
    
    %----------------------------------------------------------------------
    
    % Store complex image data to file
    % --> original data are stored before being manipulated (dB, clipping)
       
    if go.store == 1
                       
        mkdir([Output_path,'/Maps/Image_data']);
        %save ../Results/Maps/Image_data/Image_Matrix All_reflections 
        save([Output_path '/Maps/Image_data/Image_Matrix'], 'All_reflections') 
        
        msgbox('The image data has been stored to file "Image_Matrix.mat".','Image data','help');
    end
    
    %----------------------------------------------------------------------
    
    % In case of coherent summation: calculate absolute value 
    if go.coh == 1
               
        % optional: add phase noise
        if go.noise == 1
            All_reflections = add_phase_noise(All_reflections,go.phase_noise,num_r,num_a);
            Single = add_phase_noise(Single,go.phase_noise,num_r,num_a);
        end
        
        All_reflections = abs(All_reflections);
        Single = abs(Single);
    
        if go.bounce_level > 1
            
            % optional: add phase noise
            if go.noise == 1
                Double = add_phase_noise(Double,go.phase_noise,num_r,num_a);
            end
            
            Double = abs(Double); % double ...
        end

        if go.bounce_level > 2
            
            % optional: add phase noise
            if go.noise == 1
                Triple = add_phase_noise(Triple,go.phase_noise,num_r,num_a);
            end
            
            Triple = abs(Triple); % triple ...
        end

        if go.bounce_level > 3
            
            % optional: add phase noise
            if go.noise == 1
                Fourfold = add_phase_noise(Fourfold,go.phase_noise,num_r,num_a);
            end
            
            Fourfold = abs(Fourfold); % fourfold ...
        end

        if go.bounce_level > 4
            
            % optional: add phase noise
            if go.noise == 1
                Fivefold = add_phase_noise(Fivefold,go.phase_noise,num_r,num_a);
            end
            
            Fivefold = abs(Fivefold); % fivefold ...
        end
    end
        
    %--------------------------------------------------------------------------

    % Optional: Smoothing

    if go.filt == 1

        % Smooth images by 3x3 Binomial-Filter

        Single = Filter_Ref(Single);

        if exist('Double','var') == 1
            Double = Filter_Ref(Double);
        end

        if exist('Triple','var') == 1
            Triple = Filter_Ref(Triple);
        end

        if exist('Fourfold','var') == 1
            Fourfold = Filter_Ref(Fourfold);
        end

        if exist('Fivefold','var') == 1
            Fivefold = Filter_Ref(Fivefold);
        end

        if exist('All_reflections','var') == 1
            All_reflections = Filter_Ref(All_reflections);
        end

    end
    
    %----------------------------------------------------------------------
    
    % Optional: Add system response function
           
    if go.sinc_sim == 1
        
        % Single Bounce
        if exist('Single','var') == 1
            Single = add_impulse_resp(Single,num_a,num_r,go.a_pix,go.r_pix,go.a_res,go.r_res);
        end  

        % Double Bounce
        if exist('Double','var') == 1
            Double = add_impulse_resp(Double,num_a,num_r,go.a_pix,go.r_pix,go.a_res,go.r_res);
        end  

        % Triple Bounce
        if exist('Triple','var') == 1
            Triple = add_impulse_resp(Triple,num_a,num_r,go.a_pix,go.r_pix,go.a_res,go.r_res);
        end  

        % Fourfold Bounce
        if exist('Fourfold','var') == 1
            Fourfold = add_impulse_resp(Fourfold,num_a,num_r,go.a_pix,go.r_pix,go.a_res,go.r_res);
        end  

        % Fivefold Bounce
        if exist('Fivefold','var') == 1
            Fivefold = add_impulse_resp(Fivefold,num_a,num_r,go.a_pix,go.r_pix,go.a_res,go.r_res);
        end  
                        
        % Reflectivity Map
        if exist('All_reflections','var') == 1
            All_reflections = add_impulse_resp(All_reflections,num_a,num_r,go.a_pix,go.r_pix,go.a_res,go.r_res);
        end        
    end
        
    %----------------------------------------------------------------------
     
    % Logarithmic scaling or clipping
    
    if go.dBnodB == 1
        
        % Logarithmic scaling

        Single = 10*log(Single); 
        All_reflections = 10*log(All_reflections);

        if exist('Double','var') == 1
            Double = 10*log(Double); 
        end
        if exist('Triple','var') == 1
            Triple = 10*log(Triple); 
        end
        if exist('Fourfold','var') == 1 
            Fourfold = 10*log(Fourfold); 
        end
        if exist('Fivefold','var') == 1
            Fivefold = 10*log(Fivefold); 
        end

        % Detect Maximum Intensity to be displayed

        % Minimum value [dB]
        if strcmp(go.db_min,'Min') == 1 % minimum value to be taken
            im_min = -100;
        else % numerical value chosen by operator
            im_min = str2double(go.db_min);
        end

        % Maximum value [dB]
        if strcmp(go.db_max,'Max') == 1 % maximum value to be taken
            if go.bounce_level > 1
                im_max = max(max(All_reflections));
            else
                im_max = max(max(Single));
            end
        else % numberical value chosen by operator
             im_max = str2double(go.db_max);
        end

        % Necessary correction for logarithmically scaled matrices
        % Set minimum within image to db_min
        % --> remove entries = -Inf
        All_reflections(All_reflections < im_min) = im_min;
        All_reflections(All_reflections > im_max) = im_max;
        
        % Show histogram of intensities
        figure;
        hist_min = min(min(All_reflections));
        hist_max = max(max(All_reflections));
        hist_diff = hist_max-hist_min;
        hist_step = hist_min:hist_diff/25:hist_max;
        hist(All_reflections(:),hist_step); title('Distribution of intensities (dB)');
        xlabel('Intensity'); ylabel('Occurrence');
        set(gcf,'Name','Histogram of Intensities (dB)','Numbertitle','off');
        grid;   

    else
        
        % Clipping yes/no
        if strcmp(go.clip,'Max') ~= 1

            % Maximum
            max_clip = str2double(go.clip);
            
            % Avoid threshold bigger than given intensity maximum
            if max_clip <= max(max(All_reflections)) && max_clip >= 0

                % Perform clipping

                Single(Single > max_clip) =  max_clip;
                All_reflections(All_reflections > max_clip) = max_clip;

                if exist('Double','var') == 1
                    Double(Double > max_clip) =  max_clip; 
                end
                if exist('Triple','var') == 1
                    Triple(Triple > max_clip) =  max_clip; 
                end
                if exist('Fourfold','var') == 1 
                    Fourfold(Fourfold > max_clip) =  max_clip;
                end
                if exist('Fivefold','var') == 1
                    Fivefold(Fivefold > max_clip) =  max_clip;
                end
                
                % Clipping applied
                clipping_status = 1;
            else
                % No clipping applied
                clipping_status = 0;
            end
        end

        % Extrema intensities
        im_min = min(min(All_reflections));
        im_max = max(max(All_reflections));

        % Show histogram of amplitudes
        figure;
        hist_min = min(min(All_reflections));
        hist_max = max(max(All_reflections));
        hist_diff = hist_max-hist_min;
        hist_step = hist_min:hist_diff/25:hist_max;
        hist(All_reflections(:),hist_step); title('Distribution of intensities (clipping)');
        xlabel('Intensity'); ylabel('Occurrence');
        set(gcf,'Name','Histogram of Intensities (clipping)','Numbertitle','off');
        grid;
    end
    
    % Size of images to be displayed
    r_p = 1:1:size(Single,1); % rows --> range
    a_p = 1:1:size(Single,2); % columns --> azimuth
          
    %--------------------------------------------------------------------------
    
    % 5.) Plots

    %--------------------------------------------------------------------------

    % A.) Geometrical distribution of reflection contributions

    % Display scatterers in azimuth-range plane (2D)

    % Index variable
    index = 1;

    figure;
    set(gcf,'Name','Phase Centers (2D)','Numbertitle','off')
    hold on;

    if size(A,1) > 0
        P0 = plot(Az(A(1:5:length(A))), Ra_g(A(1:5:length(A))),'.b');
        set(P0,'markersize',12);
        L1{index,1} = 'single bounce';
        index = index+1;
    end

    if go.bounce_level > 1 && size(B,1) > 0
        P0 = plot(Az(B (1:3:length(B))), Ra_g(B(1:3:length(B))),'.g');
        set(P0,'markersize',12);
        L1{index,1} = 'double bounce';
        index = index+1;
    end

    if go.bounce_level > 2 && size(C,1) > 0
        P0 = plot(Az(C), Ra_g(C),'.r');
        set(P0,'markersize',12);
        L1{index,1} = 'triple bounce';
        index = index+1;
    end

    if go.bounce_level > 3 && size(D,1) > 0
        P0 = plot(Az(D), Ra_g(D),'.m');
        set(P0,'markersize',12);
        L1{index,1} = 'fourfold bounce';
        index = index+1;
    end

    if go.bounce_level > 4 && size(E,1) > 0
        P0 = plot(Az(E), Ra_g(E),'.c');
        set(P0,'markersize',12);
        L1{index,1} = 'fivefold bounce';
    end
    
    % display in ground range 
    if r_geom == 1
        xlabel('Azimuth','fontsize',14); ylabel('Ground Range','fontsize',14); title('Distribution of Intensities','fontsize',14);
    end

    % display in slant range
    if r_geom == 0
        xlabel('Azimuth','fontsize',14); ylabel('Slant Range','fontsize',14); title('Distribution of Intensities','fontsize',14);
    end
       
    if range_dir == 1
                
        % Change label of axes
        % --> required due to rotation of figure content
        set(gca,'YDir','reverse');
        set(gca,'XDir','reverse');
    end

    axis equal; axis image;
    
    % Adapt image extent to reflectivity map
    axis([go.a_min go.a_max go.r_min go.r_max]);

    % Legend
    legend(L1,'Location','NorthEastOutside');    
    
    hold off;
    
    set(gca,'FontSize',14);
    
    %--------------------------------------------------
    % Save image to file
    
    % Create new folder
    if exist([Output_path,'/Maps'],'dir') ~= 7
        mkdir([Output_path,'/Maps']);
    end
    
    % Create subfolders
    mkdir([Output_path,'/Maps/Figures']);
    mkdir([Output_path,'/Maps/Figures/Ref_Maps']);
    mkdir([Output_path,'/Maps/Figures/Ref_Maps/JPG']);
    mkdir([Output_path,'/Maps/Figures/Ref_Maps/FIG']);
    mkdir([Output_path,'/Maps/Frames']);
    mkdir([Output_path,'/Maps/Frames/Ref_Maps']);
    
    % Store image in folder
    mkdir([Output_path,'/Maps/Signal_Distribution']);
    saveas(P0,[Output_path,'/Maps/Signal_Distribution/Phase_Centers_2D.png'],'png');
    saveas(P0,[Output_path,'/Maps/Signal_Distribution/Phase_Centers_2D.fig'],'fig');
    
    %--------------------------------------------------------------------------

    % B.) Reflectivity maps

    %--------------------------------------------------------------------------
    
    % Single bounce
    if size(A,1) > 0
        
        Bounce = 'Single Bounce';
        R_Map(Single,im_min,im_max,Bounce,a_p,r_p,go.dBnodB);
    end

    %--------------------------------------------------------------------------

    % Double bounce
    if go.bounce_level > 1 && size(B,1) > 0

        Bounce = 'Double Bounce';
        R_Map(Double,im_min,im_max,Bounce,a_p,r_p,go.dBnodB);  
    end

    %--------------------------------------------------------------------------

    % Triple bounce
    if go.bounce_level > 2 && size(C,1) > 0

        Bounce = 'Triple Bounce';
        R_Map(Triple,im_min,im_max,Bounce,a_p,r_p,go.dBnodB);   
    end

    %--------------------------------------------------------------------------

    % Fourfold bounce
    if go.bounce_level > 3 && size(D,1) > 0

        Bounce = 'Fourfold Bounce';
        R_Map(Fourfold,im_min,im_max,Bounce,a_p,r_p,go.dBnodB);    
    end

    %--------------------------------------------------------------------------

    % Fivefold bounce
    if go.bounce_level > 4 && size(E,1) > 0

        Bounce = 'Fivefold Bounce';
        R_Map(Fivefold,im_min,im_max,Bounce,a_p,r_p,go.dBnodB);
    end

    %--------------------------------------------------------------------------

    % Display all bounces in one plot
    if go.bounce_level > 1

        Bounce = 'All Reflections';
        R_Map(All_reflections,im_min,im_max,Bounce,a_p,r_p,go.dBnodB);
        Index_All = gcf;
    end
            
    % If specular images have to be provided
    if go.specular_flag == 1
        
        % Create subfolders
        mkdir([Output_path,'/Maps/Figures/Specular']);
        mkdir([Output_path,'/Maps/Figures/Specular/JPG']);
        mkdir([Output_path,'/Maps/Figures/Specular/FIG']);
        mkdir([Output_path,'/Maps/Frames/Specular']);
                
        %--------------------------------------------------------------------------

        % C.) Specular maps

        %--------------------------------------------------------------------------

        % Single bounce

        if go.bounce_level > 0 && size(A,1) > 0    
            if max(max(Specular_Single)) == 1 % if specular components detected

                Bounce = 'Specular Analysis - Single Bounce';
                S_Map(Specular_Single,Bounce,a_p,r_p);

            end
        end

        % Double bounce

        if go.bounce_level > 1 && size(B,1) > 0    
            if max(max(Specular_Double)) == 1 % if specular components detected

                Bounce = 'Specular Analysis - Double Bounce';
                S_Map(Specular_Double,Bounce,a_p,r_p);

            end
        end
        %--------------------------------------------------------------------------

        % Triple bounce
        
        if go.bounce_level > 2 && size(C,1) > 0    
            if max(max(Specular_Triple)) == 1 % if specular components detected

                Bounce = 'Specular Analysis - Triple Bounce';
                S_Map(Specular_Triple,Bounce,a_p,r_p);    

            end
        end
        %--------------------------------------------------------------------------

        % Fourfold bounce

        if go.bounce_level > 3 && size(D,1) > 0 
            if max(max(Specular_Fourfold)) == 1 % if specular components detected

                Bounce = 'Specular Analysis - Fourfold Bounce';
                S_Map(Specular_Fourfold,Bounce,a_p,r_p);  

            end
        end

        %--------------------------------------------------------------------------

        % Fivefold bounce

        if go.bounce_level > 4 && size(E,1) > 0  
            if max(max(Specular_Fivefold)) == 1 % if specular components detected

                Bounce = 'Specular Analysis - Fivefold Bounce';
                S_Map(Specular_Fivefold,Bounce,a_p,r_p);

            end
        end
        %--------------------------------------------------------------------------

        % All specular bounces
        if sum(sum(Specular_All)) ~= 0

            Bounce = 'Specular Analysis - Multiple Bounce';
            S_Map(Specular_All,Bounce,a_p,r_p);

        end
    end
    
    %--------------------------------------------------------------------------
    
    % Warn message
    
    % Case: clipping failed
    if exist('clipping_status','var') == 1
        
        if clipping_status == 0
            msgbox('No clipping applied due to inappropriate threshold selection.','Clipping failed');
        end
    end
end    
    
