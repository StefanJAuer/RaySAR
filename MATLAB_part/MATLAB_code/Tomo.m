function Tomo(a_min,a_max,r_min,r_max,e_samp,min_el,max_el,angle)

% function file for displaying tomographic information

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input:
% - interval in azimuth [m] (a_min ---> a_max; width of profile poiting in range)
% - interval in range [m] (r_min --> r_max; width of profile poiting in azimuth)
% - angle: angle of incidence (default: -1 [deg])
% - e_samp: sampling stepwidth in elevation [m]
% - min_el, max_el: limits of signal samples of interest in elevation [m]

% output: at most three plots
% - distribution of signal response in elevation direction
% - profile in azimuth direction
% - profile in range direction

%--------------------------------------------------------------------------

% global parameters

global Az Ra El Intens Tr_L check_az check_ra check_el check_geom Output_path;

% Az: azimuth coordinates of reflection contributions [m]
% Ra: range coordinates of -''- [m]
% El: elevation coordinates of -''- [m]
% Intens: intensity of -''- [value between 0 and 1]
% Tr_L: bounce level of -''- [value between 1 and 5]
% check_az: display slice along slant range axis at fixed azimuth interval? no --> value 0, yes --> value 1
% check_ra: display slice along azimuth axis at fixed slant range interval? no --> value 0, yes --> value 1
% check_el: display slice along elevation axis? no --> value 0, yes --> value 1
% check_geom: way of displaying height information; 
%             value 0 --> elevation with respect to center of synthetic aperture in elevation direction [m]
%             value 1 --> height over ground [m]
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------
% slant or ground range

% rotation angle if slices have to be displayed in "height over ground"
ang_comp = (pi/2-angle)*(-1);

%--------------------------------------------------------------------------

% Detect contributions within azimuth and range intervals

% start at minimum in elevation
    
% Fix minimum in elevation
if strcmp(min_el,'Min')
    min_el = min(El);
else
    min_el = str2double(min_el);
end

% Fix maximum in elevation
if strcmp(max_el,'Max')
    max_el = max(El);
else
    max_el = str2double(max_el);
end
            
% azimuth
Ind_az_p1 = find(Az >= a_min);
Ind_az_p2 = find(Az <= a_max);
Ind_az = intersect(Ind_az_p1,Ind_az_p2); 

% range
Ind_ra_p1 = find(Ra >= r_min);
Ind_ra_p2 = find(Ra <= r_max);
Ind_ra = intersect(Ind_ra_p1,Ind_ra_p2);

% elevation
Ind_el_p1 = find(El >= min_el); 
Ind_el_p2 = find(El <= max_el); 
Ind_el = intersect(Ind_el_p1,Ind_el_p2);
    
%--------------------------------------------------------------------------

% Display amplitude distribution along the elevation direction
if check_el == 1 % check if plot has to be displayed
     
    %---------------------------------
    
    % Result: slice in elevation direction
    
    % axis left to right: elevation [m]
    % axis bottom up: intensity
        
    % Detect contributions within azimuth-range cell
            
    % Intersection
    Ind_res = intersect(Ind_az,Ind_ra); % Intersect
        
    %---------------------------------
    
    % Create folder for storing the plot
    if exist([Output_path,'/Maps/Figures/Elevation_Slice'],'dir') ~= 7
        mkdir([Output_path,'/Maps/Figures/Elevation_Slice']);  
    end
        
    % Separate different bounce levels 
    Ind_si = Ind_res(Tr_L(Ind_res) == 1); % single bounce
    Ind_do = Ind_res(Tr_L(Ind_res) == 2); % double
    Ind_tr = Ind_res(Tr_L(Ind_res) == 3); % triple
    Ind_fo = Ind_res(Tr_L(Ind_res) == 4); % fourfold
    Ind_fi = Ind_res(Tr_L(Ind_res) == 5); % fivefold
    
    %---------------------------------      
    
    % test variable
    check = 0;
    
    % Transform elevation data to zero height plane if necessary
    if check_geom == 0  
        % Display results in Slant-Range - Elevation plane
        El_f = El;
        bottom_El = 0;
    else
            
        % detect minimum elevation in subset
        bottom_El = min(El(Ind_si));
        
        if isempty(bottom_El) == 1
            msgbox('Transformation in height over ground failed. No ground information available in shadow regions. A workaround still has to be provided :-).','Transformation impossible','warning');
            
            % test variable
            check = 1;
        else
            index_bot = find(El(:,1) == bottom_El);
        
            % Display results in height over ground
            El_f = sin(ang_comp).*Ra + cos(ang_comp).*El;       
            bottom_El = sin(ang_comp)*Ra(index_bot(1),1) + cos(ang_comp)*El(index_bot(1),1);
        end
    end
    
    if check == 0
    
        % Shift to the ground
        El_f = El_f-bottom_El;

        % detect minimum elevation in subset
        min_el_f = min(El_f(Ind_si));
        max_el_f = max(El_f(Ind_si));

        % spanwidth in elevation
        span = max_el_f - min_el_f;

        % number of steps
        num_el = ceil(span/e_samp);

        % Create vectors

        % a.) create elevation vector
        el_vec = min_el_f : e_samp : min_el_f+e_samp*num_el;

        % b.) bounce layers
        el_single = zeros(1,num_el);
        el_double = zeros(1,num_el);
        el_triple = zeros(1,num_el);
        el_fourfold = zeros(1,num_el);
        el_fivefold = zeros(1,num_el);

        % start loop --> move along elevation axis (stepwidth: e_samp)

        for i = 1:num_el+1

            % limits for step in elevation direction
            el_low = min_el_f+(i-1)*e_samp-e_samp/2;
            el_high = min_el_f+i*e_samp-e_samp/2;

            %----------------------------------------------
            % detect intensity contributions

            % single bounce
            Ind_part1 = find(El_f(Ind_si) >= el_low);
            Ind_part2 = find(El_f(Ind_si) < el_high);
            Ind = intersect(Ind_part1,Ind_part2);

            % add to single bounce map
            el_single(i) = sum(Intens(Ind_si(Ind)));

            clear Ind_part1 Ind_part2 Ind;

            %--------------------------------------

            % double bounce
            Ind_part1 = find(El_f(Ind_do) >= el_low);
            Ind_part2 = find(El_f(Ind_do) < el_high);
            Ind = intersect(Ind_part1,Ind_part2);

            % add to double bounce map
            el_double(i) = sum(Intens(Ind_do(Ind)));       

            clear Ind_part1 Ind_part2 Ind;

            %--------------------------------------

            % triple bounce
            Ind_part1 = find(El_f(Ind_tr) >= el_low);
            Ind_part2 = find(El_f(Ind_tr) < el_high);
            Ind = intersect(Ind_part1,Ind_part2);

            % add to triple bounce map
            el_triple(i) = sum(Intens(Ind_tr(Ind)));

            clear Ind_part1 Ind_part2 Ind;
            %--------------------------------------

            % fourfold bounce
            Ind_part1 = find(El_f(Ind_fo) >= el_low);
            Ind_part2 = find(El_f(Ind_fo) < el_high);
            Ind = intersect(Ind_part1,Ind_part2);

            % add to fourfold bounce map
            el_fourfold(i) = sum(Intens(Ind_fo(Ind)));

            clear Ind_part1 Ind_part2 Ind;
            %--------------------------------------

            % fivefold bounce
            Ind_part1 = find(El_f(Ind_fi) >= el_low);
            Ind_part2 = find(El_f(Ind_fi) < el_high);
            Ind = intersect(Ind_part1,Ind_part2);

            % add to fivefold bounce map
            el_fivefold(i) = sum(Intens(Ind_fi(Ind)));

            clear Ind_part1 Ind_part2 Ind;
        end

        % detect maximum amplitude
        max_int = 0;

        if max(el_single) > max_int
            max_int = max(el_single);
        end

        if max(el_double) > max_int
            max_int = max(el_double);
        end

        if max(el_triple) > max_int
            max_int = max(el_triple);
        end

        if max(el_fourfold) > max_int
            max_int = max(el_fourfold);
        end

        if max(el_fivefold) > max_int
            max_int = max(el_fivefold);
        end

        %--------------------------------------------------------------
        % display results

        % Single bounce

        figure; hold on;
        index = 1;

        % width of bars: adjacent bars touch each other --> value 1
        bar_width = 1;

        if sum(el_single) ~= 0

            % Single bounce
            Elev_Sl = bar(el_vec,el_single./max_int,bar_width,'b'); hold on;
            L1{index} = 'Single Bounce';
            index = index+1;
        end

        if sum(el_double) ~= 0

            % Double bounce
            Elev_Sl = bar(el_vec,el_double./max_int,bar_width,'g'); 
            L1{index} = 'Double Bounce';
            index = index+1;
        end    

        if sum(el_triple) ~= 0

            % Triple bounce
            Elev_Sl = bar(el_vec,el_triple./max_int,bar_width,'r');
            L1{index} = 'Triple Bounce';
            index = index+1;
        end

        if sum(el_fourfold) ~= 0

            % Fourfold bounce
            Elev_Sl = bar(el_vec,el_fourfold./max_int,bar_width,'m');
            L1{index} = 'Fourfold Bounce';
            index = index+1;
        end

        if sum(el_fivefold) ~= 0

            % Fivefold bounce
            Elev_Sl = bar(el_vec,el_fivefold./max_int,bar_width,'c'); 
            L1{index} = 'Fivefold Bounce';
        end


        % Plot labels
        if check_geom == 0
            xlabel('Elevation [m]');
        else
            xlabel('Height over ground [m]');
        end

        ylabel('Amplitude (Normalized)'); title('Reflection Contributions');
        set(gcf,'Name','Elevation - Single + Multiple Bounce','Numbertitle','off');  
        legend(L1,'Location','NorthEastOutside'); grid;
        hold off;

        % Store image in folder
        % Store plot in folder
        if check_geom == 0
            saveas(Elev_Sl,[Output_path,'/Maps/Figures/Elevation_Slice/Elevation_Slice_El.png'],'png');
            saveas(Elev_Sl,[Output_path,'/Maps/Figures/Elevation_Slice/Elevation_Slice_El.fig'],'fig');
        else
            saveas(Elev_Sl,[Output_path,'/Maps/Figures/Elevation_Slice/Elevation_Slice_HOG.png'],'png');
            saveas(Elev_Sl,[Output_path,'/Maps/Figures/Elevation_Slice/Elevation_Slice_HOG.fig'],'fig');
        end
    end
end

% Display elevation information along selected range interval
if check_ra == 1 % check if plot has to be displayed
    
    %---------------------------------
    
    % Result: slice in azimuth direction
    
    % axis left to right: azimuth [m]
    % axis bottom up: elevation [m]
    
    %---------------------------------
    
    % Create folder for storing the plot
    if exist([Output_path,'/Maps/Figures/Azimuth_Slice'],'dir') ~= 7
        mkdir([Output_path,'/Maps/Figures/Azimuth_Slice']);  
    end
    
    % All contributions in azimuth direction fulfilling defined elevation thresholds
    Ind_ra2 = intersect(Ind_ra,Ind_el);
        
    % Separate different bounce levels
    Ind_si2 = Ind_ra2(Tr_L(Ind_ra2) == 1); % single bounce
    Ind_do2 = Ind_ra2(Tr_L(Ind_ra2) == 2); % double
    Ind_tr2 = Ind_ra2(Tr_L(Ind_ra2) == 3); % triple
    Ind_fo2 = Ind_ra2(Tr_L(Ind_ra2) == 4); % fourfold
    Ind_fi2 = Ind_ra2(Tr_L(Ind_ra2) == 5); % fivefold
    
    % Transform elevation data to zero height plane if required
    if check_geom == 0 
        El_f = El;
        bottom_El_1 = 0;
    else
        % Display results in Ground - Elevation plane
        El_f = sin(ang_comp).*Ra + cos(ang_comp).*El;
        bottom_El_1 = min(El_f(Ind_si2));
    end
    
    % Shift to ground level
    if bottom_El_1 < 0
        El_f = El_f-bottom_El_1;
    else
        El_f = El_f+bottom_El_1;
    end
        
    % display result
    figure;
    
    % index for legend entries
    index_1 = 1;
    
    % Single Bounce
    P1 = plot(Az(Ind_si2),El_f(Ind_si2),'.b'); hold on; grid; axis equal; 
  
    % Store information for legend (to be displayed in plot)
    L2{index_1} = 'Single Dounce';
    index_1 = index_1 + 1;
    
    if isempty(Ind_do2) == 0
        % Double Bounce
        P1 = plot(Az(Ind_do2),El_f(Ind_do2),'og'); 
        set(P1,'MarkerSize',8);
        set(P1,'MarkerFaceColor','g');
        L2{index_1} = 'Double Bounce';
        index_1 = index_1 + 1;         
    end
    
    if isempty(Ind_tr2) == 0
        % Triple Bounce
        P1 = plot(Az(Ind_tr2),El_f(Ind_tr2),'or'); 
        set(P1,'MarkerSize',8);
        set(P1,'MarkerFaceColor','r');
        L2{index_1} = 'Triple Bounce';
        index_1 = index_1 + 1;
    end
    
    if isempty(Ind_fo2) == 0
        % Fourfold Bounce
        P1 = plot(Az(Ind_fo2),El_f(Ind_fo2),'om'); 
        set(P1,'MarkerSize',8);
        set(P1,'MarkerFaceColor','m');
        L2{index_1} = 'Fourfold Bounce';
        index_1 = index_1 + 1;
    end
    
    if isempty(Ind_fi2) == 0
        % Fivefold Bounce
        P1 = plot(Az(Ind_fi2),El_f(Ind_fi2),'oc'); 
        set(P1,'MarkerSize',8);
        set(P1,'MarkerFaceColor','c');
        L2{index_1} = 'Fivefold Bounce';
    end
        
    % Label for y-axis
    if check_geom == 0
        ylabel('Elevation [m]');
    else  
        ylabel('Height over ground [m]');
    end
    
    title(['Slice in azimuth direction: slant range interval between ' num2str(r_min,'%10.2f') ' m and ' num2str(r_max,'%10.2f') ' m']);
    xlabel('Azimuth [m]'); 
    set(gcf,'Name','Azimuth direction','Numbertitle','off');
    legend(L2,'Location','NorthEastOutside');
    hold off;   
    
    % Store plot in folder
    if check_geom == 0
        saveas(P1,[Output_path,'/Maps/Figures/Azimuth_Slice/Azimuth_Slice_EL.png'],'png');
        saveas(P1,[Output_path,'/Maps/Figures/Azimuth_Slice/Azimuth_Slice_EL.fig'],'fig');
    else
        saveas(P1,[Output_path,'/Maps/Figures/Azimuth_Slice/Azimuth_Slice_HOG.png'],'png');
        saveas(P1,[Output_path,'/Maps/Figures/Azimuth_Slice/Azimuth_Slice_HOG.fig'],'fig');
    end
end

% Display elevation information along selected azimuth interval
if check_az == 1 % check if plot has to be displayed
    
    %---------------------------------
    
    % Result: Slice in range direction
    
    % axis left to right: range [m]
    % axis bottom up: elevation [m]
    
    %---------------------------------
    
    % Create folder for storing the plot
    if exist([Output_path,'/Maps/Figures/Range_Slice'],'dir') ~= 7
        mkdir([Output_path,'/Maps/Figures/Range_Slice']);  
    end
    
     % All contributions in range direction fulfilling defined elevationthresholds
    Ind_az2 = intersect(Ind_az,Ind_el);
        
    % Separate different bounce levels
    Ind_si3 = Ind_az2(Tr_L(Ind_az2) == 1); % single bounce
    Ind_do3 = Ind_az2(Tr_L(Ind_az2) == 2); % double
    Ind_tr3 = Ind_az2(Tr_L(Ind_az2) == 3); % triple
    Ind_fo3 = Ind_az2(Tr_L(Ind_az2) == 4); % fourfold
    Ind_fi3 = Ind_az2(Tr_L(Ind_az2) == 5); % fivefold
    
    % Provide data for display
    % --> check_geom = 0: x-axis: slant-range; y-axis: elevation with respect to center 
    %                     of synthetic aperture (elevation
    % --> check_geom = 1: x-axis: ground range; y-axis: elevation over ground
    
    % Apply elevation correction if needed   
    if check_geom == 0
        bottom_El_2 = 0;
    else
        % detect minimum elevation in subset
        bottom_El_2 = min(El(Ind_si3));
        index_bot = find(El(:,1) == bottom_El_2);
    end
       
    % Transform elevation data to zero height plane if required
    if check_geom == 0  
        % Display results in Slant-Range - Elevation plane
        Ra_f = Ra;
        El_f = El;
    else
        % Display results in Ground - Elevation plane
        Ra_f = cos(ang_comp).*Ra - sin(ang_comp).*El;
        El_f = sin(ang_comp).*Ra + cos(ang_comp).*El;
        bottom_El_2 = sin(ang_comp)*Ra(index_bot(1),1) + cos(ang_comp)*El(index_bot(1),1);
    end
      
    % display result
    figure;
    
    % index for legend entries
    index_2 = 1;
    
    % Single Bounce
    P2 = plot(Ra_f(Ind_si3),El_f(Ind_si3)-bottom_El_2,'.b'); hold on; axis equal; grid;
    
    % Store information for legend (to be displayed in plot)
    L3{index_2} = 'Single Bounce';
    index_2 = index_2 + 1;
    
    if isempty(Ind_do3) == 0
        % Double Bounce
        P2 = plot(Ra_f(Ind_do3),El_f(Ind_do3)-bottom_El_2,'og'); 
        set(P2,'MarkerSize',8);
        set(P2,'MarkerFaceColor','g');
        L3{index_2} = 'Double Bounce';
        index_2 = index_2 + 1;
    end
    
    if isempty(Ind_tr3) == 0
        % Triple Bounce
        P2 = plot(Ra_f(Ind_tr3),El_f(Ind_tr3)-bottom_El_2,'or'); 
        set(P2,'MarkerSize',8);
        set(P2,'MarkerFaceColor','r');
        L3{index_2} = 'Triple Bounce';
        index_2 = index_2 + 1;
    end
    
    if isempty(Ind_fo3) == 0
        % Fourfold Bounce
        P2 = plot(Ra_f(Ind_fo3),El_f(Ind_fo3)-bottom_El_2,'om'); 
        set(P2,'MarkerSize',8);
        set(P2,'MarkerFaceColor','m');
        L3{index_2} = 'Fourfold Bounce';
        index_2 = index_2 + 1;
    end
    
    if isempty(Ind_fi3) == 0
        % Fivefold Bounce
        P2 = plot(Ra_f(Ind_fi3),El_f(Ind_fi3)-bottom_El_2,'oc'); 
        set(P2,'MarkerSize',8);
        set(P2,'MarkerFaceColor','c');
        L3{index_2} = 'Fivefold Bounce';
    end
    
    % Label for y-axis
    if check_geom == 0
        ylabel('Elevation [m]'); 
        xlabel('Slant Range [m]');
    else
         ylabel('Height over ground [m]');
         xlabel('Ground distance [m]');
    end
    
    title(['Slice in slant range: azimuth interval between ' num2str(a_min,'%10.2f') ' m and ' num2str(a_max,'%10.2f') ' m']);
    
    set(gcf,'Name','Range direction','Numbertitle','off');
    legend(L3,'Location','NorthEastOutside');    
    hold off;
    
    % Store plot in folder
    if check_geom == 0
        saveas(P2,[Output_path,'/Maps/Figures/Range_Slice/Range_Slice_EL.png'],'png');
        saveas(P2,[Output_path,'/Maps/Figures/Range_Slice/Range_Slice_EL.fig'],'fig');
    else
        saveas(P2,[Output_path,'/Maps/Figures/Range_Slice/Range_Slice_HOG.png'],'png');
        saveas(P2,[Output_path,'/Maps/Figures/Range_Slice/Range_Slice_HOG.fig'],'fig');
    end
end