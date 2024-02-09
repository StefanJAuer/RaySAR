function Write_Intersect_Pts(a_min,a_max,r_min,r_max,c_size,trans_POV)

% function file for writing 3D intersection points to a CAD file

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input:
% - azimuth limits of selected pixel (a_min ---> a_max)
% - range limits of selected pixel (r_min --> r_max)
% - c_size: size of cubes to be modeled [m]
% - trans_POV: string containing POV-Ray transformation/scaling information to be
%              compensated

% output: 
% - geometry file representing the 3D distribution of intersection points
% corresponding to the seleced pixel (Format: Wavefront, .obj).

%--------------------------------------------------------------------------

% global parameters

global Az Ra Intens Tr_L S_X S_Y S_Z Output_path;

% Az: azimuth coordinates of reflection contributions [m]
% Ra: range coordinates of -''- [m]
% Tr_L: bounce level of -''- [value between 1 and 5]
% S_X, S_Y, S_Z: coordinates of intersection points detected at surfaces in 3D model scene [m]
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% 1.) Select data of interest

% Check whether intersection points have been loaded
if size(S_X,1) == 0
    
    % Error
    msgbox('No intersection points available. Please check input data.','Lack of Intersection Points','Error');
    
else

    % azimuth
    a_ROI_1 = find(Az >= a_min);
    a_ROI_2 = find(Az <= a_max);

    a_ROI = intersect(a_ROI_1,a_ROI_2);
    clear a_ROI_1 a_ROI_2

    % range
    r_ROI_1 = find(Ra >= r_min);
    r_ROI_2 = find(Ra <= r_max);

    r_ROI = intersect(r_ROI_1,r_ROI_2);
    clear r_ROI_1 r_ROI_2
    
    % intersect azimuth and range
    p_ROI = intersect(a_ROI,r_ROI); % contributions corresponding to selected pixel
    
    % only accept intersection points "visible" to the SAR 
    % --> zero amplitude at hidden points is not accepted
    I_ROI_temp = find(Intens > 0);
    p_ROI = intersect(p_ROI,I_ROI_temp);
                   
    % matrix marking intersection points of interest   
    tr_ROI = Tr_L(p_ROI); % vector containing bounce level information
    clear a_ROI r_ROI
                                  
    % Check whether at least one intersection point is included in p_ROI
    if size(p_ROI,1) == 0
        
        % Warning
        msgbox('No intersection points found for selected pixel.','Lack of Intersection Points','Warn');
        
    else
        
        %--------------------------------------------------------------------------

        % Get further intersection points in case of multiple reflections
        % --> so far: in case of double bounce only the index of the first intersection
        % point is included in p_ROI
        
        bounce_max = max(max(Tr_L(p_ROI)));
                          
        % find missing intersection points starting at bounce level 2
        for j = 2:bounce_max
         
            % find entries marking bounce level of interest
            ind_p = find(Tr_L(p_ROI) == j);
                                                                                  
            % Store indices of lower bounces correponding to bounce level of interest
            for k = 1:(j-1)
                                
               % number of elements to be added
               num_el = length(ind_p);
                                  
               % Select
               p_ROI_temp = zeros(num_el,1); 
               tr_ROI_temp = zeros(num_el,1);
               p_ROI_temp(:,1) = p_ROI(ind_p,1)-k; % get intersection point of former bounce
               tr_ROI_temp(:,1) = j; % assign bounce level of current bounce maximum
                                
               % add indices to vector p_ROI
               p_ROI = [p_ROI;p_ROI_temp];
               % add bounce level information to tr_ROI
               tr_ROI = [tr_ROI;tr_ROI_temp];
            
            end            
        end
                              
        %--------------------------------------------------------------------------
        % Define final matrix containing points and bounce level information
        
        M_ROI = zeros(length(p_ROI),4);
        M_ROI(:,1) = S_X(p_ROI(:,1),1);
        M_ROI(:,2) = S_Y(p_ROI(:,1),1);
        M_ROI(:,3) = S_Z(p_ROI(:,1),1);
        M_ROI(:,4) = tr_ROI;
                
        % Remove identical points
        M_ROI = remove_red_inter(M_ROI);
                                                
        %--------------------------------------------------------------------------

        % 2.) Remove rotations, translations and scaling performed in POV Ray Editor

        % set flag to 1
        flag = 1;

        % Create temporay matrix for point transformations
        M_ROI_temp = M_ROI(:,1:3);

        if isempty(trans_POV) == 0

            % Get transformation parameters for compensating POV Ray transformations
            [M_trans,flag] = Check_Trans(trans_POV,flag);
            
            % Conduct necessary transformations 
            for i = 1:size(M_trans,1)

                if M_trans(i,4) == 0

                    % Rotation
                    rot_ang = M_trans(i,1:3)*(pi/180)*(-1); % [rad]
                    M_ROI_temp = rotation_M(M_ROI_temp,rot_ang); 
                end

                if M_trans(i,4) == 1

                    % Translation
                    transf = M_trans(i,1:3); % [m]
                    M_ROI_temp = translate_M(M_ROI_temp,transf);  
                end

                if M_trans(i,4) == 2

                    % Scaling
                    sc_f = M_trans(i,1:3); % [dimensionless]
                    M_ROI_temp = scale_M(M_ROI_temp,sc_f); 
                end
            end
        end

        % Correct for change of coordinate system: z --> -z
        % comment: axis is mirrowed in POV Ray coordinate system (right system -->
        % left system)
        M_ROI_temp(:,3) = M_ROI_temp(:,3)*(-1);
                
        % Update coordinates in matrix M_ROI
        M_ROI(:,1:3) = M_ROI_temp;
                
        %--------------------------------------------------------------------------

        % 3.) Write data to geometry file

        if flag == 1

            if exist([Output_path,'/Maps/Intersect'],'dir') ~= 7
               mkdir([Output_path,'/Maps/Intersect']); 
            end

            if exist([Output_path,'/Maps/Intersect/Inters_Pts.obj'],'file') == 2
               delete([Output_path,'/Maps/Intersect/Inters_Pts.obj']); 
            end

            if exist([Output_path,'/Maps/Intersect/Model_Color.mtl'],'file') == 2
                delete([Output_path,'/Maps/Intersect/Model_Color.mtl']);
            end

            % Create new or recreate
            fid = fopen([Output_path,'/Maps/Intersect/Inters_Pts.obj'],'w');
            fid2 = fopen([Output_path,'/Maps/Intersect/Model_Color.mtl'],'w');

            % header
            fprintf(fid,'mtllib Model_Color.mtl\n');

            % Create material file
            fprintf(fid2,'newmtl Single_Bounce\n Ns 10\n Ka 0.000000 0.000000 1.000000\n Kd 0.000000 0.000000 1.000000\n Ks 0.000000 0.000000 1.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
            fprintf(fid2,'newmtl Double_Bounce\n Ns 10\n Ka 0.000000 1.000000 0.000000\n Kd 0.000000 1.000000 0.000000\n Ks 0.000000 1.000000 0.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
            fprintf(fid2,'newmtl Triple_Bounce\n Ns 10\n Ka 1.000000 0.000000 0.000000\n Kd 1.000000 0.000000 0.000000\n Ks 1.000000 0.000000 0.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
            fprintf(fid2,'newmtl Fourfold_Bounce\n Ns 10\n Ka 1.000000 0.000000 1.000000\n Kd 1.000000 0.000000 1.000000\n Ks 1.000000 0.000000 1.000000\n Ni 1.000000\n d 1.000000\n illum 2\n\n');
            fprintf(fid2,'newmtl Fivefold_Bounce\n Ns 10\n Ka 0.000000 1.000000 1.000000\n Kd 0.000000 1.000000 1.000000\n Ks 0.000000 1.000000 1.000000\n Ni 1.000000\n d 1.000000\n illum 2');

            % loop over all points
            for i = 1:size(M_ROI,1)

                % limits on x-axis
                x_min = M_ROI(i,1)-c_size;
                x_max = M_ROI(i,1)+c_size;

                % limits on y-axis
                y_min = M_ROI(i,2)-c_size;
                y_max = M_ROI(i,2)+c_size;

                % limits on z-axis
                z_min = M_ROI(i,3)-c_size;
                z_max = M_ROI(i,3)+c_size;

                m_ind = [1 2 6 5; 4 1 5 8; 4 3 7 8; 3 2 6 7; 4 3 2 1; 7 6 5 8]+(i-1)*8;

                % Writing process
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_min,z_max);
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_min,z_max);
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_max,z_max);
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_max,z_max);
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_min,z_min);
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_min,z_min);
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_min,y_max,z_min);
                fprintf(fid,'v %7.3f %7.3f %7.3f\n',x_max,y_max,z_min);

                if M_ROI(i,4) == 1
                    fprintf(fid,'usemtl Single_Bounce\n');
                end

                if M_ROI(i,4) == 2
                    fprintf(fid,'usemtl Double_Bounce\n');
                end

                if M_ROI(i,4) == 3
                    fprintf(fid,'usemtl Triple_Bounce\n');
                end

                if M_ROI(i,4) == 4
                    fprintf(fid,'usemtl Fourfold_Bounce\n');
                end

                if M_ROI(i,4) == 5
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

            msgbox({'Export of intersection points as 3D box models accomplished (see folder Intersect including: .obj file, .mtl file).'},'3D Model Export','Help');

        end
    end
end
    





