function[M_sens,flag]= Check_Sensor(sensor_POV,flag)

% function for analyzing a string containing information about the sensor geometry

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% input: 
% - sensor_POV: string containing sensor commands of syntax x<a b c>
%   (separated by semicolon)

%    --> x: position information ("location", "look_at") 
%    --> a, b, c: coordinates [m]
% - flag: flag for checking errors

% output:
% - M_sens [2 x 4]: matrix describing imaging geometry of virtual SAR sensor
%   --> columns 1 - 3: rotation angles [deg], scaling factors [dimensionless], shift values for translation [m]
%   --> column 4: index indicating the type of position (0: location, 1: look_at)
% - flag: extraction successful? --> 0 = no, 1 = yes


% A.) Count number of entries

% Get indices representing ";" in given string
str_ind_temp = strfind(sensor_POV,';'); % Seek for semicolons

if length(str_ind_temp) ~= 1
    msgbox('Error: Number of entries has to be two (location + look_at). Please check entries in input mask.','Erronenous Entries','Error');
else

    l_str = length(str_ind_temp)+1;
    str_ind = zeros(l_str); % create empty dummy matrix
    str_ind(2:l_str) = str_ind_temp; % index 0 added as additional starting position

    % Define default matrix M_sens
    M_sens = zeros(2,4);

    %--------------------------------------------------------------------------

    % B.) Fill matrix M_sens with entries

    for j = 1:2

        % Next starting position
        index_test = str_ind(j)+1;

        % move on until first entry
        while sensor_POV(index_test) == ' '
            index_test = index_test + 1;
        end

        %----------------------------------------------------------------------
        % 1.) Get position type

        % Check for entry "location"
        if sensor_POV(index_test:index_test+2) == 'loc'
            type = 0;
        end

        % Check for entry "look_at
        if sensor_POV(index_test:index_test+2) == 'loo'
            type = 1;
        end

        %----------------------------------------------------------------------
        % 2.) First entry: parameter a

        % move on until first brace sign
        while sensor_POV(index_test) ~= '<'
            index_test = index_test + 1;
        end

        % Get Index of first number
        index_start = index_test + 1;

        % move on until first comma
        while sensor_POV(index_test) ~= ','
            index_test = index_test + 1;
        end

        % Get Index of last number
        index_end = index_test - 1;

        % Extract parameter a
        a = str2double(sensor_POV(index_start:index_end));

        %----------------------------------------------------------------------
        % 3.) Second entry: parameter b

        index_test = index_test + 1;
        index_start = index_test;

        % move on until second comma
        while sensor_POV(index_test) ~= ','
            index_test = index_test + 1;
        end

        % Get Index of last number
        index_end = index_test - 1;

        % Extract parameter b
        b = str2double(sensor_POV(index_start:index_end));

        %----------------------------------------------------------------------
        % 4.) Third entry: parameter c

        index_test = index_test + 1;
        index_start = index_test;

        % move on until final brace
        while sensor_POV(index_test) ~= '>'
            index_test = index_test + 1;
        end

        % Get Index of last number
        index_end = index_test - 1;

        % Extract parameter c
        c = str2double(sensor_POV(index_start:index_end));

        %----------------------------------------------------------------------
        % 5.) Fill matrix with entries

        if exist('type','var') == 0

            msgbox('Error: position type is undefined. Please check entries in input mask.','Unknown Type','Error');

            % Set flag to zero --> error occurred
            flag = 0;

        else
            
            M_sens(j,:) = [a b c type];    
            clear type;
        end
    end
end