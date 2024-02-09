function[M_trans,flag]= Check_Trans(trans_POV,flag)

% function for analyzing string containing transformation information

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% input: 
% - trans_POV: string containing transformation commands of syntax x<a b c>
%   (separated by semicolon)

%    --> x: transformation method ("rotate", "translate", "scale") 
%    --> a, b, c: rotation angles [degree], scaling factors [dimensionless], shift values for translation [m]
% - flag: flag for checking errors


% output:
% - M_trans [N x 4]: matrix containing transformation data
%   --> columns 1 - 3: rotation angles [degree], scaling factors [dimensionless], shift values for translation [m]
%   --> column 4: index indicating the type of transformation (0: rotation, 1: translation, 2: scaling)
% - flag: extraction successful? --> 0 = no, 1 = yes

%--------------------------------------------------------------------------

% A.) Count number of entries

% Get indices representing ";" in given string
str_ind_temp = strfind(trans_POV,';'); % Seek for semicolons
l_str = length(str_ind_temp)+1;
str_ind = zeros(l_str); % create empty dummy matrix
str_ind(2:l_str) = str_ind_temp; % index 0 added as additional starting position

% Define default matrix M_trans
M_trans = zeros(l_str,4);

%--------------------------------------------------------------------------

% B.) Fill matrix M_trans with entries

for j = 1:l_str
      
    % Next starting position
    index_test = str_ind(j)+1;
    
    % move on until first entry
    while trans_POV(index_test) == ' '
        index_test = index_test + 1;
    end
       
    %----------------------------------------------------------------------
    % 1.) Get transformation type
    if trans_POV(index_test) == 'r'
        type = 0;
    end
    
    if trans_POV(index_test) == 't'
        type = 1;
    end
    
    if trans_POV(index_test) == 's'
        type = 2;
    end
    
    %----------------------------------------------------------------------
    % 2.) First entry: parameter a
    
    % move on until first brace sign
    while trans_POV(index_test) ~= '<'
        index_test = index_test + 1;
    end
    
    % Get Index of first number
    index_start = index_test + 1;
    
    % move on until first comma
    while trans_POV(index_test) ~= ','
        index_test = index_test + 1;
    end
    
    % Get Index of last number
    index_end = index_test - 1;
    
    % Extract parameter a
    a = str2double(trans_POV(index_start:index_end));
    
    %----------------------------------------------------------------------
    % 3.) Second entry: parameter b
    
    index_test = index_test + 1;
    index_start = index_test;
    
    % move on until second comma
    while trans_POV(index_test) ~= ','
        index_test = index_test + 1;
    end
    
    % Get Index of last number
    index_end = index_test - 1;
    
    % Extract parameter b
    b = str2double(trans_POV(index_start:index_end));
    
    %----------------------------------------------------------------------
    % 4.) Third entry: parameter c
    
    index_test = index_test + 1;
    index_start = index_test;
    
    % move on until final brace
    while trans_POV(index_test) ~= '>'
        index_test = index_test + 1;
    end
    
    % Get Index of last number
    index_end = index_test - 1;
    
    % Extract parameter c
    c = str2double(trans_POV(index_start:index_end));
    
    %----------------------------------------------------------------------
    % 5.) Fill matrix with entries
    
    if exist('type','var') == 0
        
        msgbox('Error: type of transformation is undefined. Please check entries in input mask.','Unknown Type','Error');
        
        % Set flag to zero --> error occurred
        flag = 0;
        
    else
    
        % caution: 
        % - reverse direction --> compensate transformations 
        % - consider POV Ray coordinate system --> y-axis of POV Ray
        %   corresponds to z-Axis of model coordinate system

        M_trans(l_str+1-j,:) = [-a -b -c type];    
        clear type;
    end
end

