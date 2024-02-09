function [M_red] = remove_red_inter(M)

% function for removing redundant intersection points considering bounce
% level information
% --> equal points with different bounce level information are preserved

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------
% input:
% - M: N x 4 matrix containing coordinates of intersection points (columns 1 to 3) and
%   bounce level information for each point (column 4)

% output:
% - M_red: point cloud where each point of certain bounce level only appears once 

% Number of points
s = size(M);

% Define empty matrices
M_temp1 = zeros(s(1)+1,s(2));
M_temp2 = zeros(s(1)+1,s(2));

% Sorting
M_temp1(1:s(1),:) = sortrows(M,3); % sort according to z-values
M_temp1(1:s(1),:) = sortrows(M,2); % sort according to y-values
M_temp1(1:s(1),:) = sortrows(M,1); % sort according to x-values

% Fill second matrix by shifting M_temp1
M_temp2(2:s(1)+1,:) = M_temp1(1:s(1),:);

% Calculate difference between shifted matrices
M_Diff = abs(M_temp2-M_temp1);
M_Diff(:,1:3) = M_Diff(:,1:3).^2; 

% Condition 1: Spatial distance 
Vec_Diff = sqrt(sum(M_Diff(:,1:3),2));

% Find rows fulfilling pre-condition (difference: delta [m] --> spatial distance between points)
ind_1 = find(Vec_Diff > 0);

% Condition 2: equality of bounce level
ind_2 = find(M_Diff(:,4) ~= 0);

% Indices of points to be conserved
ind = union(ind_1,ind_2);

% Reselect points of interest
M_red = M_temp1(ind(1:length(ind)-1),:);