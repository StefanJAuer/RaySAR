function M_sc = scale_M(M,sc_f)

% function file for scaling matrix

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% input: 
% M: point matrix to be scaled [m]
% sc_f: vector containing scaling factors (along x,y,z-axis) [dimensionless]

% output:
% M_sc: scaled point matrix [m]
%--------------------------------------------------------------------------

% Get size of input matrix
size_M = size(M);

% Create empty matrix M_sc
M_sc = zeros(size_M(1),size_M(2));

% Perform scaling
M_sc(:,1) = M(:,1).*sc_f(1);
M_sc(:,2) = M(:,2).*sc_f(2);
M_sc(:,3) = M(:,3).*sc_f(3);