function M_tr = translate_M(M,transf)

% function file for translating matrix

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% input: 
% M: matrix to be translated
% transf: vector containing transformation information (along x,y,z-axis) [m]

% output:
% M_tr: scaled matrix
%--------------------------------------------------------------------------

% Get size of input matrix
size_M = size(M);

% Create empty matrix M_sc
M_tr = zeros(size_M(1),size_M(2));

% Perform scaling
M_tr(:,1) = M(:,1)+transf(1);
M_tr(:,2) = M(:,2)+transf(2);
M_tr(:,3) = M(:,3)+transf(3);