function M_rot = rotation_M(M,rot_ang)

% function file for rotating matrix

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% function for rotating a point cloud in 3D space

% input: 
% M: point matrix to be rotated
% rot_ang: vector containing rotation angles (rotation about x,y,z-axis) [rad]

% output:
% M_rot: rotated point matrix

% rotation order:
% 1.) rotation about x-axis, 2.) about y-axis, 3.) about z-axis
%--------------------------------------------------------------------------

% Capture rotation angles
angle_x = rot_ang(1); % unit [rad]
angle_y = rot_ang(2);
angle_z = rot_ang(3);

% Create rotation matrices

% Rotation matrix 1 (about x-axis)
R_x = [1 0 0; 0 cos(angle_x) -sin(angle_x); 0 sin(angle_x) cos(angle_x)];

% Rotation matrix 2 (about y-axis)
R_y = [cos(angle_y) 0 sin(angle_y); 0 1 0; -sin(angle_y) 0 cos(angle_y)];

% Rotation matrix 3 (about z-axis)
R_z = [cos(angle_z) -sin(angle_z) 0; sin(angle_z) cos(angle_z) 0; 0 0 1];

%--------------------------------------------------------------------------

% Perform rotations
M_rot = ((M*R_x)*R_y)*R_z; 
%--------------------------------------------------------------------------