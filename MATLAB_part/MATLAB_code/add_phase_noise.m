function [M_new] = add_phase_noise(M,phase_noise,num_r,num_a)

% function file for adding phase noise to complex image

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% input:
% M: complex image matrix
% phase_noise: maximum phase noise [rad]
% num_r: number of pixels in range
% num_a: number of pixels in azimuth

% output:
% M_new: noisy complex image matrix

% phase of input image
phi_1 = atan2(imag(M),real(M));

% amplitude of input image
A = sqrt(real(M).^2+imag(M).^2);

% add phase noise to original phase
phase_matrix = rand(num_r,num_a)*phase_noise;
phi_2 = phi_1+phase_matrix;

% calculate new signal
M_new = A.*cos(phi_2)+1i.*A.*sin(phi_2);

