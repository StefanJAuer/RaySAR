function [Im_out] = im_trans(Im_in,g_min,g_max,num_g)

% function file for scaling image values to gray values 

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% Input: 
% - Im_in: input image [dB]
% - g_min: minimum intensity value to be displayed [dB]
% - g_max: maximum intensity value to be displayed [dB]
% - num_g: number of gray values to defined

% Output:
% - Image scaled to defined interval of gray values

% Maximum interval
interval_im = g_max - g_min;

% Find step width for grey values
step_width = num_g/interval_im;

% Shift to zero --> lowest gray value will be 0
Im_in = Im_in - g_min;

% Transform to gray values
Im_out = round(Im_in.*step_width);



