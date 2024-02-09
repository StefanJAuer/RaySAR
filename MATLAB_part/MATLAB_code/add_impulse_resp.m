function [Im_conv] = add_impulse_resp(Im,num_a,num_r,azimuth_pixel,range_pixel,res_a,res_r)

% function adding SAR system impulse response

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% input:
% - Im: image matrix to be convolved
% - num_a: image size in azimuth [pixels]
% - num_r: image size in range [pixels]
% - azimuth_pixel: pixel spacing in azimuth [m]
% - range_pixel: pixel spacing in range [m]
% - res_a: azimuth resolution [m]
% - res_r: range resolution [m]

% output:
% - Im_conv: image matrix convolved with 2D SAR system response

%--------------------------------------------------------------------------

% Signal properties
fs_a = 1/azimuth_pixel; % samples per meter in azimuth [Hz]
fs_r = 1/range_pixel; % samples per meter in range [Hz]
Dur_a = num_a*azimuth_pixel; % duration of signal in azimuth [m] 
Dur_r = num_r*range_pixel; % duration of signal in range [m] 

f_a = -fs_a/2 : 1/Dur_a : -fs_a/2+(num_a-1)*1/Dur_a; % frequency steps in azimuth [Hz]
f_r = -fs_r/2 : 1/Dur_r : -fs_r/2+(num_r-1)*1/Dur_r; % frequency steps in range [Hz]

% Periods in azimuth and range
T_a = res_a;
T_r = res_r;

% Maximum frequencies of system transfer function
% --> cut off frequency for hamming window
fg_a = 1/T_a; % [Hz] 
fg_r = 1/T_r; % [Hz] 

%--------------------------------------------------------------------------

% 2D-System transfer function

% Separate definition in azimuth and range
RECT_a = rectpuls(f_a,fg_a);
RECT_b = rectpuls(f_r,fg_r);

A = zeros(1,length(RECT_b));
index_A = (RECT_b == 1);
A(index_A)=1;
B = zeros(1,length(RECT_a));
index_B = (RECT_a == 1);
B(index_B)=1;

% 2D-System transfer function in azimuth and range (span 2D image space
% using vectors A and B)
CONV_MAP_sort = A'*B;

% Unsort in frequency domain
CONV_MAP_temp = fftshift(CONV_MAP_sort,1);
CONV_MAP_temp = fftshift(CONV_MAP_temp,2);

%--------------------------------------------------------------------------

% Create hamming window

% Hamming function
alpha = 0.54;
HA_f = alpha+(1-alpha)*cos(2*pi.*(f_a./fg_a));
HR_f = alpha+(1-alpha)*cos(2*pi.*(f_r./fg_r));

% 2D-hamming window
HAM = HR_f'*HA_f;

% Unsort in frequency domain
HAM_temp = fftshift(HAM,1);
HAM_temp = fftshift(HAM_temp,2);

%--------------------------------------------------------------------------

% Image processing (performed in frequency domain)

% Fourier Transform of original image data
IM = fft2(Im);

% Multiplication with SAR system transfer function
IM_CONV = IM.*CONV_MAP_temp;

% Apply hamming window
IM_CONV = IM_CONV.*HAM_temp;

% scaling factors
m = size(IM_CONV,1);
n = size(IM_CONV,2);

% Inverse Fourier Transform
Im_conv = ifft2(IM_CONV,m,n);
Im_conv = abs(Im_conv); % absolute value
%--------------------------------------------------------------------------