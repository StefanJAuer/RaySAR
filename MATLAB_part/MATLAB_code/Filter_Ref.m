function[im_filt] = Filter_Ref(im)

% function file for filtering image by a Binomial filter

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

% input: 
% - im: image to be filtered (data type: double)

% output:
% - im_filt: filtered image

%--------------------------------------------------------------------------

% A.) Define filter matrix

% Get size of image
s = size(im);

% Predefine dummy matrices
h = zeros(s(1),s(2));

% Fill first elements with entries
h(1:3,1:3) = (1/16)*[1 2 1; 2 4 2; 1 2 1]; % Binomial Filter

%--------------------------------------------------------------------------

% B.) Apply smoothing

% Fourier Transformation
IM = fft2(im); 

% Filter
H = fft2(h); 

% Multiplication in frequency domain
% (corresponds to convolution in time domain)
IM_FILT = IM.*H;

% Inverse Fourier Transform
im_filt = abs(ifft2(IM_FILT));

%--------------------------------------------------------------------------


