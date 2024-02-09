# function files for simulating SAR images

# Stefan Auer, DLR-IMF, 2017

import os
import numpy as np
import time
import scipy as sp
from scipy import weave  # For inline c code
from scipy.weave import converters
from subprocess import Popen, PIPE, check_call
import pandas as pd


def Gen_Refl_Map(go, ang, r_geom, range_dir, Output_path):

    # function-file for creating reflectivity maps, maps for specular analysis and intensity distribution maps

    #--------------------------------------------------------------------------

    # Stefan Auer
    # Remote Sensing Technology Institute
    # German Aerospace Center (DLR)
    # 2017

    #--------------------------------------------------------------------------

    # from raysar_linux_function

    global Az, Ra, Intens, Tr_L

    # local parameters:

    # Elements within class "go":
    # - go.a_pix: azimuth pixel spacing of reflectivity map to be generated [m]
    # - go.r_pix: range pixel spacing of reflectivity map to be generated [m]
    # - go.a_min, go.a_max, go.r_min, go.r_max: image limits on azimuth and range axis [m]
    # - go.db_min, go.db_max: power range to be displayed [dB in string format]
    # - go.clip: threshold for clipping of image intensities [string]
    # - go.filt: flag for applying binomial filter [no: value 0; yes: value 1]
    # - go.bounce_level: maximum bounce level to be analyzed [value between 1 and 5]
    # - go.specular_flag: flag for preparing specular maps [no --> value 0; yes --> value 1]
    # - go.sinc_sim: flag for convolving the simulated image with the SAR
    #   system response [no --> value 0; yes --> value 1]
    # - go.a_res,go.r_res: resolution in azimuth and range (3dB width of 2D impulse response)
    # - go.coh: flag for coherent summation of signal contributions (0 = absolute values only; 1 = summation of complex values)
    # - go.dBnodB: flag for deciding between dB amplitude scale and clipping [value 0: clipping, value 1: dB]
    # - go.store: flag for storing image data to file [value 0: do not store pixel values, value 1: store image matrix to .mat file]

    # Output:
    # - geometrical distribution of scatterers in the azimuth-range plane
    # - reflectivity map and separate layers for different reflection levels
    #   (clipped or dB)
    # - specular map: binary map marking pixels containing specular signal
    #   contributions (angle tolerance: +-10)

    #--------------------------------------------------------------------------

    # Az: azimuth coordinates of reflection contributions [m]
    # Ra: range coordinates of -''- [m]
    # Intens: intensity of -''- [value between 0 and 1]
    # Tr_L: bounce level of -''- [value between 1 and 5]
    # ang: incidence angle of simulated signal [rad]
    # r_geom: image geometry in range direction; 0: slant range, 1: ground range
    # Output_path: absolute path to designated folder where simulation results are to be stored

    #--------------------------------------------------------------------------

    from scipy import signal

    print 'Generate SAR image layers.'

    # Time specification
    SchrittNr = 0
    timeM = np.zeros(10)
    timeM[SchrittNr] = time.time()

    # Delete older result folder
    path_results = os.path.join(Output_path, 'Maps')

    if os.path.exists(path_results):
        check_call('rm -rf %s' % path_results, shell=True)

    if r_geom == 0:
        # Slant range geometry
        Ra_g = Ra
    else:
        # Change to ground range geometry
        Ra_g = Ra * (1 / np.sin(ang))  # ground range

    # Number of lines
    s = np.size(Az)

    # Signal wavelength (X-Band)
    wavelength = 0.031  # [m]

    #--------------------------------------------------------------------------

    # 1.) Set image properties

    #--------------------------------------------------------------------------

    # Set boundaries of image

    # minimum in azimuth
    if go.a_min == 'Min':
        go.a_min = round(min(Az))
    else:
        go.a_min = float(go.a_min)

    # maximum in azimuth
    if go.a_max == 'Max':
        go.a_max = round(max(Az))
    else:
        go.a_max = float(go.a_max)

    # minimum in range
    if go.r_min == 'Min':
        go.r_min = round(min(Ra_g))
    else:
        go.r_min = float(go.r_min)

    # maximum in range
    if go.r_max == 'Max':
        go.r_max = round(max(Ra_g))
    else:
        go.r_max = float(go.r_max)

    # Maximum bounce level to be displayed
    if go.bounce_level == 'Max':
        go.bounce_level = max(Tr_L)
    else:
        go.bounce_level = float(go.bounce_level)

    #--------------------------------------------------------------------------

    # Save image limits globally for tomographic analysis
    az_min = go.a_min  # Minimum in Azimuth [m]
    az_max = go.a_max  # Maximum in -''- [m]
    ra_min = go.r_min  # Minimum in Slant-Range/Ground-Range [m]
    ra_max = go.r_max  # Maximum in -''- [m]
    # Maximum bounce level to be analyzed [value between 1 and 5]
    bounce_rem = go.bounce_level

    #--------------------------------------------------------------------------

    # Intervals on axes
    diff_a = go.a_max - go.a_min  # interval azimuth
    diff_r = go.r_max - go.r_min  # interval slant range

    num_a = np.ceil(diff_a / go.a_pix)  # azimuth
    num_r = np.ceil(diff_r / go.r_pix)  # range

    #--------------------------------------------------------------------------

    # 2.) Predefine size of image matrices

    #--------------------------------------------------------------------------

    # Reflectivity Maps

    # For gathering single bounce contributions
    Single = np.zeros((num_r, num_a), 'float')

    if go.bounce_level > 1:
        Double = np.zeros((num_r, num_a), 'float')  # double ...

    #--------------------------------------------------------------------------

    # 3.) Fill values into image matrices

    #--------------------------------------------------------------------------

    # Set parameters for inline C-processing
    a_min = go.a_min
    a_max = go.a_max
    r_min = go.r_min
    r_max = go.r_max
    a_pix = go.a_pix
    r_pix = go.r_pix
    bounce_level = go.bounce_level
    num_a = int(num_a)
    num_r = int(num_r)

    # Define inline code for C-processing
    # --> summing of signal contributions in matrix

    code = '''

	int k, row_pix_float, column_pix_float;
	int row_pix, column_pix;
	double pi, Az_temp, Ra_temp, Tr_L_temp, Intens_temp, Phi_temp, Phi, Cycles, Signal;

	// Loop over all contributions
	for (k=0; k<s; k++)
	{

		// Assign signal contributions to image layers

		Az_temp = Az(k);
		Ra_temp = Ra_g(k);
		Tr_L_temp = Tr_L(k);
		Intens_temp = Intens(k);

		//----------------------------------------------------------------------

		// If contribution is located within defined margins
		if (Az_temp >= a_min && Az_temp <= a_max)
		{
			if (Ra_temp >= r_min && Ra_temp <= r_max)
			{

				//--------------------------------------------------------------

				// B.)Pixel center
				// --> rounded coordinate in azimuth + 0.5
				// --> rounded coordinate in range + 0.5

				//e.g. for first pixel: row_pix = 0.5, column_pix = 0.5

				row_pix_float = floor((Ra_temp-r_min)/r_pix)+0.5;
				column_pix_float = floor((Az_temp-a_min)/a_pix)+0.5;

				//--------------------------------------------------------------

				// D.) Add contribution to pixels

				// Round pixel coordinates for obtaining matrix indices
				column_pix_float = round(column_pix_float);
				row_pix_float = round(row_pix_float);

				// Change variable type to int
				column_pix = (int) column_pix_float;
				row_pix = (int) row_pix_float;

				if (row_pix <= num_r && column_pix <= num_a)
				{
					// Absolute
					Signal = Intens_temp;

					// Single bounce
					if (Tr_L_temp == 1)
					{
						Single(row_pix-1,column_pix-1) = Single(row_pix-1,column_pix-1)+Signal;
					}

					// Double bounce
					if (Tr_L_temp == 2 && bounce_level > 1)
					{
						Double(row_pix-1,column_pix-1) = Double(row_pix-1,column_pix-1)+Signal;
					}
				}
			}
		}
	}
	'''

    err = weave.inline(code, ['Az', 'Ra_g', 'Tr_L', 'Intens', 'a_min', 'a_max', 'r_min', 'r_max', 'num_r', 'num_a', 'a_pix', 'r_pix', 'bounce_level', 's', 'Single', 'Double'],
                       type_converters=converters.blitz)

    # Time specification
    SchrittNr += 1
    timeM[SchrittNr] = time.time()

    print "\nTime for generating image matrices:", np.round(timeM[SchrittNr] - timeM[SchrittNr - 1], 2), "seconds"

    #--------------------------------------------------------------------------

    # Sum up signals

    # Summing up signals according to chosen bounce maximum
    # in case of imaginary numbers: maps are added coherently --> absolute
    # values are calculated afterwards

    All_reflections = Single

    if go.bounce_level > 1:
        # Sum of signals
        All_reflections = All_reflections + Double

    #----------------------------------------------------------------------

    # Time specification
    SchrittNr += 1
    timeM[SchrittNr] = time.time()

    print "\nTime for summing up matrices:", np.round(timeM[SchrittNr] - timeM[SchrittNr - 1], 2), "seconds"

    # Optional: Smoothing

    if go.filt == 1:

        # Predefine dummy matrices
        h = np.zeros((3, 3), 'float')

        # Fill first elements of filter matrix with entries of Binomial Filter
        h[0, 0] = 1.0 / 16.0
        h[0, 1] = (1.0 / 16.0) * 2.0
        h[0, 2] = 1.0 / 16.0
        h[1, 0] = (1.0 / 16.0) * 2.0
        h[1, 1] = (1.0 / 16.0) * 4.0
        h[1, 2] = (1.0 / 16.0) * 2.0
        h[2, 0] = 1.0 / 16.0
        h[2, 1] = (1.0 / 16.0) * 2.0
        h[2, 2] = 1.0 / 16.0

        # Smooth images by 3x3 Binomial-Filter

        Single = signal.convolve2d(Single, h, mode='same')

        if 'Double' in locals():
            Double = signal.convolve2d(Double, h, mode='same')

        if 'All_reflections' in locals():
            All_reflections = signal.convolve2d(
                All_reflections, h, mode='same')

    # Time specification
    SchrittNr += 1
    timeM[SchrittNr] = time.time()

    print "\nTime for matrix filtering:", np.round(timeM[SchrittNr] - timeM[SchrittNr - 1], 2), "seconds"

    #----------------------------------------------------------------------

    # Optional: Add system response function

    # if go.sinc_sim == 1:
    #
    # # Single Bounce
    # if 'Single' in locals():
    # 	Single = add_impulse_resp(Single, num_a, num_r, go.a_pix, go.r_pix, go.a_res, go.r_res)
    #
    # # Double Bounce
    # if 'Double' in locals():
    # 	Double = add_impulse_resp(Double, num_a, num_r, go.a_pix, go.r_pix, go.a_res, go.r_res)
    #
    # # Reflectivity Map
    # if 'All_reflections' in locals():
    # 	All_reflections = add_impulse_resp(All_reflections, num_a, num_r, go.a_pix, go.r_pix, go.a_res, go.r_res)

    #----------------------------------------------------------------------

    # Logarithmic scaling or clipping

    if go.dBnodB == 1:

        # Logarithmic scaling

        Single = 10 * np.log(Single)
        All_reflections = 10 * np.log(All_reflections)

        if 'Double' in locals():
            Double = 10 * np.log(Double)

        # Detect Maximum Intensity to be displayed

        # Minimum value [dB]
        if go.db_min == 'Min':
            im_min = -100
        else:  # numerical value chosen by operator
            im_min = float(go.db_min)

        # Maximum value [dB]
        if go.db_max == 'Max':  # maximum value to be taken
            if go.bounce_level > 1:
                im_max = All_reflections.max()
            else:
                im_max = Single.max()

        else:  # numberical value chosen by operator
            im_max = float(go.db_max)

        # Necessary correction for logarithmically scaled matrices
        # Set minimum within image to db_min
        # --> remove entries = -Inf

        Single[Single < im_min] = im_min
        Single[Single > im_max] = im_max

        All_reflections[All_reflections < im_min] = im_min
        All_reflections[All_reflections > im_max] = im_max

        if 'Double' in locals():
            Double[Double < im_min] = im_min
            Double[Double > im_max] = im_max

    else:

        # Clipping yes/no
        if go.clip != 'Max':

            # Maximum
            max_clip = float(go.clip)

            # Avoid threshold bigger than given intensity maximum
            if max_clip <= All_reflections.max() and max_clip >= 0:

                # Find outliers in image layers
                index_single_clip = np.where(Single > max_clip)
                index_single_clip = index_single_clip[0]

                index_all_reflections_clip = np.where(
                    All_reflections > max_clip)
                index_all_reflections_clip = index_all_reflections_clip[0]

                Single[index_single_clip, :] = max_clip
                All_reflections[index_all_reflections_clip, :] = max_clip

                if 'Double' in locals():

                    index_double_clip = np.where(Double > max_clip)
                    index_double_clip = index_double_clip[0]

                    Double[index_double_clip, :] = max_clip

                # Clipping applied
                clipping_status = 1
            else:
                # No clipping applied
                clipping_status = 0

        # Extreme intensities
        im_min = All_reflections.min()
        im_max = All_reflections.max()

    # Size of images to be displayed
    r_p = np.linspace(0, np.size(Single, 0) - 1,
                      np.size(Single, 0))  # rows --> range
    a_p = np.linspace(0, np.size(Single, 1) - 1,
                      np.size(Single, 1))  # columns --> azimuth

    # Time specification
    SchrittNr += 1
    timeM[SchrittNr] = time.time()

    print "\nTime for logarithmic scaling:", np.round(timeM[SchrittNr] - timeM[SchrittNr - 1], 2), "seconds"

    #--------------------------------------------------------------------------

    # 5. Figures

    #--------------------------------------------------
    # Save image to file

    # Create new folder
    if os.path.exists(path_results) != 1:
        check_call('mkdir %s' % path_results, shell=True)

    # Create sub-folders
    check_call('mkdir %s' % path_results + '/Frames', shell=True)
    check_call('mkdir %s' % path_results + '/Frames/Ref_Maps', shell=True)

    #--------------------------------------------------------------------------

    # B.) Reflectivity maps

    #--------------------------------------------------------------------------

    # Single bounce

    Bounce = 'Single Bounce'
    R_Map(Single, im_min, im_max, Bounce, range_dir, Output_path)

    #--------------------------------------------------------------------------

    # Double bounce
    if go.bounce_level > 1:

        Bounce = 'Double Bounce'
        R_Map(Double, im_min, im_max, Bounce, range_dir, Output_path)

    #--------------------------------------------------------------------------

    # Display all bounces in one plot
    if go.bounce_level > 1:

        Bounce = 'All Reflections'
        R_Map(All_reflections, im_min, im_max, Bounce, range_dir, Output_path)

    # Time specification
    SchrittNr += 1
    timeM[SchrittNr] = time.time()

    print "\nTime for generating image frames:", np.round(timeM[SchrittNr] - timeM[SchrittNr - 1], 2), "seconds"

    #--------------------------------------------------------------------------


def R_Map(Im_Ma, im_min, im_max, Bounce, range_dir, Output_path):

    # function for displaying reflectivity maps

    #--------------------------------------------------------------------------

    # Stefan Auer
    # Remote Sensing Technology Institute
    # German Aerospace Center (DLR)
    # 2017

    #--------------------------------------------------------------------------

    # local parameters

    # input:
    # - Im_Ma: image matrix containing pixel amplitudes
    # - im_min, im_max: power range in amplitude to be displayed (dB values or clipping threshold)
    # - Bounce: string containing bounce level information
    # - range_dir: flag for marking range direction (0: bottom up, 1: top down)
    # - Output_path: path to output folder

    # output: none

    #--------------------------------------------------------------------------

    # global parameters

    # r_geom: image geometry in range direction; 0: slant range, 1: ground range
    # range_dir: range direction within simulated images; value 0: bottom up; value 1: top down
    # Output_path: absolute path to designated folder where simulation results are to be stored

    #--------------------------------------------------------------------------

    import scipy.misc

    # Scale image entries to interval of gray values
    num_g = 255  # maximum gray value --> 255
    Im_Ma = im_trans(Im_Ma, im_min, im_max, num_g)

    if range_dir == 0:
        Im_Ma = np.flipud(Im_Ma)  # flip up down
    else:
        Im_Ma = np.fliplr(Im_Ma)  # flip left right

    # Store image in folder
    export_string = Output_path + '/Maps/Frames/Ref_Maps/' + \
        Bounce + '_Fr.tif'  # export string
    sp.misc.imsave(export_string, Im_Ma)  # store image matrix to tif file


def im_trans(Im_in, g_min, g_max, num_g):

    # function file for scaling image values to gray values

    # Stefan Auer
    # Remote Sensing Technology Institute
    # German Aerospace Center (DLR)
    # 2017

    # Input:
    # - Im_in: input image [dB]
    # - g_min: minimum intensity value to be displayed [dB]
    # - g_max: maximum intensity value to be displayed [dB]
    # - num_g: number of gray values to defined

    # Output:
    # - Image scaled to defined interval of gray values

    # Maximum interval
    interval_im = g_max - g_min

    # Find step width for grey values
    step_width = num_g / interval_im

    # Shift to zero --> lowest gray value will be 0
    Im_in = Im_in - g_min

    # Transform to gray values
    Im_out = Im_in * step_width
    Im_out = np.round(Im_out)

    return Im_out


def raysar_function(contribution_filename, Output_path, simulation_parameter_filename):

    # function file for initiating image generation

    # Stefan Auer
    # Remote Sensing Technology Institute
    # German Aerospace Center (DLR)
    # 2017

    global Az, Ra, Intens, Tr_L

    print 'Interpret contribution file.'

    # Initialize parameter class
    go = SimSettings()

    text_file = open(simulation_parameter_filename, 'r')
    lines = text_file.readlines()

    # Get number of entries
    num_max = len(lines)

    for i in range(0, num_max):
        line = lines[i]
        split_result = line.split('=', 2)

        if i == 0:
            go.a_pix = float(split_result[1])  # pixel spacing azimuth

        if i == 1:
            go.r_pix = float(split_result[1])  # pixel spacing range

        if i == 2:
            go.a_min = float(split_result[1])  # mimimum azimuth

        if i == 3:
            go.a_max = float(split_result[1])  # maximum azimuth

        if i == 4:
            go.r_min = float(split_result[1])  # minimum range

        if i == 5:
            go.r_max = float(split_result[1])  # maximum range

        if i == 6:
            go.bounce_level = float(split_result[1])  # reflection level

        if i == 7:
            go.db_min = float(split_result[1])  # dB Min

        if i == 8:
            go.db_max = 'Max'  # dB Max

        if i == 9:
            range_dir = float(split_result[1])  # Range direction

        if i == 10:
            go.filt = float(split_result[1])  # Flag for binomial filtering

        if i == 11:
            # Flag for range geometry (-> ground range)
            r_geom = float(split_result[1])

        if i == 12:
            ang = float(split_result[1]) * (np.pi / 180.0)  # Incidence angle

    # Further settings
    go.clip = 'Max'
    go.dBnodB = 1
    go.coh = 0

    data_df = pd.read_csv(contribution_filename, delimiter=' ', header=None)

    # Re-Name data columns
    data_df.columns = ['Az', 'Ra', 'El', 'Intens', 'Tr_L', 'Sp', 'Blah']

    Az = data_df.Az.values  # azimuth coordinates
    Ra = data_df.Ra.values  # range coordinates
    Intens = data_df.Intens.values  # intensities
    Tr_L = data_df.Tr_L.values  # reflection level

    # Generate SAR image(s)
    Gen_Refl_Map(go, ang, r_geom, range_dir, Output_path)

    return 1


class SimSettings(object):

    # Class for steering SAR simulation
    a_min = ''
    a_max = ''
    r_min = ''
    r_max = ''
    db_min = ''
    db_max = ''
    bounce_level = ''
    clip = ''
    filt = ''
    specular_flag = ''
    sinc_sim = ''
    a_res = ''
    r_res = ''
    coh = ''
    dBnodB = ''


# def add_impulse_resp(Im, num_a, num_r, azimuth_pixel, range_pixel, res_a, res_r):
#
# 	# function file for adding SAR system impulse response to reflectivity map
#
# 	#--------------------------------------------------------------------------
#
# 	# Stefan Auer
#   # Remote Sensing Technology Institute
#   # German Aerospace Center (DLR)
#   # 2017
#
# 	#--------------------------------------------------------------------------
#
# 	# input:
# 	# - Im: image matrix to be convolved
# 	# - num_a: image size in azimuth [pixels]
# 	# - num_r: image size in range [pixels]
# 	# - azimuth_pixel: pixel spacing in azimuth [m]
# 	# - range_pixel: pixel spacing in range [m]
# 	# - res_a: azimuth resolution [m]
# 	# - res_r: range resolution [m]
#
# 	# output:
# 	# - Im_conv: image matrix convolved with 2D SAR system response
#
# 	#--------------------------------------------------------------------------
#
# 	# Signal properties
# 	fs_a = 1/azimuth_pixel  # samples per meter in azimuth [Hz]
# 	fs_r = 1/range_pixel  # samples per meter in range [Hz]
# 	Dur_a = num_a*azimuth_pixel  # duration of signal in azimuth [m]
# 	Dur_r = num_r*range_pixel  # duration of signal in range [m]
#
# 	f_a = np.linspace(-fs_a/2, (fs_a/2)-1/Dur_a, num_a)
# 	f_r = np.linspace(-fs_r/2, (fs_r/2)-1/Dur_r, num_r)
#
# 	#f_a = -fs_a/2 : 1/Dur_a : (fs_a/2)-1/Dur_a  # frequency steps in azimuth [Hz]
# 	#f_r = -fs_r/2 : 1/Dur_r : (fs_r/2)-1/Dur_r  # frequency steps in range [Hz]
#
# 	# Periods in azimuth and range
# 	T_a = res_a
# 	T_r = res_r
#
# 	# Maximum frequencies of system transfer function
# 	# --> cut off frequency for hamming window
# 	fg_a = 1/T_a  # [Hz]
# 	fg_r = 1/T_r  # [Hz]
#
# 	#--------------------------------------------------------------------------
#
# 	# 2D-System transfer function
#
# 	# Separate definition in azimuth and range
# 	RECT_a = rectpuls(f_a, fg_a)
# 	RECT_b = rectpuls(f_r, fg_r)
#
# 	A = np.zeros((1, len(RECT_b)))
# 	index_A = (RECT_b == 1)
# 	A(index_A)=1
# 	B = np.zeros((1, len(RECT_a)))
# 	index_B = (RECT_a == 1)
# 	B(index_B)=1
#
# 	# 2D-System transfer function in azimuth and range (span 2D image space
# 	# using vectors A and B)
# 	CONV_MAP_sort = A'*B
#
# 	# Unsort in frequency domain
# 	CONV_MAP_temp = np.fft.fftshift(CONV_MAP_sort, 1)
# 	CONV_MAP_temp = np.fft.fftshift(CONV_MAP_temp, 2)
#
# 	#--------------------------------------------------------------------------
#
# 	# Create hamming window
#
# 	# Hamming function
# 	alpha = 0.54
# 	HA_f = alpha+(1-alpha)*np.cos(2* np.pi*(f_a/fg_a))
# 	HR_f = alpha+(1-alpha)*np.cos(2* np.pi*(f_r/fg_r))
#
# 	# 2D-hamming window
# 	HAM = HR_f'*HA_f
#
# 	# Unsort in frequency domain
# 	HAM_temp = np.fft.fftshift(HAM, 1)
# 	HAM_temp = np.fft.fftshift(HAM_temp, 2)
#
# 	#--------------------------------------------------------------------------
#
# 	# Image processing (performed in frequency domain)
#
# 	# Fourier Transform of original image data
# 	IM = np.fft.fft2(Im)
#
# 	# Multiplication with SAR system transfer function
# 	IM_CONV = IM*CONV_MAP_temp
#
# 	# Apply hamming window
# 	IM_CONV = IM_CONV*HAM_temp
#
# 	# scaling factors
# 	m = np.size(IM_CONV, 1)
# 	n = np.size(IM_CONV, 2)
#
# 	# Inverse Fourier Transform
# 	Im_conv = np.fft.ifft2(IM_CONV, m, n)
# 	Im_conv = abs(Im_conv)  # absolute value
#
# 	return Im_conv
