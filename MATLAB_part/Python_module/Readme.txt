In this folder you find a Python (dummy) implementation of the core function Gen_Refl_Map (see corresponding Gen_Refl_Map.m in the MATLAB package), which generates 2D SAR images
out of signal contributions derived from the ray tracing process. 

Please note that the code only covers bounce levels 2 to start with. Accordingly, bounce levels 3-5 should be added if necessary (five is the current limit covered in the ray tracing 
process). Just cross-compare the Python code to the MATLAB code and add the missing parts (line 251: summation of signal contributions, line 279: summation of all contributions, 
line 313: convolution, lines 387 and 416: dB scaling and clipping, line 474: visualize and export reflectivity maps).

Moreover, please note that scipy.weave inline C code is unfortunately not supported by Python versions > 3. 

Finally, if you don't have access to MATLAB, also check out the compiled MATLAB package in order to test the functionalities of RaySAR.

Kind regards
Stefan Auer

