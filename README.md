# Off-Grid DS-Spread Channel Estimation using OFDM Modulation

 - We estimate the DS-spread off-grid channel using OFDM modulation in this code.
 - For this, the main file is: *off_grid_CE_OFDM.m*
 - As mentioned in the paper, we have two methods for channel estimation (CE):
		 1. First-order approximation based VB method, i.e., FVB in the paper- for this, the function is *FVB.m*
		 2.  Second-order approximation based VB method, i.e., SVB in the paper- for this, the function is *FVB.m*


# Supporting Files

We have the following supporting files to run the code:
-The file named *gen_channel_mat.m* generates the DS-spread channel matrix.
- The file named *gen_dictionary_mat.m* generates the dictionary matrix for CE
- The file named *omp.m* estimates the sparse channel vector using the OMP algorithm


# FVB Function

The function *FVB.m* generates the DS-spread channel parameters: complex channel gain, delay, and scale parameters of each channel path using the first-order approximation method as mentioned in the paper
# SVB Function

The function *SVB.m* generates the DS-spread channel parameters: complex channel gain, delay, and scale parameters of each channel path using the second-order approximation method as mentioned in the paper

# Results

The code will give the normalised mean-square error (NMSE) vs. signal-to-noise ratio (SNR) in the estimated channel matrix
