Sensitivity Analysis
====================

Experiment to understand the sensitivity of the 107 parameters in the tremor propagation model.
Initially created by Paul Curtis, expounded and added upon by Steven Zinn.

How to use:
-----------
After downloading this project, you will notice there are two subfolders: `functions` and `mats`. In the functions subfolder there are various function files that are used
extensively in the main code for plotting and performing the sensitivity of the tremor propagation model. To use the functions, it is imperative that you ensure that they are
in the matlab project's path. To do this:
1. Right click on the `functions` folder
2. Hover over `Add to path` in the drop down
3. Click on `Selected Folders`

The mat files are already integrated. No changes needed.

To perform a sensitivity analysis of the tremor propagation model including the 4th submodel (translation of joint displacement to hand position), go to `AIE_X_Sensitivity.m` and 
`SIE_X_Sensitivity.m`. These take about 5 minutes to run. To plot the data, go to `plot_x.mlx`. 

Files [TO UPDATE]:
------
* AIE_Means.mat - Data file for the all input excitation MEAN sensitivity with the 3rd submodel implemented
* AIE_Plot.m - Function to plot all input excitation sensitivity model with the 3rd submodel
* AIE_Sensitivities.mat - Data file for the all input excitation sensitivity with the 3rd submodel implemented
* AIE_X_Senstivities.mat - Data file for the aie sensitivity with the 4th submodel implemented
* AIE_X_Sensitivity.m - Code to run all input excitation sensitivity calculation with the 4th submodel implemented
* All_Input_Excitation_Sensitivity.m - Code to run all input excitation sensitivity calculation with the 3rd submodel implemented
* G.mat - Data file with transfer function of tremor model up to 3rd submodel
* G_Hist.mat - Data file with comparative data of numerical differentiation vs standard symbolic differentiation calculation of sensitivities
* Jacobian.mat - Jacobian (3x7 matrix) to implement in tremor model to create the 4th submodel (conversion to hand position in xyz)
* Parameters.m - Parameters and other useful variable data for sensitivity calculation
* plot_sensitivities.mlx - Live script to plot the sensitivities
* plot_x.mlx - Live script to plot the sensitivity data of the 4th submodel
* SIE_Means.mat - Data file for the single input excitation MEAN sensitivity with the 3rd submodel implemented
* SIE_Means_Plot.m - Function to plot SIE Mean calculation
* SIE_Plot.m - Function to plot SIE calculation
* SIE_Sensitivities.mat - Data file for the SIE sensitivities with the 3rd submodel implemented
* Single_Input_Excitation_Sensitivity.m - Code to run single input excitation sensitivity calculation with the 3rd submodel
* SIE_X_Sensitivity.m - Code to run single input excitation sensitivity calculation with the 4th submodel implemented
* SIE_X_Sensitivities.mat - Data file with single input excitation sensitivity calculations with the 4th submodel
* SIE_X_Means.mat - Data file with SIE means with 4th submodel

Information for the ignorant:
-----------------------------
The `X` indicates the vector x, which represents the 3 dimensions in space. If there is no `X`, then we are using the tremor model that ends in joint displacement. 
If the `X` is present in the name of a file, then we are using the tremor model that ends in hand position.

These are the link-lengths used for the jacobian for a 50th percentile male (with height of 1.763 m):
Length upper arm (Lua): 0.2853 m
Length forearm (Lfa): 0.2723 m
Length hand: 0.0873 m
Distance from axis (wrist) to center of mass of hand (Lh): 0.069 m

TODO:
-----
* Check posture 1 (with Dr. Charles)
* Check work with Dr. Charles
* Update README.md