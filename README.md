Sensitivity Analysis
====================

Experiment to understand the sensitivity of the 107 parameters (now 110 with additions as discussed below) in the tremor propagation model. 
Initially created by Paul Curtis, expounded and added upon by Steven Zinn. The additions to Paul Curtis' work include:
* Numerical differentiation implementation rather than using symbolic differentiation
* Addition of a 4th submodel with Jacobian through the RVC toolbox to translate joint displacement to hand position in the xyz. This added 3 more parameters [length of hand, length of forearm, length of upper arm]
* Plots to analyze sensitivities of 110 parameters

How to use:
-----------
After downloading this project, you will notice there are two subfolders: `functions` and `mats`. In the functions subfolder there are various 
function files that are used extensively in the main code for plotting and performing the sensitivity of the tremor propagation model. To use 
the functions, it is imperative that you ensure that they are in the matlab project's path. To do this:
1. Right click on the `functions` folder
2. Hover over `Add to path` in the drop down
3. Click on `Selected Folders`

The mat files are already integrated. No changes needed.

To perform a sensitivity analysis of the tremor propagation model including the 4th submodel (translation of joint displacement to hand position), go to `AIE_X_Sensitivity.m` 
and `SIE_X_Sensitivity.m`. These take about 5 minutes to run. To plot the data, go to `plot_x.mlx`. 

Files:
------
The main files to be aware of:
* AIE_X_Sensitivity - This code runs the all input excitation for the 4th and last submodel
* AIE_Sensitivity - This code runs the aie for the 3rd submodel (joint displacement)
* Parameters - Holds important parameter information and some useful variables for global use
* plot_joint_displacement - INCOMPLETE PROJECT to plot the 3rd submodel
* plot_x - Plots the 4th submodel
* SIE_Sensitivity - This code runs the single input excitation for the 3rd submodel
* SIE_X_Sensitivity - This code runs the sie for the 4th submodel

Information for the ignorant:
-----------------------------
The `X` indicates the vector x, which represents the 3 dimensions in space. If there is no `X`, then we are using the tremor model that ends in joint displacement. 
If the `X` is present in the name of a file, then we are using the tremor model that ends in hand position. We obtained the first posture (posture 1) that fills in the
angles (q1-q7) in the jacobian from `Fundamental Principles of Tremor Propagation in the Upper Limb by Andrew D. Davidson and Steven K. Charles`. You can find the angles
specifically at `Fig. 1`. Change paths and names at your own risk. The files are interconnected and require upkeep if you change the name or path in one location. 
Check the `How to use` above so you have the functions in your matlab path.

These are the link-lengths used for the jacobian for a 50th percentile male (with height of 1.763 m):
- Length upper arm (Lua): 0.2853 m
- Length forearm (Lfa): 0.2723 m
- Length hand: 0.0873 m
- Distance from axis (wrist) to center of mass of hand (Lh): 0.069 m

TODO:
-----
* Check posture 1 (with Dr. Charles)
* Check work with Dr. Charles