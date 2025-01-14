Sensitivity Analysis
====================

Experiment to understand the sensitivity of the 120 parameters in the tremor propagation model. 
Initially created by Paul Curtis, expounded and added upon by Steven Zinn. The additions to Paul Curtis' work include:
* Numerical differentiation implementation rather than using symbolic differentiation
* Addition of a 4th submodel with jacobian through the RVC7 toolbox to translate joint displacement to hand velocity in the xyz. This added 3 more parameters [length of hand, length of forearm, length of upper arm]
* Plots to analyze sensitivities of 120 parameters

How to use:
-----------
After downloading this project, you will notice there are two subfolders: `functions` and `mats`. In the functions subfolder there are various 
function files that are used extensively in the main code for plotting and performing the sensitivity of the tremor propagation model. To use 
the functions, it is imperative that you ensure that they are in the matlab project's path. They should be added in the `Parameters.m`. The following steps will explain 
how to run the code to perform the analysis:

1. Go to `main.m`
2. If you do not have the jacobian already, uncomment line 23 (where `create_robot_arm` is located)
3. Under `Step 2` in `main.m`, change the postures you want to perform an analysis on if desired.
4. Run the file. It should take between 30 mins to an hour, depending on your computer system for the analysis to run to completion
5. Go to `AIE_Plot.m`.
6. Run the code.

And your DONE!

Files:
------
The main files to be aware of:
* `main.m`
* `create_robot_arm.m` - We utilize R7 robot arm to create jacobian for the 4th submodel
* `functions/AIE_Sensitivity_Analysis.m` - This code runs the all input excitation for the last submodel
* `Parameters.m` - Holds important parameter information and some useful variables for global use
* `AIE_Plot.m` - Plots the Sensitivity Analysis for AIE

Information for the ignorant:
-----------------------------


Issues:
-------
No computational issues.