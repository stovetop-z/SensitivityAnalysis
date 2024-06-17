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

Repeat the above steps for the `RVC2` folder with the modified step #3 to `Selected Folders and Subfolders`.

The mat files are already integrated. No changes needed.

To perform a sensitivity analysis of the tremor propagation model including the 4th submodel (translation of joint displacement to hand position), first run `create_robot_arm.m`. Then, go to `AIE_X_Sensitivity.m` 
and `SIE_X_Sensitivity.m`. These take about 5 minutes to run. To plot the data, use `plot_x.mlx`. 

Files:
------
The main files to be aware of:
* create_robot_arm - We utilize R7 robot arm to create jacobian for the 4th submodel
* AIE_X_Sensitivity - This code runs the all input excitation for the 4th and last submodel
* AIE_Sensitivity - This code runs the aie for the 3rd submodel (joint displacement)
* Parameters - Holds important parameter information and some useful variables for global use
* plot_joint_displacement - INCOMPLETE PROJECT to plot the 3rd submodel
* plot_x - Plots the 4th submodel
* SIE_Sensitivity - This code runs the single input excitation for the 3rd submodel
* SIE_X_Sensitivity - This code runs the sie for the 4th submodel

Information for the ignorant:
-----------------------------
The `X` indicates the vector x, which represents the 3 dimensions in space. If there is no `X`, then we are using the tremor model that ends in joint displacement `Q`. 
If the `X` is present in the name of a file, then we are using the tremor model that ends in hand position. We obtained the first posture (posture 1) that fills in the
angles (q1-q7) in the jacobian from `Fundamental Principles of Tremor Propagation in the Upper Limb by Andrew D. Davidson and Steven K. Charles`. You can find the angles
specifically at `Fig. 1`. Change paths and names at your own risk. The files are interconnected and require upkeep if you change the name or path in one location. 
Check the `How to use` above so you have the functions in your matlab path.

These are the link-lengths used for the jacobian for a 50th percentile male (with height of 1.763 m):
- Length upper arm (Lua): 0.2853 m
- Length forearm (Lfa): 0.2723 m
- Distance from axis (wrist) to center of mass of hand (Lh): 0.069 m

Issues:
-------
Paul's and mine sensitivities do not seem to be aligning. By checking the transfer function with all the default parameters in place, we see the follow differences
between mine and his (Gs/Gp):
    1.9955   16.1657    1.1338    4.5483    0.8395    0.8245    1.0101    1.0000    1.0000    1.0461    1.7599    0.4812    0.4747    0.6915    0.3794
    1.0000    1.0000    1.0000    1.0000    2.3142    0.6132    1.0000    1.0000    1.0000    1.0577   15.4800    0.0933    0.1042    0.1050    0.0773
    0.2214    0.4847    0.4197    0.2388    2.1308    0.2659    1.0702    1.0000    1.0000    3.0161   10.0433   18.7752   10.0067   10.0311    7.4156
    0.5718    0.6889    0.5378    0.6016    1.0439    1.0329    1.0109    1.0000    1.0000    0.9852    0.6342    0.7955    1.1881    0.8502    0.9230
    1.0000    1.0000    1.0000    1.0000   14.3922   40.3710    1.0000    1.0000    1.0000    2.7243   38.4368    1.8237    1.2763    1.8302    1.0472
    1.0000    1.0000    1.0000    1.0000   13.3787    0.4414    1.0000    1.0000    1.0000    2.7311    1.0887    0.9883    1.0847    1.0814    0.8244
    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000

Below is my transfer function:
    0.0020    0.0056    0.0003    0.0046    0.0590    0.0350    0.0388    0.0357    0.0679    0.0556    0.0228    0.0175    0.0182    0.0344    0.0102
    0.0568    0.1669    0.0091    0.1327    0.0655    0.0209    0.0404    0.0025    0.0048    0.0042    0.0294    0.0626    0.0791    0.0759    0.0110
    0.0387    0.1937    0.0132    0.0871    0.0909    0.0253    0.0766    0.0022    0.0041    0.0114    0.0337    0.3596    0.4437    0.4170    0.0648
    0.0041    0.0064    0.0008    0.0073    0.2083    0.1237    0.1076    0.0994    0.1889    0.1454    0.0260    0.0775    0.2092    0.2218    0.0959
    0.0832    0.2181    0.0147    0.1823    2.3547    1.3568    0.1305    0.0799    0.1518    0.3356    1.7633    0.2812    0.4002    0.4313    0.1895
    0.1893    0.5073    0.0676    0.1171    0.3200    0.0674    0.2503    0.0282    0.0536    0.0616    0.1458    1.9773    2.4654    2.3109    0.3523
    0.0147    0.0280    0.0033    0.0261    0.3582    0.2221    0.2922    0.2783    0.5289    0.4637    0.3244    0.2437    1.4341    1.1388    0.5777

Below is Paul's
    0.0010    0.0003    0.0003    0.0010    0.0703    0.0424    0.0384    0.0357    0.0679    0.0531    0.0130    0.0363    0.0383    0.0497    0.0270
    0.0568    0.1669    0.0091    0.1327    0.0283    0.0341    0.0404    0.0025    0.0048    0.0040    0.0019    0.6705    0.7597    0.7225    0.1429
    0.1747    0.3995    0.0314    0.3649    0.0427    0.0952    0.0715    0.0022    0.0041    0.0038    0.0034    0.0192    0.0443    0.0416    0.0087
    0.0072    0.0093    0.0016    0.0121    0.1995    0.1198    0.1064    0.0994    0.1889    0.1476    0.0411    0.0974    0.1761    0.2609    0.1039
    0.0832    0.2181    0.0147    0.1823    0.1636    0.0336    0.1305    0.0799    0.1518    0.1232    0.0459    0.1542    0.3136    0.2357    0.1810
    0.1893    0.5073    0.0676    0.1171    0.0239    0.1527    0.2503    0.0282    0.0536    0.0225    0.1340    2.0006    2.2729    2.1369    0.4273
    0.0147    0.0280    0.0033    0.0261    0.3582    0.2221    0.2922    0.2783    0.5289    0.4637    0.3244    0.2437    1.4341    1.1388    0.5777

The differences are striking.

Even after changing Paul's code slightly to calculate the transfer function with his parameter's and it yielded the same results as mine (as expected). This leads
me to believe there are different operations done for numerical vs symbolic operations in matlab, possibly having to do with taking the inverse of the matrices,
which is already a complicated algorithm to do numerically.