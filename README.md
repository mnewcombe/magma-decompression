# magma-decompression
MATLAB code for using water concentration gradients in olivine or cpx to constrain syneruptive magma decompression rates. This example is set up to fit a water concentration gradient measured along the crystallographic 'a' axis of an olivine phenocryst from the 1977 eruption of Seguam volcano. Please cite Newcombe et al. 2020, in revision at Journal of Volcanology and Geothermal Research.

Instructions for running H-in-olivine Monte Carlo error analysis:

%{ 
input files: 
seguam_ol1.txt is my data file. The first column is radial distance in 
microns, the second column is water concentration in ppm.

Solex_PHrelation_Seguam.txt contains the degassing path I used for Seguam
(calculated using VolatileCalc or Solex). The first column is pressure in bars, the
second column is water concentration in wt%. 

master_script is the file where you should define all setup parameters

olivineMC.m is the function that performs the diffusion calculation and
least-squares error calculation

olivineMCplot.m is the same as olivineMC.m except that it produces a plot
at the end
