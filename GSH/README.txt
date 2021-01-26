README: Gravity Software TUDelft

NOTE: This software is not tested thoroughly for commercial use. If used in research please give credits to the developer and his University.

written by Bart Root, 06-nov-2014, Delft University of Technology:

The following software package contains MATLAB scripts to do Global Spherical  Harmonic Analyses and Synthesis for Crust1.0. The setup is build as follows:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
There are two executables

run.m

and 

plotResults_full.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run.m is divided in to parts:

1. GSHA -> performed by model_SH_analysis.m
2. GSHS -> performed by model_SH_synthesis.m

Both modules can be used separate. THE GSHA module uses a importfile = inputModel.m. Here the density model can be defined by using .gmt files. GMT files must have a certain format, which can be inspected in the given tutorial files (dir: Data). The unit in .gmt files are SI/1000. So, km or g/cm3.

The GSHS module will return a structure ‘data’. In this structure are all the different results situated. These results can be view with plotResults_full.m

The run.m file can be modified if different area or resolution is wanted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Structure:

Run.m uses: 

	- inputModel.m
	- model_SH_synthesis.m
		- gravityModule.m
			- visu2plm_ww.m
		- gravityModule_full.m (not tested yet, slow!!!!)
			- visu2plm_ww.m
	- model_SH_analysis.m
		- import_layer.m
			- gmt2matrix.m
		- layer_SH_analysis.m
			- gsha_crust.m
				-cs2sc.m
			- sc2vecml.m
			- geocradius.m (MATLAB func. toolbox aero)
					and only when geoid is WGS84
			- neumann.m (but is not used when wls is on)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Another useful script is the visual_gmtfile.m located in the Tools directory.
This function will read in a gmt file and plot the file. It allows you to check 
easily the file and if an error is present in the file it will not run, but then the 
GSH code will also not work. 

Another common error is a 'single' value NaN in the Gmt file, due to interpolation error.
This needs to be checked, because the GSH code will also stop working, when this is the case.

If you have any comments, feedback, or recommendations, please let me know: b.c.root@tudelft.nl

Thank you and enjoy!