This code is for the simulated populations.

It requires the python module SIMUPOP.

When run, it will prompt for a set of simulation parameters and a run name.

It will produce a set of replicate samples from simulated populations in a folder containing a number of files:


<one level up from folder containing SimulateDemoNe.py>/Samples/<run_name>/sampleX.dat

	For X in 0 to N-1, these are fstat formated genotypes.


<one level up from  folder containing SimulateDemoNe.py>/Samples/<run_name>/param.csv
	This file retains the parameters used in the simulation.


<one level up from folder containing SimulateDemoNe.py>/Samples/<run_name>/runInfo.csv
	This file saves specific information relevant to each replicate population. 
	In particular, the demographic Ne for the 4 generations prior to the sample generation.
