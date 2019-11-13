# Analysis_CLAS6

This file is meant as a guide through how to understand this analysis framework for double charged pion electroproduction off the proton. The purpose of this analysis is to extract the polarization observables as well as the electron beam helicity asymmetry for this reaction. 

There is another repository called *analysis* on this GitHub, but please feel free to ignore that one. While it accomplishes goals up until raw yields, it is by no means a complete analysis and is significantly clunkier. 

## How to Navigate

There are several analysis programs in here
1. analysis_1 
  1. This program is a first pass through the "raw" rootfile. It performs event selection and does all cuts relevant to that process
  1. Its source code is located in "/analysis_clas6/src"
  1. Its outputs are located in "/analysis_clas6/bin"
1. analysis_2
  1. This program is the second pass through the TTrees created in analysis_1
  1. This program extracts the four vectors from the event selected TTree and outputs the final cross sections
  1. Its source code is located in "/analysis_clas6/analysis_2/src"
  1. Its output is located in "/analysis_clas6/analysis_2/bin"
1. Golden Run Determination
  1. This program was used to determine the list of Golden Runs for both e1-6 and e1f
  1. It takes in the individual root files from the designated run, finds each run's average current, and provides a list of runs that have currents which fall within the region of interest. It also outputs plots showing each run's average current
  1. Its source code is located in "analysis_clas6/Golden\Run\Determination\/src"
  1. Its output is located in "analysis_clas6/Golden\Run\Determination\/bin"
1. parameter_fitting
  1. All the fitting for this analysis was performed in python, and then exported to the other C++ files for implementation due to python's superior fitting algorithms in comparison to C++ and R00t. 
  1. It takes in a variety of things, but *at present 11/13/19* has yet to be built 
1. Simulation
  1. Simulation files are all located in "analysis_clas6/sim/"
  1. The files here are the ones used to perform simulations for both e1-6 and e1f data sets
