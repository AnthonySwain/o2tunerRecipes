# Brief Overview of the optimisations for the O2pdgspecific build

## Prerequisites
The AliceO2 build (O2pdgspecific) with the modified stepping function can be found here: https://github.com/AnthonySwain/AliceO2/tree/ParticleSpecificMaps. The relevant readme is found in the "Steer" directory.

The macros found in this repo - https://github.com/AnthonySwain/AliceO2MacroDev - are used as a submodule. . The filepaths in the config files may (will) need changing to point towards the macros found in this repo.
You will need to run: 
git submodule init
git submodule update


Some optimisations make use of voxelmaps made using VecGeom - https://gitlab.cern.ch/VecGeom/VecGeom/-/blob/89a05d148cc708d4efc2e7b0eb6e2118d2610057/VecGeom/base/FlatVoxelHashMap.h#L207

This framework is built upon O2Tuner - https://github.com/AliceO2Group/o2tuner 

## Commands
o2tuner -w voxels_test -s cylinder_xy -c /home/answain/alice/o2tunerRecipes/ParticleSpecificMaps/cylinder_xy/config.yaml

This will run the cylinder_xy optimisation found in "cylinder_xy/" in a folder created in the CWD called "voxels_test".

## Note
Some optimisations use voxel maps, using VecGeom. Some simply work out the geometry in the stepping function and stop the transport the particles that way - it should be obvious which optimisations do which. 

## Optimisations

***ExampleCSVInputs*** contains examples of input CSV files that the optimisations below take as input. 

***zdc_skip*** is a variable in all config files - it indicates whether the ZDC detectors should be included in the simulations. 

***Events, Jobs & Trials*** should be given in the relevant config file.

***VoxelMap dimensions*** should be given in the relevant config file.

### cylinder_xy
Finds the optimal mapping for a cylinder given a certain region with Zmin,Zmax and maximum radius as constraints given.

### fine_tuning
Takes a set of cylinders from CSV input (hopefully found through cylidner_xy) and fits these together so that the Z axis is
filled with cylidners and then tries to optimise further to find the optimal Z extents and radii. 

### fine_tuning_fixed_Z
Same as fine_tuning but only the radius of the cylinders is tuned, not the Z extents. 
This also splits the cylinders given into more cylinders depending on the optimisation chosen (hyper_fine_splitting does this)
 I.e if 5 cylinders are given to define the cuts, this is turned into how many specified (i.e. 50) which are then optimised to perform a "hyper_fine" 
 tuning. The code isn't the most clean so to to perform the hyper_fine tuning careful consideration of the csv filepath is needed. 

### fine_tuning_no_voxel_map
Same as fine_tuning but does not use voxelmaps to speed up the optimisation process (creating voxel maps takes a while). Input is a CSV file which depics the Z extent and radii of cylndrical cuts.

### Test_Custom_Map
Tests csv of cylinders, or set of voxel maps, and runs the map against a reference simulation to test how good it is. 
To Note: 
Needs to specify in config whether it is a voxel map or csv of cylinders with the NoVoxelMap variable.
The csv of cylinders is called "cylinders_data.csv".


### AddZDC
Takes a mapping which was created without the ZDC detectors (because ZDC takes so much computational power) and adds a specified amount of cylinders
either side of the current mapping up to the ZDC detectors, then optimsies **_only_** these cylinders further. 

