
'''
Uses a bayesian type approach to find try and find an optimal voxel map for the inner barrel
'''


import sqlite3
import sys
from os.path import join, abspath
from os import environ
import random
import os
import numpy as np 
import shutil
import re
import math
import imp

from o2tuner.system import run_command
from o2tuner.utils import annotate_trial
from o2tuner.optimise import optimise
from o2tuner.io import parse_json, dump_json, dump_yaml, parse_yaml, exists_file
from o2tuner.optimise import needs_cwd
from o2tuner.config import resolve_path


# Get environment variables we need to execute some cmds
O2_ROOT = environ.get("O2_ROOT") #Just the path to the O2 og enviroment
MCSTEPLOGGER_ROOT = environ.get("MCSTEPLOGGER_ROOT")


@needs_cwd #Ran in its on directory (cwd = current working directory)
def objective(trial, config):
    """
    The central objective function for the optimisation
    """

    # construct voxel hash map
    # 1. how many voxels in x, y, z
    #    could use the same logic as here: https://gitlab.cern.ch/bvolkel/VecGeom/-/blob/master/VecGeom/base/FlatVoxelHashMap.h#L163 to only have one 1D list
    #    e.g. write true/false (or 0, 1) to a file which will then be read by ROOT to make the actual HashMap and store it in a ROOT file
    # 2. which to switch on

    #Imports functions
    absolute_filepath = "/home/answain/alice/o2tunerRecipes/voxels/optimise.py"
    optimise = imp.load_source("optimise", absolute_filepath)

    #Gets whats needed from teh config file
    penalty_below = config["penalty_below"]
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    
    min = config["min"]
    length = config["length"]
    map_creation_macro = config["map_creation_macro_fullpath"]
    
    save_file_line_by_line = config["voxels_sampled_file"]
    save_root_hashmap_file = config["hashmap_file"]

    #Creates the hashmap (first by writing 0,1s to a .txt file (sample_voxels) then reading this .txt file (create_hash_map))
    optimise.sample_voxels(trial, nx * ny * nz, save_file_line_by_line,map_creation_macro,
                           save_root_hashmap_file,min,length,nx,ny,nz)
    #optimise.create_hash_map(config["CreateHashMapFromTxtMacroFullPath"], save_file_line_by_line, nx, ny, nz, save_root_hashmap_file)

    #Run it 
    rel_steps_avg, rel_hits_avg = optimise.run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    annotate_trial(trial, "rel_steps", rel_steps_avg)
    annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return optimise.compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"], penalty_below)
