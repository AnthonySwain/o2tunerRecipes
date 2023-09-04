'''
Optimisation for a single cylinder, Zmin/max and Rmax should be defined in the yaml file to help speed-up the process of finding the desired cylinder. 
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
def cylinder_xy(trial, config):
    """
    Works from outside in with a radial hashmap (cylindrical with the cylindrical axis corresponding to the beam axis)
    Everything outside the cylinder is a blackhole. 
    """

    #Imports functions
    absolute_filepath = "/home/answain/alice/o2tunerRecipes/ParticleSpecificMaps/optimise.py"
    optimise = imp.load_source("optimise", absolute_filepath)

    #Get neccessary information from the config file
    penalty_below = config["penalty_below"]
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    save_root_hashmap_file = config["hashmap_file"]

    #Min/Max parameters for the cylinder
    Rmax = config["Rmax"]
    Zmax = config["Zmax"]
    Zmin = config["Zmin"]

    #Chooses the parameters and creates the radial hashmap. 

    Zchoice1 = trial.suggest_float("Zchoice1",Zmin,Zmax)
    Zchoice2 = trial.suggest_float("Zchoice2",Zmin,Zmax)

    if Zchoice1 > Zchoice2:
        minZchosen = Zchoice2
        maxZchosen = Zchoice1

    else:
        minZchosen = Zchoice1
        maxZchosen = Zchoice2
            
    Rchosen = trial.suggest_float("MaxR", 0,Rmax)

    optimise.CreateRadialHashMap(trial, config["CreateRadialHashMapFullPath"], nx, ny, nz, save_root_hashmap_file,Rchosen,Zextent=[minZchosen,maxZchosen])

    optimise.add_map_to_txt(config["txt_of_maps"],config["hashmap_file"],config["particle_map"])
    # rng = np.random.default_rng()
    # batch_id = rng.integers(0, batches)

    #Run
    rel_steps_avg, rel_hits_avg = optimise.run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    annotate_trial(trial, "rel_steps", rel_steps_avg)
    annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return optimise.compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"],penalty_below)
