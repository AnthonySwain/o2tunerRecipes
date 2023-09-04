'''
For testing custom voxel maps against a reference simulation
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
import pandas as pd

from o2tuner.system import run_command
from o2tuner.utils import annotate_trial
from o2tuner.optimise import optimise
from o2tuner.io import parse_json, dump_json, dump_yaml, parse_yaml, exists_file
from o2tuner.optimise import needs_cwd
from o2tuner.config import resolve_path

# Get environment variables we need to execute some cmds
O2_ROOT = environ.get("O2_ROOT") #Just the path to the O2 og enviroment
MCSTEPLOGGER_ROOT = environ.get("MCSTEPLOGGER_ROOT")

@needs_cwd
def Test_Custom_Map(trial,config):
    """
    Takes a custom voxel map and uses it against a reference simulation to see how good the mapping is. Re-uses old code so it does create a new
    CSV which can be seen as in-efficient but it works. 
    """

    #Imports functions
    absolute_filepath = "/home/answain/alice/o2tunerRecipes/ParticleSpecificMaps/optimise.py"
    optimise = imp.load_source("optimise", absolute_filepath)

    #Get neccessary information from the config file
    penalty_below = config["penalty_below"]

    #Reads the CSV file into a panda and then breaks it into 3 different seperate dataframes 
    csv_filepath = config['csv_filepath_hyper_fine']
    data = pd.read_csv(csv_filepath)

    csv_filepath_write = config['csv_filepath_write']

    Zmin_values = []
    Zmax_values = []
    innerR_values = []


    cylinder_counter = 0

    for index, row in data.iterrows():
  
        #if no previous cylinder
        #get the Z max of the cylinder
        Zmax = row['Zmax'] #get from dataframe
        Zmin = row['Zmin'] #get from dataframe
        R_data = row['radius']  #get from data frame

        
        Zmax_values.append(Zmax)
        Zmin_values.append(Zmin)
        innerR_values.append(R_data)

        #Writes the data to csv
        to_check = row['To_check']
        pdgs_list= row[5:].dropna().tolist()  # Drop NaN values and convert to list
        append_to_csv(Zmin,Zmax,R_data,to_check,pdgs_list, csv_filepath_write)

        cylinder_counter +=1
    

    #Run
    rel_steps_avg, rel_hits_avg = optimise.run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    optimise.annotate_trial(trial, "rel_steps", rel_steps_avg)
    optimise.annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return optimise.compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"],penalty_below)

def append_to_csv(Zmin, Zmax, radius, to_check, pdgs_list, csv_filepath):
    # Load existing CSV if it exists, otherwise create a new DataFrame
    try:
        existing_df = pd.read_csv(csv_filepath)
        #max_lengths = [len(pdgs_list), existing_df['pdgs_list'].apply(lambda x: len(eval(x))).max()]
        #max_pdgs_length = max(max_lengths, default=0)

    except FileNotFoundError:
        existing_df = pd.DataFrame()
        #max_pdgs_length = len(pdgs_list)

    # Determine the maximum length of pdgs_list among all rows
    

    # Fill pdgs_list with NaNs to match the maximum length
    #pdgs_list += [np.nan] * (max_pdgs_length - len(pdgs_list))

    # Create a DataFrame for the new row
    new_data = {
        'Zmin': [Zmin],
        'Zmax': [Zmax],
        'radius': [radius],
        'to_check': [to_check]
    }

    for i, pdg in enumerate(pdgs_list):
        new_data[f'pdg_{i+1}'] = pdg


    new_df = pd.DataFrame(new_data)
    # Append the new data to the existing DataFrame
    updated_df = pd.concat([existing_df, new_df], ignore_index=True)

    # Write the updated DataFrame to the CSV file
    updated_df.to_csv(csv_filepath, index=False)#,header=None)