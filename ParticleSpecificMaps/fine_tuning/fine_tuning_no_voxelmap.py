'''
This is pretty much the same as fine_tuning but it does not use voxel maps.
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
def fine_tuning_cylinders(trial,config):
    """
    Takes an optimal set of cylinders found so far (be it found using this algorithm of cylinder_xy for individual cylinders) from a csv file.
    Splits the cylinders into 3 groups: Main Barrel, Positive Z, Negative Z (positive/negative Z describe their position on the Z axis) (these groups will be indicated in the CSV file)

    Then for Positive Z / Negative Z, works from the outside in towards the main barrel, for each cylinder choosing Zmax/Zmin/radius with some leeway from 
    what was given in the CSV file (e.g. can pick from +-50 either side of the given values (50 is arbitary, will be different for Z vs radius))

    The outer Z value for a cylinder is defined by the previous's cylinder inner Z value to connect the cylinders. 

    Then we arrive at the inner barrel which will be defined by the radius given (will also be picked with some lee-way) as well as
    the inner Z values chosen from the detectors either side. 

    #Then passes a CSV of the cylinders onto the simulation
    """

    #Imports functions
    absolute_filepath = os.path.abspath(os.path.join(os.path.dirname(__file__), config["optimisation_framework_filepath"]))
    optimise = imp.load_source("optimise", absolute_filepath)

    #Get neccessary information from the config file
    penalty_below = config["penalty_below"]

    #Reads the CSV file into a panda and then breaks it into 3 different seperate dataframes 
    csv_filepath = config['csv_filepath_cylinders']
    data = pd.read_csv(csv_filepath)

    csv_filepath_write = config['csv_filepath_write']

    categories = ['positiveZ','negativeZ','mainBarrel']

    positiveZ_data = data[data['DetectorPart'] == categories[0]]
    negativeZ_data = data[data['DetectorPart'] == categories[1]]
    main_barrel_data = data[data['DetectorPart'] == categories[2]]

    positiveZ_data_sorted = positiveZ_data.sort_values(by='Zmax',ascending=False)
    negativeZ_data_sorted = negativeZ_data.sort_values(by='Zmin',ascending=True)

    #sorting the dataframes 

    #Then some code to sort in order of Zmax high -> low for positiveZ_data
    # and some code for Zmin low -> high for negativeZ_data

    
    leewayZ_axis_percent = config['ZaxisLeeway_percent'] #add something in config for this
    leewayRadial_percent = config['RadialLeeWay_percent'] #add something in config for this 

    Zmin_values_positiveZ = []
    Zmax_values_positiveZ = []
    innerR_values_positiveZ = []
    file_paths_positiveZ = []
    cylinder_counter_positiveZ = 0 # first cylinder = 0


    #Goes over from negative direction and positive direction, choosing parameters 
    for index, row in positiveZ_data_sorted.iterrows():
        if cylinder_counter_positiveZ == 0: 
            #if no previous cylinder
            #get the Z max of the cylinder
            Zmax_chosen = row['Zmax'] #get from dataframe
            Zmin = row['Zmin'] #get from dataframe
            R_data = row ['R'] #get from data frame

            Zmin_chosen = trial.suggest_float(f"minZ_positive{cylinder_counter_positiveZ}",Zmin*(1-leewayZ_axis_percent/100), Zmin*(1+leewayZ_axis_percent/100))
            Rchosen = trial.suggest_float(f"MaxR_positive{cylinder_counter_positiveZ}", R_data*(1 - leewayRadial_percent/100), R_data*(1+leewayRadial_percent/100))
           
        
        else:
            Zmin = row['Zmin'] #get from dataframe
            R_data = row['R'] #get from dataframe

            Zmax_chosen = Zmin_values_positiveZ[cylinder_counter_positiveZ -1]
            Zmin_chosen = trial.suggest_float(f"minZ_positive{cylinder_counter_positiveZ}",Zmin*(1-leewayZ_axis_percent/100), Zmin*(1+leewayZ_axis_percent/100))
            Rchosen = trial.suggest_float(f"MaxR_positive{cylinder_counter_positiveZ}"  ,R_data*(1 - leewayRadial_percent/100), R_data*(1+leewayRadial_percent/100))


        Zmax_values_positiveZ.append(Zmax_chosen)
        Zmin_values_positiveZ.append(Zmin_chosen)
        innerR_values_positiveZ.append(Rchosen)

    
        #Writes the data to csv
        to_check = row['To_check']
        pdgs_list= row[5:].dropna().tolist()  # Drop NaN values and convert to list
        append_to_csv(Zmin_chosen,Zmax_chosen,Rchosen,to_check,pdgs_list, csv_filepath_write)

        cylinder_counter_positiveZ +=1
    

    Zmin_values_negativeZ = []
    Zmax_values_negativeZ = []
    innerR_values_negativeZ = []
    file_paths_negativeZ = []
    cylinder_counter_negativeZ = 0

    #
    for index, row in negativeZ_data_sorted.iterrows():

        
        if cylinder_counter_negativeZ == 0: 
            #if no previous cylinder
            #get the Z max of the cylinder
            Zmin_chosen = row['Zmin'] #get from dataframe
            Zmax = row['Zmax'] #get from dataframe
            R_data = row['R'] #get from data frame

            Zmax_chosen = trial.suggest_float(f"minZ_negative{cylinder_counter_negativeZ}",Zmax*(1+leewayZ_axis_percent/100), Zmax*(1-leewayZ_axis_percent/100))
            Rchosen = trial.suggest_float(f"MaxR_negative{cylinder_counter_negativeZ}", R_data*(1 - leewayRadial_percent/100), R_data*(1+leewayRadial_percent/100))
           
        
        else:
            Zmax = row['Zmax'] #get from dataframe
            R_data = row['R'] #get from dataframe
            
            Zmin_chosen = Zmin_values_negativeZ[cylinder_counter_negativeZ -1]
            Zmax_chosen = trial.suggest_float(f"minZ_negative{cylinder_counter_negativeZ}",Zmax*(1+leewayZ_axis_percent/100), Zmax*(1-leewayZ_axis_percent/100))
            Rchosen = trial.suggest_float(f"MaxR_negative{cylinder_counter_negativeZ}", R_data*(1 - leewayRadial_percent/100), R_data*(1+leewayRadial_percent/100))
        
        Zmax_values_negativeZ.append(Zmax_chosen)
        Zmin_values_negativeZ.append(Zmin_chosen)
        innerR_values_negativeZ.append(Rchosen)
    
        #Writes the data to csv
        to_check = row['To_check']
        pdgs_list= row[5:].dropna().tolist()  # Drop NaN values and convert to list
        append_to_csv(Zmin_chosen,Zmax_chosen,Rchosen,to_check,pdgs_list, csv_filepath_write)

        cylinder_counter_negativeZ +=1
    

    Rmain_barrel = main_barrel_data.iloc[0]['R']
    to_check_main_barrel = main_barrel_data.iloc[0]['To_check']

    R_mainbarrel_chosen = trial.suggest_float("RmainBarrel", Rmain_barrel*(1- leewayRadial_percent/100), Rmain_barrel*(1+ leewayRadial_percent/100)) 
    Z_barrel_max = Zmin_values_positiveZ[-1]
    Z_barrel_min = Zmax_values_negativeZ[-1]

    pdgs_list= main_barrel_data.iloc[0][5:].dropna().tolist()  # Drop NaN values and convert to list
    append_to_csv(Z_barrel_min,Z_barrel_max,R_mainbarrel_chosen,to_check_main_barrel,pdgs_list, csv_filepath_write)

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