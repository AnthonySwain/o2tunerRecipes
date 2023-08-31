'''
Fine tuning a set of cylinders i.e give their Z extents and inner radii and this will try to optimise further. 
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

    #Reads the CSV file into a panda and then breaks it into 3 different seperate dataframes 
    csv_filepath = config['csv_filepath']
    data = pd.read_csv(csv_filepath)

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
        filepath = f"positiveZmap{cylinder_counter_positiveZ}.root"
        file_paths_positiveZ.append(filepath)


        #then go create a cylinder, queen
        optimise.CreateRadialHashMap(trial, config["CreateRadialHashMapFullPath"], nx, ny, nz, filepath,Rchosen,Zextent=[Zmin_chosen,Zmax_chosen])
        
        
        values = row[4:].dropna().tolist()  # Drop NaN values and convert to list
        optimise.add_map_to_txt(config["txt_of_maps"],filepath,values)
        #function here which adds the mapfilepath and particles to delete to a .txt file

        cylinder_counter_positiveZ +=1
    
    Zmin_values_negativeZ = []
    Zmax_values_negativeZ = []
    innerR_values_negativeZ = []
    file_paths_negativeZ = []
    cylinder_counter_negativeZ = 0

    #check this to make sure all the positive and negatives are the right way round.
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
        filepath = f"negativeZmap{cylinder_counter_negativeZ}.root"
        file_paths_negativeZ.append(filepath)


        #then go create a cylinder, queen
        optimise.CreateRadialHashMap(trial, config["CreateRadialHashMapFullPath"], nx, ny, nz, filepath,Rchosen,Zextent=[Zmin_chosen,Zmax_chosen])

        #function here which adds the mapfilepath and particles to delete to a .txt file
        values = row[4:].dropna().tolist()  # Drop NaN values and convert to list

        optimise.add_map_to_txt(config["txt_of_maps"],filepath,values)

        cylinder_counter_negativeZ +=1
    

  

    Rmain_barrel = main_barrel_data.iloc[0]['R']
    R_mainbarrel_chosen = trial.suggest_float("RmainBarrel", Rmain_barrel*(1- leewayRadial_percent/100), Rmain_barrel*(1+ leewayRadial_percent/100)) 
    main_barrel_filepath = ["MainBarrelMap.root"]

    Z_barrel_max = Zmin_values_positiveZ[-1]
    Z_barrel_min = Zmax_values_negativeZ[-1]

    optimise.CreateRadialHashMap(trial, config["CreateRadialHashMapFullPath"], nx, ny, nz, main_barrel_filepath[0],R_mainbarrel_chosen,Zextent=[Z_barrel_min,Z_barrel_max])
    #function here which adds the mapfilepath and particles to delete to a .txt file
    values = main_barrel_data.iloc[0][4:].dropna().tolist()  # Drop NaN values and convert to list

    optimise.add_map_to_txt(config["txt_of_maps"],main_barrel_filepath[0],values)

    #now we have everything, add all the hashmaps
    #maybe put this as a function or something
    AddMapsMacroPath = config['AddMapsMacroFilePath'] #add it here
    maps_filepaths = main_barrel_filepath + file_paths_negativeZ + file_paths_positiveZ

    #add_all_maps(AddMapsMacroPath, maps_filepaths,save_root_hashmap_file,nx,ny,nz) #Not needed with new implementation

    #Run
    rel_steps_avg, rel_hits_avg = optimise.run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    optimise.annotate_trial(trial, "rel_steps", rel_steps_avg)
    optimise.annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return optimise.compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"],penalty_below)

def add_all_maps(AddMapsMacroPath,map_filepaths,final_save_loc,Nx,Ny,Nz):
    #goes through all of the hashmaps and adds them together. 

    for i in range(len(map_filepaths)-1):
 
       # add the first 2 maps, then add the third onto those, then 4th ect...
        if i == 0:
            MapPath1 = map_filepaths[i]
            MapPath2 = map_filepaths[i + 1]
        else:
            MapPath1 = MapSaveLoc
            MapPath2 = map_filepaths[i + 1]

        if i == (len(map_filepaths)-2):
            MapSaveLoc = final_save_loc
            print(MapSaveLoc)
        
        else:
            MapSaveLoc = f"AddedHashMaps{i}.root" 
        

        
        add_maps = f"root -l -b -q '{AddMapsMacroPath}(\"{MapPath1}\",\"{MapPath2}\",{Nx},{Ny},{Nz},\"{MapSaveLoc}\")'"

        #Runs the macro
        _, hashmap_file = run_command(add_maps, log_file="hits.dat")

 