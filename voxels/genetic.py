''' Genetic Algorithm to find the best voxelmaps of a generation and breed them. '''

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


#from o2tuner.system import run_command
#from o2tuner.config import resolve_path

@needs_cwd #Ran in its on directory (cwd = current working directory)
def genetic(trial,config):
    '''
    Genetic algorithm implementation for the optimisation
    The dependency 'next_gen' has already kept, bred, and, mutated the maps (or created the first generation)
    This just needs to find the .txt files containing the maps and create the .root hashmap and run the simulation
    '''

    #Imports functions
    absolute_filepath = "/home/answain/alice/o2tunerRecipes/voxels/optimise.py"
    optimise = imp.load_source("optimise", absolute_filepath)

    #Gets what's needed from the config file
    penalty_below = config["penalty_below"]
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    save_root_hashmap_file = config["hashmap_file"]
    move_txt_to = config["voxels_sampled_file"]

    #Get trial numb
    trial_numb = trial.number
    print(trial_numb) 

    '''Gets the hashmap for this trial from the folder of hashmaps that 'next_gen' dependency created
    and sticks it in the cwd and then creates the .root voxelmap from it'''
    hash_map_txt_filepath = os.path.join(os.path.dirname(os.path.dirname(os.getcwd())),f"next_gen/HashMaps/hashmap_{trial_numb}.txt")
    optimise.relocate_file(hash_map_txt_filepath,f"{os.getcwd()}/{move_txt_to}")
    optimise.create_hash_map(config["CreateHashMapFromTxtMacroFullPath"],move_txt_to,nx,ny,nz,save_root_hashmap_file)

    #Run
    rel_steps_avg, rel_hits_avg = optimise.run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    optimise.annotate_trial(trial, "rel_steps", rel_steps_avg)
    optimise.annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return optimise.compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"],penalty_below)


def saveastxt(filename, input_list):
    '''Saves a mapping to a .txt file'''
    concatenated_string = ''.join(map(str, input_list))
    try:
        with open(filename, 'w') as file:

            file.write(concatenated_string)

    except Exception as e:
        print("Error while saving the file:", e)

def mutate(mutation_rate,voxel):
    '''
    mutates the mapping
    mutation_rate = mean %chance of a value flip for each element
    voxel = list of values'''
    mean = mutation_rate/100
    num_samples = 1
    how_many_mutations = math.ceil(np.random.exponential(mean*len(voxel), num_samples))

    which_elements_mutated = np.random.randint(0, len(voxel), size=how_many_mutations)
    
    #Flip the voxel value!
    for i in which_elements_mutated:
            #bitflip! 
            if voxel[i]:
                voxel[i] = '0'
            else: 
                voxel[i] = '1'

    return voxel

def breed(mutation_rate, voxel_filepath1, voxel_filepath2, save_loc,parentnumb1,parentnumb2):
    '''Reads 2 voxelfilepaths given, breeds them and saves them to a place.
    Can call mutate() during the function to add mutations. 
    mutation_rate = average %chance of a single element mutating
    voxel_filepath1/2 = filepaths of the parents
    save_loc = where to save the children
    parentnumb1/2 = just some unique values to identify the children / give them a unique name
    'oh childof0120.txt have you grown!' 
    '''

    #Read each voxelmap and make them a np.array
    with open(voxel_filepath1, 'r') as voxelmap1:
        voxel1 = voxelmap1.read().strip()

    voxel1= np.array(list(voxel1), dtype=int)
    
    with open(voxel_filepath2, 'r') as voxelmap2:
        voxel2 = voxelmap2.read().strip()

    voxel2= np.array(list(voxel2), dtype=int)


    #the maps should be the same length and this will be assumed. 
    length1 = len(voxel1)
    length2 = len(voxel2)

    #Breed
    crossover_point = random.randint(1, length1 - 1)
    child1 = np.append(voxel1[:crossover_point],voxel2[crossover_point:])
    child2 = np.append(voxel2[:crossover_point], voxel1[crossover_point:])

    #Mutate!
    child1 = mutate(mutation_rate,child1)
    child2 = mutate(mutation_rate,child2)

    #Where to save the children
    cwd = os.getcwd()

    save_path1 = join(cwd,save_loc,f"child1of{str(parentnumb1)}{str(parentnumb2)}.txt")
    save_path2 = join(cwd,save_loc,f"child2of{str(parentnumb1)}{str(parentnumb2)}.txt")
   
    #Saves the children (our hero!)
    np.savetxt(save_path1, child1, fmt='%d', delimiter='', newline='')
    np.savetxt(save_path2, child2, fmt='%d', delimiter='', newline='')

    #Returns the filepaths of the children
    return [save_path1,save_path2]

def GetBestLossFunctions(FilePathToDB, PathToTrials,keep_percent): 
    '''
    Open a DB file, finds the best loss functions of a generation and gets the corresponding
    filepaths to the voxels.txt file for each trial and returns these as a list

    FilePathToDB = says what it does on the tin
    PathToTrials = path to where the trials of the previous generation were ran
    keep_percent - the percentage of trials to keep and additionaly breed. (i.e the top x% is kept, but they are also bred to get some kids) 
    '''
    
    '''Some SQL to get a list of filepaths of the best % of loss functions'''
    voxel_filepaths = []
    db = sqlite3.connect(FilePathToDB) #Open DB
    cursor = db.cursor()

    table = "trial_values"
    query = f'SELECT COUNT(*) FROM {table}'
        
    cursor.execute(query)
    no_rows = cursor.fetchone()[0]
    
    value_column_name = 'value'
    trial_id_column = 'trial_id'

    query = f'SELECT {value_column_name}, {trial_id_column} FROM {table} ORDER BY {value_column_name} ASC LIMIT {no_rows * keep_percent}'

    cursor.execute(query)
    result = cursor.fetchall()
    
    loss_values = []
    trial_ids = []

    for row in result:
        loss_values.append(row[0])
        trial_ids.append(row[1])
    

    folder_column = 'value_json'
    trial_user_attributes = "trial_user_attributes"
    
    query = f'SELECT {folder_column} FROM {trial_user_attributes} WHERE trial_id = ? AND key = ?'
    
    for i in trial_ids:
        cursor.execute(query,(i,'cwd'))
        folder = (cursor.fetchall()[0][0])
        folder = folder.replace('"','')
        folder = folder.replace('.','')
        folder = folder.replace('/','')
        
        path = join(PathToTrials,folder,"voxels.txt")
        voxel_filepaths.append(path)

    db.close() #Close DB

    return voxel_filepaths

def relocate_files(filepaths, destination_folder):
    '''
    Relocates (by copy & paste) files to a new location
    filepaths = list of filepaths we want to relocate
    destination_folder = where to locate them to
    '''
    for i in range(len(filepaths)):
        if filepaths[i].endswith(".txt"):
            try:
                new_filepath = os.path.join(os.getcwd(),destination_folder, f"hashmap_{i}.txt")
                shutil.copy(filepaths[i], new_filepath)
                print(f"Copied {filepaths[i]} to {new_filepath}")
            except Exception as e:
                print(f"Error copying {filepaths[i]}: {e}")

def next_gen(inspectors, config):
    '''
    Calculates and saves the next generation of hashmaps for the genetic algorithm
    '''

    #Gets required info from .yaml file
    #FilePathToDB = config['DBFilepath']
    keep_percent = config['KeepPercent']
    mutation_rate = config['mutation_rate_percentage']
    save_loc_breed = config['saveBredMaps']
    save_loc_all = config['saveNextGen'] #some way of working out the generation number! 
    population = config['population']
    fill_percent = config['fill_percent']
    nx = config['n_voxels_x']
    ny = config['n_voxels_y']
    nz = config['n_voxels_z']


    #Checks if there was a previous generation 
    current_generation = current_gen()
    if (current_generation == None) or (current_generation == 0):
        previous_generation = 0
        current_generation = 0
    else:
        previous_generation = current_generation - 1 

    #Creates a directory for the next gen of hashmaps 
    cwd = os.getcwd()
    if (os.path.exists(cwd+save_loc_all)==False):
        os.mkdir(cwd+"/"+save_loc_all)

    #If there was a previous generation, make the next one
    if (current_generation > 0):
        print(f"Current generation = {current_generation}. Creating the next_gen of hashmaps from generation {previous_generation}.")
        
        #Creates directory for the children
        if (os.path.exists(cwd+save_loc_breed)==False):
            os.mkdir(cwd+"/"+save_loc_breed)

        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
        prev_gen = join(base_dir,f"generation_{previous_generation}")
        PathToPrevTrials = join(prev_gen,"voxels_test","genetic")
        FilePathToDB = join(PathToPrevTrials,"genetic.db")
        
        #Gets the best loss functions from the prev generation
        elite_filepaths = GetBestLossFunctions(FilePathToDB,PathToPrevTrials,keep_percent)

        #Breed the best! (breeds all the best filepaths with eachother)
        bred_filepaths = []
        for i in range(len(elite_filepaths)):
            for j in range(len(elite_filepaths)):
                if (j <= i):
                    continue
                print(f"Breeding {i} with {j}.")
                bred_filepaths += breed(mutation_rate, elite_filepaths[i],elite_filepaths[j],save_loc_breed,i,j)
                
        #Gets all the filepaths of the new hashmap.txt files in the same list
        new_voxel_map_filepaths = elite_filepaths + bred_filepaths

        #Relocates them all to the same single folder. 
        #Should just put the children in this place in the first place to prevent unneccessary copying / savetime
        relocate_files(new_voxel_map_filepaths, save_loc_all)

        #Also should put a function which copy pastes the reference and base into the new generational folder:) 

    #Create a new generation if there wasn't a previous one. 
    else:
        print("This is the very first generation, intialising.")
        initialise_generation(population,nx*ny*nz,fill_percent, save_loc_all)

    return(True)
    
def current_gen():
    '''
    Scans one above the current working directory for generation_x folders, 
    the one with the highest value x is taken as the previous geneartion 
    (returns 0 if there was none -> means to initialise a new gen)
    '''

    target_prefix = "generation"
    parent_directory = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd()))) #don't judge me

    matching_folders = find_folders_with_prefix(parent_directory, target_prefix)
    
    current_gen = extract_highest_number(matching_folders)
    return current_gen 

def extract_highest_number(strings):
    '''
    Finds the highest number at the end of a string in a list of strings
    '''
    highest_number = None

    for string in strings:
        match = re.search(r'\d+$', string)
        if match:
            number = int(match.group())
            if highest_number is None or number > highest_number:
                highest_number = number

    return highest_number

def find_folders_with_prefix(parent_folder, prefix):
    '''
    Finds folders that start with 'prefix'
    '''
    matching_folders = []

    for folder_name in os.listdir(parent_folder):
        folder_path = os.path.join(parent_folder, folder_name)
        if os.path.isdir(folder_path) and folder_name.startswith(prefix):
            matching_folders.append(folder_path)

    return matching_folders

def initialise_generation(population, length, fill_percent, folder):
    '''
    Initialises a new generation of maps

    population = how many maps to create
    length = length of the binary string
    fill = % of map that will be 1's 
    '''

    #Better to work in absoloute filepaths 
    #folder = "/HashMaps"
    folder = os.getcwd()+folder
    if (os.path.exists(folder)==False):
        os.mkdir(folder)

    for i in range(population):
        array = generate_array_with_percentage(length,fill_percent)
        filepath = f"{folder}/hashmap_{i}.txt"
        saveastxt(filepath,array)

def generate_array_with_percentage(size, percentage):
    '''
    Generates an array of 0/1s with the average number of 1s total being size*percentage 
    '''
    fill = np.random.exponential(size*percentage/100, 1)
    num_ones = int(fill)
    num_zeros = size - num_ones

    ones_array = np.ones(num_ones, dtype=int)
    zeros_array = np.zeros(num_zeros, dtype=int)

    rng = np.random.default_rng()
    

    combined_array = np.concatenate((ones_array, zeros_array))
    rng.shuffle(combined_array)

    return combined_array


