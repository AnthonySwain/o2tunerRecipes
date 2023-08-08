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


#from o2tuner.system import run_command
#from o2tuner.config import resolve_path

def saveastxt(filename, input_list):
    concatenated_string = ''.join(map(str, input_list))
    try:
        with open(filename, 'w') as file:

            file.write(concatenated_string)

    except Exception as e:
        print("Error while saving the file:", e)

def mutate(mutation_rate,voxel):
    '''mutation_rate = %chance of a value flip for each element'''
    mean = mutation_rate
    std_dev = mutation_rate
    num_samples = 1
    how_many_mutations = np.random.Generator.normal(mean, std_dev, num_samples)
    
    which_elements_mutated = np.random.randint(0, len(voxel), size=how_many_mutations)
    
    #Flip the voxel value!
    for i in which_elements_mutated:
            if voxel[i]:
                voxel[i] = False
            else: 
                voxel[i] = True

    return voxel

def breed(mutation_rate, voxel_filepath1, voxel_filepath2, save_loc,parentnumb1,parentnumb2):
    '''Reads 2 voxelfilepaths given, breeds them and saves them to a place.
    Can call mutate() during the function to add mutations. '''

    #Read each voxelmap. 
    with open(voxel_filepath1, 'r') as voxelmap1:
        voxel1 = voxelmap1.read().strip()

    with open(voxel_filepath2, 'r') as voxelmap2:
        voxel2 = voxelmap1.read().strip()

    #these should be the same length and this will be assumed. 
    length1 = len(voxel1)
    length2 = len(voxel2)

    crossover_point = random.randint(1, length1 - 1)
    child1 = voxel1[:crossover_point] + voxel2[crossover_point:]
    child2 = voxel2[:crossover_point] + voxel1[crossover_point:]

    #Mutate!
    child1 = mutate(mutation_rate,child1)
    child2 = mutate(mutation_rate,child2)

    #Where to save the children
    save_path1 = save_loc + "/child1of" + str(parentnumb1) + str(parentnumb2) + ".txt" 
    save_path2 = save_loc + "/child2of" + str(parentnumb1) + str(parentnumb2) + ".txt"

    #Saves the children
    saveastxt(child1,save_path1)
    saveastxt(child2,save_path2)

    #Returns the filepaths of the children
    return [save_path1,save_path2]

def GetBestLossFunctions(FilePathToDB, keep_percent): #keep_percent - the percentage of trials to keep and breed. 
    '''Open a DB file, find the best loss functions of a generation and get the corresponding
    filepaths to the voxels.txt file for each trial'''

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
        voxel_filepaths.append(folder+"/voxels.txt")

    db.close() #Close DB

    return voxel_filepaths

def relocate_files(filepaths, destination_folder):
    '''Relocates (by copy & paste) files to a new location'''
    for i in range(len(filepaths)):
        if filepaths[i].endswith(".txt"):
            try:
                filename = os.path.basename(filepaths[i])
                new_filepath = os.path.join(destination_folder, f"hashmap{i}.txt")
                shutil.copy(filepaths[i], new_filepath)
                print(f"Copied {filepaths[i]} to {new_filepath}")
            except Exception as e:
                print(f"Error copying {filepaths[i]}: {e}")

def next_gen(config):
    '''Calculates and saves the next generation of hashmaps'''

    #Insert config into above from yaml file
    #FilePathToDB = config['DBFilepath']
    keep_percent = config['KeepPercent']
    mutation_rate = config['mutation_rate']
    save_loc_breed = config['saveBredMaps']
    save_loc_all = config['saveNextGen'] #some way of working out the generation number! 
    save_loc_breed = "/children"
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

    if (os.path.exists(save_loc_breed)==False):
        os.mkdir(save_loc_breed)

    if (os.path.exists(save_loc_all)==False):
        os.mkdir(save_loc_all)

    #If there was a previous generation, make the next one
    if (current_generation > 0):
        FilePathToDB = f"../generation_{previous_generation}/voxels_test/genetic/genetic.db"
        elite_filepaths = GetBestLossFunctions(FilePathToDB,keep_percent)
        #works
    
        #Breed them! 
        bred_filepaths = []
        for i in range(len(elite_filepaths)):
            for j in range(len(elite_filepaths)):
                
                if (j <= i):
                    continue

                bred_filepaths += breed(mutation_rate, elite_filepaths[i],elite_filepaths[j],save_loc_breed,i,j)
                
        #Gets all the filepaths of the new hashmap.txt files
        new_voxel_map_filepaths = elite_filepaths + bred_filepaths

        #Relocates them all to the same single folder. 
        relocate_files(new_voxel_map_filepaths, save_loc_all)

        #Also should put a function which copy pastes the reference and base into the new generational folder:) 

    #Create a new generation if there wasn't a previous one. 
    else:
        initialise_generation(population,nx*ny*nz,fill_percent, save_loc_all)
        pass
    
def current_gen():
    '''Scans one above the current working directory for generation_x folders, the one with the highest value x is taken as the previous geneartion  (returns 0 if there was none -> means to initialise a new gen)'''
    target_prefix = "generation"
    parent_directory = os.path.dirname(os.getcwd())

    matching_folders = find_folders_with_prefix(parent_directory, target_prefix)
    
    prev_gen = extract_highest_number(matching_folders)
    return prev_gen 

def extract_highest_number(strings):
    highest_number = None

    for string in strings:
        match = re.search(r'\d+$', string)
        if match:
            number = int(match.group())
            if highest_number is None or number > highest_number:
                highest_number = number

    return highest_number

def find_folders_with_prefix(parent_folder, prefix):
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
    
    if (os.path.exists(folder)==False):
        os.mkdir(folder)

    for i in range(population):
        array = generate_array_with_percentage(length,fill_percent)
        filepath = f"{folder}/hashmap{i}.txt"
        saveastxt(filepath,array)

  
def generate_array_with_percentage(size, percentage):
    num_ones = int(size * (percentage / 100))
    num_zeros = size - num_ones

    ones_array = np.ones(num_ones, dtype=int)
    zeros_array = np.zeros(num_zeros, dtype=int)

    combined_array = np.concatenate((ones_array, zeros_array))
    np.random.shuffle(combined_array)

    return combined_array


