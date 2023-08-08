''' Genetic Algorithm to find the best voxelmaps of a generation and breed them. '''

import sqlite3

import sys
from os.path import join, abspath
from os import environ
import random
import os
import numpy as np 

def saveastxt(list,save_path):
    HashMapToSave = ''.join(str(i) for i in list)

    with open(save_path, "w") as file:
        file.write(HashMapToSave)


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


def genetic_optimisation():
    #Insert config into above from yaml file
    #FilePathToDB = config['DBFilepath']
    #keep_percent = config['KeepPercent']
    #mutation_rate = config['mutation_rate']
    #save_loc = config['saveNewMaps']
    
    save_loc = "/children"

    if (os.path.exists(save_loc)==False):
        os.mkdir(save_loc)

    FilePathToDB = "/home/answain/alice/o2tunerRecipes/voxels/iterate_layers_xy.db"
    keep_percent = 0.2
    mutation_rate = 0.2

    #works
    elite_filepaths = GetBestLossFunctions(FilePathToDB,keep_percent)
    
    #Breed them! 
    bred_filepaths = []
    for i in range(len(elite_filepaths)):
        for j in range(len(elite_filepaths)):
            if (j <= i):
                continue

            bred_filepaths.append(breed(mutation_rate, elite_filepaths[i],elite_filepaths[j],save_loc,i,j))
            
    
    new_voxel_map_filepaths = elite_filepaths + bred_filepaths
    
  


genetic_optimisation()

