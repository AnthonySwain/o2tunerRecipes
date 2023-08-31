'''
Macro to perform a "hyperfine splitting"
Say how many cylinders are wanted (more than the current optimisation)
Creates a CSV file using the current optimisation but with a larger number of cylinders whose radius can be varied! :) 
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


def hyper_fine_splitting(csv_file_path, number_of_cylinders, new_csv_filepath):
    data = pd.read_csv(csv_file_path)
    
    sorted_data = data.sort_values(by='Zmax',ascending=False)

    Zmax = data.max()['Zmax'] 
    Zmin = data.min()['Zmin']

    print(Zmax," <- Zmax | Zmin -> ", Zmin)
    
    delta_Z = (Zmax - Zmin) / (number_of_cylinders)


    
    Z_min_list = []
    Z_max_list = []
    R_list = []

    for i in range(number_of_cylinders):
        Z_loc = Zmin + delta_Z(0.5 + i)
   
        
        # Filter the rows where Zmin < Z < Zmax
        filtered_rows = data[(data['Zmin'] < Z_loc) & (data['Zmax'] > Z_loc)]


        if not filtered_rows.empty:
            print("Rows corresponding to Zmin < Z < Zmax:")
            print(filtered_rows) #should only be one
            filtered_rows_sort = filtered_rows.sort_values(by='R',ascending=False)
            print(filtered_rows_sort)
            R = filtered_rows_sort['R'][0]
            print("R = ", R)    
        else:
            print("No rows found for the given condition.")

        Z_min_list += Zmin + delta_Z(i)
        Z_max_list += Zmin + delta_Z(i+1)
        R_list += R

    To_check_list = []

    for i in range(len(R_list)):
        To_check_list += "All"
    
    
    new_data = { 'Zmin' : Z_min_list,
            'Zmax' : Z_max_list,
            'radius' : R_list,
            'to_check' : To_check_list}
    
    df = pd.DataFrame(new_data)

    df.to_csv(new_csv_filepath,index=False)




def main():
    csv_filepath = ""
    no_cylinders = 50
    save_new_csv_filepath = ""
    hyper_fine_splitting(csv_filepath,no_cylinders,save_new_csv_filepath)

main() 