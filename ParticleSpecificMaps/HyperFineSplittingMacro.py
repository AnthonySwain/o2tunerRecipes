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


def hyper_fine_splitting(csv_file_path, number_of_cylinders_to_add):
    '''Reads CSV fi'''
    if (number_of_cylinders_to_add <= 0):
        '''error'''
        
        pass
    
    data = pd.read_csv(csv_file_path)
    
    sorted_data = data.sort_values(by='Zmax',ascending=False)

    Zmax = data['Zmax'][0] 
    Zmin = data['Zmax'][-1]

    print(Zmax," <- Zmax | Zmin -> ", Zmin)
    
    delta_Z = (Zmax - Zmin) / (number_of_cylinders_to_add)

    for i in range(number_of_cylinders_to_add):
        Z_loc = Zmin + delta_Z(0.5 + i)
        

    
    pass



def main():
    pass

main() 