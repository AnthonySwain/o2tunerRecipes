from os.path import join, abspath
from os import environ
from platform import system as os_system
import pandas as pd

from o2tuner.system import run_command
from o2tuner.config import resolve_path

O2_ROOT = environ.get("O2_ROOT")
MCSTEPLOGGER_ROOT = environ.get("MCSTEPLOGGER_ROOT")


def add_zdc_cylinders(inspectors, config):
    ''' 
    Takes an input CSV file that was tuned without the ZDC detectors and adds a specified number of 
    cylinders onto each side upto a specified Z value (+-) with set radii
    '''

    csv_file_path = config["csv_filepath_base_data"]
    number_of_cylinders = config["number_of_cylinders"]
    new_csv_filepath = config["csv_filepath_ZDC_added"]
    Z_extent = config["Zextent"]

    data = pd.read_csv(csv_file_path)

    #Adding a new column to show the original data (i.e. the main barrel)
    # to show not to optimise these cylinders
    zdc_old = []
    for index,rows in data.iterrows():
        zdc_old.append("old")

    data['ZDC_add'] = zdc_old


    #Z extent of the final geometry / setting up to create the new cylinders
    Zmax = data.max()['Zmax'] 
    Zmin = data.min()['Zmin']
    
    half_no_cylinders = round(number_of_cylinders/2)

    delta_Z_negative = (Zmin-Z_extent[0]) / half_no_cylinders
    delta_Z_positive = (Z_extent[1]-Zmax) / half_no_cylinders

    #Fixed radii of the cylinders
    fixed_R_added = config["R_added"]

    
    Z_min_list = []
    Z_max_list = []
    R_list = []

    #Creating the cylinders either side of the main barrel
    for i in range(half_no_cylinders):
        Z_loc = Zmax + delta_Z_positive*(0.5 + i)        
        R = fixed_R_added
    

        Z_min_list.append(Zmax + delta_Z_positive*(i))
        Z_max_list.append(Zmax + delta_Z_positive*(i+1))
        R_list.append(R)

    for i in range(half_no_cylinders):
        Z_loc = Zmin - delta_Z_negative*(0.5 + i)        
        R = fixed_R_added
    

        Z_min_list.append(Zmin - delta_Z_negative*(i+1))
        Z_max_list.append(Zmin - delta_Z_negative*(i))
        R_list.append(R)


    To_check_list = []
    ZDC_added = []

    #Columns to show
    # a) That the map will delete all particles
    # b) That these cylinders should be optimised
    for i in range(len(R_list)):
        To_check_list.append("All")
        ZDC_added.append("new")
    
 
    new_data = {'Zmin' : Z_min_list,
                'Zmax' : Z_max_list,
                'radius' : R_list,
                'to_check' : To_check_list,
                'ZDC_add' : ZDC_added}
    
    df = pd.DataFrame(new_data)

    #Combined with the main barrel cylinder mapping
    df2 = df.append(data)

    #Save, ready for optimisation. 
    df2.to_csv(new_csv_filepath,index=False)

    return(True)