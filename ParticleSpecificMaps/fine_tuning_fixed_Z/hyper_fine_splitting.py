"""
Hyper_fine_splitting dependency. (Takes current input of cylinders and cuts them into finer detail.)
"""

from os.path import join, abspath
from os import environ
from platform import system as os_system
import pandas as pd

from o2tuner.system import run_command
from o2tuner.config import resolve_path

O2_ROOT = environ.get("O2_ROOT")
MCSTEPLOGGER_ROOT = environ.get("MCSTEPLOGGER_ROOT")

def hyper_fine_splitting(inspectors, config):
    
    csv_file_path = config["csv_filepath_base_data"]
    number_of_cylinders = config["number_of_cylinders"]
    new_csv_filepath = config["csv_filepath_hyper_fine"]
    

    data = pd.read_csv(csv_file_path)


    Zmax = data.max()['Zmax'] 
    Zmin = data.min()['Zmin']

    print(Zmax," <- Zmax | Zmin -> ", Zmin)
    
    delta_Z = (Zmax - Zmin) / (number_of_cylinders)


    
    Z_min_list = []
    Z_max_list = []
    radius_list = []

    for i in range(number_of_cylinders):
        Z_loc = Zmin + delta_Z*(0.5 + i)
   
        
        # Filter the rows where Zmin < Z < Zmax
        filtered_rows = data[(data['Zmin'] <= Z_loc) & (data['Zmax'] >= Z_loc)]


        if not filtered_rows.empty:
            filtered_rows_sort = filtered_rows.sort_values(by='radius',ascending=False)
            print(filtered_rows_sort)
            radius = filtered_rows_sort['radius'].iloc[0]
    
        else:
            print("No rows found for the given condition. ", Z_loc)

        Z_min_list.append(Zmin + delta_Z*(i))
        Z_max_list.append(Zmin + delta_Z*(i+1))
        radius_list.append(radius)

    To_check_list = []

    for i in range(len(radius_list)):
        To_check_list.append("All")
    
 
    new_data = {'Zmin' : Z_min_list,
                'Zmax' : Z_max_list,
                'radius' : radius_list,
                'to_check' : To_check_list}
    
    df = pd.DataFrame(new_data)
    print(df)
    df.to_csv(new_csv_filepath,index=False)

    return(True)
