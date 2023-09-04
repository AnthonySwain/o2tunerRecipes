"""
The cut tuning reference/basline run
"""

from os.path import join, abspath
from os import environ
from platform import system as os_system
import pandas as pd

from o2tuner.system import run_command
from o2tuner.config import resolve_path

O2_ROOT = environ.get("O2_ROOT")
MCSTEPLOGGER_ROOT = environ.get("MCSTEPLOGGER_ROOT")

def reference(inspectors, config):
    """
    Called after arguments have been parsed
    """
    events = config["events"]
    generator = config["generator"]
    engine = config["engine"]
    

    lib_extension = ".dylib" if os_system() == "Darwin" else ".so"
    preload = "DYLD_INSERT_LIBRARIES" if os_system() == "Darwin" else "LD_PRELOAD"

    #To skip or not to skip the ZDC detector is the question at hand! 
    zdc_skip = config['zdc_skip']
    if zdc_skip:
        cmd = f'MCSTEPLOG_NO_MAGFIELD=1 MCSTEPLOG_TTREE=1 {preload}={MCSTEPLOGGER_ROOT}/lib/libMCStepLoggerInterceptSteps{lib_extension} ' \
          f'o2-sim-serial -n {events} -g {generator} -e {engine} --skipModules ZDC'
    
    else: 
        cmd = f'MCSTEPLOG_NO_MAGFIELD=1 MCSTEPLOG_TTREE=1 {preload}={MCSTEPLOGGER_ROOT}/lib/libMCStepLoggerInterceptSteps{lib_extension} ' \
          f'o2-sim-serial -n {events} -g {generator} -e {engine}'
    
    run_command(cmd, log_file=config["o2_sim_log"])
    return True


def baseline(inspectors, config):
    reference_dir = resolve_path(config["reference_dir"])
    kine_file = join(reference_dir, "o2sim_Kine.root")
    steplogger_file = join(reference_dir, "MCStepLoggerOutput.root")
    


    #Is ZDC skip required here?
    zdc_skip = config['zdc_skip']
    if zdc_skip:
        cmd = f'o2-sim-serial -n {config["events"]} -g extkinO2 --extKinFile {kine_file} -e MCReplay --configKeyValues="MCReplayParam.allowStopTrack=true;MCReplayParam.stepFilename={steplogger_file}" --skipModules ZDC'
    
    else:
        cmd = f'o2-sim-serial -n {config["events"]} -g extkinO2 --extKinFile {kine_file} -e MCReplay --configKeyValues="MCReplayParam.allowStopTrack=true;MCReplayParam.stepFilename={steplogger_file}"'
    
    run_command(cmd, log_file=config["o2_sim_log"])

    #extract_hits_root = abspath(join(O2_ROOT, "share", "macro", "analyzeHits.C"))
    extract_hits_root = config['analyzeHitsFilePath']
    cmd_extract_hits = f"root -l -b -q {extract_hits_root}"
    run_command(cmd_extract_hits, log_file="hits.dat")
    return True

#what is "inspectors" - why is this dependency needed, doesn't it do it again in baseline 6 lines above? 
def baseline_hits(inspectors, config):
    #extract_hits_root = abspath(join(O2_ROOT, "share", "macro", "analyzeHits.C"))
    extract_hits_root = config['analyzeHitsFilePath']
    cmd_extract_hits = f"root -l -b -q {extract_hits_root}"
    run_command(cmd_extract_hits, log_file="hits.dat")
    return True

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
    R_list = []

    for i in range(number_of_cylinders):
        Z_loc = Zmin + delta_Z*(0.5 + i)
   
        
        # Filter the rows where Zmin < Z < Zmax
        filtered_rows = data[(data['Zmin'] <= Z_loc) & (data['Zmax'] >= Z_loc)]


        if not filtered_rows.empty:
            filtered_rows_sort = filtered_rows.sort_values(by='R',ascending=False)
            print(filtered_rows_sort)
            R = filtered_rows_sort['R'].iloc[0]
    
        else:
            print("No rows found for the given condition. ", Z_loc)

        Z_min_list.append(Zmin + delta_Z*(i))
        Z_max_list.append(Zmin + delta_Z*(i+1))
        R_list.append(R)

    To_check_list = []

    for i in range(len(R_list)):
        To_check_list.append("All")
    
 
    new_data = {'Zmin' : Z_min_list,
                'Zmax' : Z_max_list,
                'R' : R_list,
                'to_check' : To_check_list}
    
    df = pd.DataFrame(new_data)
    print(df)
    df.to_csv(new_csv_filepath,index=False)

    return(True)

def add_zdc_cylinders(inspectors, config):
    
    csv_file_path = config["csv_filepath_base_data"]
    number_of_cylinders = config["number_of_cylinders"]
    new_csv_filepath = config["csv_filepath_hyper_fine"]
    Z_extent = config["Zextent"]

    data = pd.read_csv(csv_file_path)

    #Adding new column to show the original data does not include new regions to cut out. 
    zdc_old = []
    for index,rows in data.iterrows():
        zdc_old.append("old")

    data['ZDC_add'] = zdc_old

    Zmax = data.max()['Zmax'] 
    Zmin = data.min()['Zmin']

    print(Zmax," <- Zmax | Zmin -> ", Zmin)
    
    half_no_cylinders = round(number_of_cylinders/2)

    delta_Z_negative = (Zmin-Z_extent[0]) / half_no_cylinders
    delta_Z_positive = (Z_extent[1]-Zmax) / half_no_cylinders

    fixed_R_added = config["R_added"]

    
    Z_min_list = []
    Z_max_list = []
    R_list = []

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

    for i in range(len(R_list)):
        To_check_list.append("All")
        ZDC_added.append("new")
    
 
    new_data = {'Zmin' : Z_min_list,
                'Zmax' : Z_max_list,
                'radius' : R_list,
                'to_check' : To_check_list,
                'ZDC_add' : ZDC_added}
    
    df = pd.DataFrame(new_data)

    df2 = df.append(data)
    df2.to_csv(new_csv_filepath,index=False)

    return(True)