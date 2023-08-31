"""
Contains all the functions used in optimisation runs (see the other py files for the set-ups)
"""

import sys
from os.path import join, abspath
from os import environ
import os
import shutil
from o2tuner.system import run_command
from o2tuner.utils import annotate_trial
from o2tuner.optimise import optimise
from o2tuner.io import parse_json, dump_json, dump_yaml, parse_yaml, exists_file
from o2tuner.optimise import needs_cwd
from o2tuner.config import resolve_path
import numpy as np
import subprocess

# Get environment variables we need to execute some cmds
O2_ROOT = environ.get("O2_ROOT") #Just the path to the O2 og enviroment
MCSTEPLOGGER_ROOT = environ.get("MCSTEPLOGGER_ROOT")


def extract_hits(path, o2_detectors, det_name_to_id):
    """
    Retrieve the number of hits per detector
    """
    hits = [None] * len(o2_detectors)
    with open(path, "r", encoding="utf8") as hit_file:
        for line in hit_file:
            fields = line.split()
            if not fields:
                continue
            if fields[0] in o2_detectors:
                pot_nan = fields[1].lower()
                if "nan" in pot_nan:
                    hits[det_name_to_id[fields[0]]] = None
                    continue
                # NOTE relying on the first field containing the number of hits, this may change if the O2 macro changes
                hits[det_name_to_id[fields[0]]] = float(fields[1])
        return hits

def extract_avg_steps(path):
    """
    Retrieve the average number of original and skipped steps over all events 
    """
    search_string = "Original number, skipped, kept, skipped fraction and kept fraction of steps"
    extract_start = len(search_string.split()) 
    steps_orig = []
    steps_skipped = []
    with open(path, "r", encoding="utf8") as step_file:
        for line in step_file:
            if search_string in line:
                line = line.split()
                steps_orig.append(int(line[extract_start]))
                steps_skipped.append(int(line[extract_start + 1]))
    if not steps_orig:
        print("ERROR: Could not extract steps")
        sys.exit(1)


    #I am not sure why this is divided by length of the list but sure
    return sum(steps_orig) / len(steps_orig), sum(steps_skipped) / len(steps_skipped)

def compute_metrics(hits_path, hits_baseline_path, steps_path, steps_baseline_path, o2_detectors,loss_data_save):
    """
    Compute the loss and return steps and hits relative to the baseline
    """
    hits_opt = extract_hits(hits_path, o2_detectors, {n: i for i, n in enumerate(o2_detectors)})
    hits_ref = extract_hits(hits_baseline_path, o2_detectors, {n: i for i, n in enumerate(o2_detectors)})
    rel_hits = [h / r if (h is not None and r is not None and r > 0) else None for h, r in zip(hits_opt, hits_ref)] #wtf is happening here 

    steps = extract_avg_steps(steps_path)
    steps_baseline = extract_avg_steps(steps_baseline_path)

    # baseline steps
    # '''
    # INTRICATE BEHAVIOUR THATS UNEXPECTED
    # steps_baseline[0] will always be the same as steps[0] as this is the same number as the reference simulation. BUT.
    # steps_baseline also has skipped steps - the optimisation will have the same replay as this, not the reference, therefore we need to comapre the number of steps in the baseline
    # versues the simulation (i.e take into account skipped steps in the baseline to give the total number of steps in the baseline)
    # '''
    steps_base = steps_baseline[0] - steps_baseline[1]
    steps_optimise = steps[0]-steps[1]
    rel_steps = steps_optimise/steps_base
    
    # Write data to the file so it can be easily found post-simulations. 
    #loss_calc_data_save_filepath = ""
    with open (loss_data_save, "a") as file:
        file.write(f"Hits/HitsBase || ")
        for i, n in enumerate(o2_detectors):
            file.write(f"{n} : {rel_hits[i]} | ")

        file.write(f"\n Steps/Steps_Base : {rel_steps}")

    return rel_steps, rel_hits

def compute_loss(rel_hits, rel_steps, rel_hits_cutoff, penalty_below):
    """
    Compute the loss and return steps and hits relative to the baseline
    """
    rel_hits_valid = [rh for rh in rel_hits if rh is not None] #Gets rid of hits that just say none (i.e detectors whos analyzeHits does not work...)
    
    loss = rel_steps**2 

    for rvh in rel_hits_valid:
        if rvh < rel_hits_cutoff:
            loss += (penalty_below * (1 - rvh))**2
        else:
            loss += (1 - rvh)**2

    #The normalisation over rel_hits_vald doesn't really achieve anything here does it now... 
    return loss 

def run_on_batch(config):
    # in the reference directory we have the MCStepoLoggerOutput.root file
    reference_dir = resolve_path(config['reference_dir']) 
    kine_file = join(reference_dir, "o2sim_Kine.root")
    steplogger_file = join(reference_dir, "MCStepLoggerOutput.root")
    # in the baseline directory we have the baseline steps and baseline hits
    baseline_dir = resolve_path(config['baseline_dir'])
    sim_log_baseline = join(baseline_dir, config["o2_sim_log"])
    baseline_hits_file = join(baseline_dir, "hits.dat")

    loss_data_save_file = config['Loss_data_save_file']

    # replay the simulation
    #Is ZDC skip required here? 
    zdc_skip = config['zdc_skip']
    NoVoxelMap = config['NoVoxelMap']
    if zdc_skip:
        if NoVoxelMap:
            cmd = f'o2-sim-serial -n {config["events"]} -g extkinO2 --extKinFile {kine_file} -e MCReplay ' \
                f'--configKeyValues="MCReplayParam.allowStopTrack=true;MCReplayParam.stepFilename={steplogger_file};GlobalSimProcs.GeoCutsCSVFile={config["csv_filepath_write"]}" --skipModules ZDC'

        else:
             cmd = f'o2-sim-serial -n {config["events"]} -g extkinO2 --extKinFile {kine_file} -e MCReplay ' \
                f'--configKeyValues="MCReplayParam.allowStopTrack=true;MCReplayParam.stepFilename={steplogger_file};GlobalSimProcs.blackholeVoxelFile={config["txt_of_maps"]}" --skipModules ZDC'

    
    
    
    else:
        if NoVoxelMap:
            cmd = f'o2-sim-serial -n {config["events"]} -g extkinO2 --extKinFile {kine_file} -e MCReplay ' \
                f'--configKeyValues="MCReplayParam.allowStopTrack=true;MCReplayParam.stepFilename={steplogger_file};GlobalSimProcs.GeoCutsCSVFile={config["csv_filepath_write"]}"'
        
        else:
            cmd = f'o2-sim-serial -n {config["events"]} -g extkinO2 --extKinFile {kine_file} -e MCReplay ' \
                f'--configKeyValues="MCReplayParam.allowStopTrack=true;MCReplayParam.stepFilename={steplogger_file};GlobalSimProcs.blackholeVoxelFile={config["txt_of_maps"]}"'
        
    

    _, sim_log = run_command(cmd, log_file="sim.log")

    # extract the hits using O2 macro and pipe to file
    #extract_hits_root = abspath(join(O2_ROOT, "share", "macro", "analyzeHits.C"))
    extract_hits_root = config['analyzeHitsFilePath'] #Uses custom analyzeHits file becasue of faulty detector analysis (hence not using a path above)
    cmd_extract_hits = f"root -l -b -q {extract_hits_root}"
    _, hit_file = run_command(cmd_extract_hits, log_file="hits.dat")

    # compute the loss and further metrics...
    return compute_metrics(hit_file, baseline_hits_file, sim_log, sim_log_baseline, config["O2DETECTORS"],loss_data_save_file)

def sample_voxels(trial, n_voxels, save_file_line_by_line,map_creation_macro_fullpath):
    """
    create a simple single line .txt file with n_voxels of 0's or 1's
    """

    #Okay - perhaps every 1000 or so parameters, they should be saved to the txt file, reset the binary list and continue
    #Also not using a goddamn .txt file 

    binary_list = np.array([])
    for nv in range(n_voxels):
        np.append(binary_list,trial.suggest_categorical(f"voxel_{nv}",[0,1]))
        
    # Chunk size for sending data
    chunk_size = 1000
    #np.savetxt(save_file_line_by_line, binary_list, fmt='%d', delimiter='', newline='')

    # Use a loop to send data in chunks through a pipe
    for i in range(0, len(binary_list), chunk_size):
        chunk = binary_list[i:i + chunk_size]

        # Pack the chunk as bytes
        data_bytes = chunk.tobytes()

        # Call the C++ program and pass the data bytes as stdin
        process = subprocess.Popen([map_creation_macro_fullpath], stdin=subprocess.PIPE)
        process.stdin.write(data_bytes)
        process.stdin.close()
        process.wait()
            
    '''
    with open(save_file_line_by_line, "w") as f:
        for nv in range(n_voxels):
            on_or_off = trial.suggest_categorical(f"voxel_{nv}", [0, 1]) 
            f.write(str(on_or_off))
    '''

def create_hash_map(macro_path, rel_txtfilepath, nx, ny, nz, rel_root_hashmap_saveloc):
    """
    Creates hashmap from saved .txt file from sample_voxels
    """
    
    CreateHashMap = f"root -l -b -q '{macro_path}(\"{rel_txtfilepath}\",\"{rel_root_hashmap_saveloc}\",{nx},{ny},{nz})'"
    _, hashmap_file = run_command(CreateHashMap, log_file="hits.dat")
    
def CreateRadialHashMap(trial, RadialMacroPath, Nx, Ny, Nz, RootHashMapSaveLoc,innerRadius,Zextent = None):
    """
    Creates radial hashmap in the XY plane (i.e a cylinder where the cylindrical axis is the Z axis)
    RadialMacroPath = path to the root macro which creates the hashmap
    Nx,Ny,Nz = number of bins in the X,Y,Z directions respectively
    RootHashMapSaveLoc = where to save the resulting hashmap
    Inner radius - inner radius of the cylinder
    Zextent = [ZMin,Zmax] (optional)
    """

    if Zextent == None:
        #Macro command
        CreateRadialHashMap = f"root -l -b -q '{RadialMacroPath}(\"{RootHashMapSaveLoc}\",{Nx},{Ny},{Nz},{innerRadius})'"

    #I.e using the cylinder_xy algorithm (fitting Z extent too)
    else: 
        minZchosen = Zextent[0]
        maxZchosen = Zextent[1]
        CreateRadialHashMap = f"root -l -b -q '{RadialMacroPath}(\"{RootHashMapSaveLoc}\",{Nx},{Ny},{Nz},{innerRadius},{minZchosen},{maxZchosen})'"

    #Runs the macro
    _, hashmap_file = run_command(CreateRadialHashMap, log_file="hits.dat")
    
def relocate_file(filepath, destination):
    '''
    Relocates (by copy & paste) .txt files to a new location
    filepath = original location
    destination = final location
    '''
    if filepath.endswith(".txt"):
        try:
            filename = os.path.basename(filepath)
            shutil.copy(filepath, destination)
        except Exception as e:
            print(f"Error copying {filepath}: {e}")

def add_map_to_txt(txt_filepath,map_filepath, particles):
    '''
    Adds a map to the txt file which contains all the information about the maps.
    '''
    particle_string = "[" + ",".join(str(p) for p in particles) + "]"
    line = f"{particle_string}|{map_filepath}\n"

    try:
        # Try to open the file in append mode
        with open(txt_filepath, 'a') as file:
            file.write(line)
    except FileNotFoundError:
        # If the file doesn't exist, create it and write the line
        with open(txt_filepath, 'w') as file:
            file.write(line)




























@needs_cwd #Ran in its on directory (cwd = current working directory)
def objective(trial, config):
    """
    The central objective function for the optimisation
    """

    # construct voxel hash map
    # 1. how many voxels in x, y, z
    #    could use the same logic as here: https://gitlab.cern.ch/bvolkel/VecGeom/-/blob/master/VecGeom/base/FlatVoxelHashMap.h#L163 to only have one 1D list
    #    e.g. write true/false (or 0, 1) to a file which will then be read by ROOT to make the actual HashMap and store it in a ROOT file
    # 2. which to switch on

    #Gets whats needed from teh config file
    penalty_below = config["penalty_below"]
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    save_file_line_by_line = config["voxels_sampled_file"]
    save_root_hashmap_file = config["hashmap_file"]

    #Creates the hashmap (first by writing 0,1s to a .txt file (sample_voxels) then reading this .txt file (create_hash_map))
    sample_voxels(trial, nx * ny * nz, save_file_line_by_line)
    create_hash_map(config["CreateHashMapFromTxtMacroFullPath"], save_file_line_by_line, nx, ny, nz, save_root_hashmap_file)

    #Run it 
    rel_steps_avg, rel_hits_avg = run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    annotate_trial(trial, "rel_steps", rel_steps_avg)
    annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"], penalty_below)
















@needs_cwd #Ran in its on directory (cwd = current working directory)
def iterate_layers_xy(trial, config):
    """
    Works from outside in with a radial hashmap (cylindrical with the cylindrical axis corresponding to the beam axis)
    Everything outside the cylinder is a blackhole. 
    """

    #Get neccessary information from the config file
    penalty_below = config["penalty_below"]
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    save_root_hashmap_file = config["hashmap_file"]
    Rmax = config["Rmax"]

    #Gets the list of layers from the .yaml file
    layers = config["search_space"]["i_layer_xy"]

    #Creates the radial hashmap. 
    CreateRadialHashMap(trial, config["CreateRadialHashMapFullPath"], nx, ny, nz, save_root_hashmap_file,Rmax,layer_selection=layers)

    # rng = np.random.default_rng()
    # batch_id = rng.integers(0, batches)

    #Run
    rel_steps_avg, rel_hits_avg = run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    annotate_trial(trial, "rel_steps", rel_steps_avg)
    annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"],penalty_below)




















@needs_cwd #Ran in its on directory (cwd = current working directory)
def cylinder_xy(trial, config):
    """
    Works from outside in with a radial hashmap (cylindrical with the cylindrical axis corresponding to the beam axis)
    Everything outside the cylinder is a blackhole. 
    """

    #Get neccessary information from the config file
    penalty_below = config["penalty_below"]
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    save_root_hashmap_file = config["hashmap_file"]

    #Min/Max parameters for the cylinder
    Rmax = config["Rmax"]
    Zmax = config["Zmax"]
    Zmin = config["Zmin"]

    #Creates the radial hashmap. 
    CreateRadialHashMap(trial, config["CreateRadialHashMapFullPath"], nx, ny, nz, save_root_hashmap_file,Rmax,Zextent=[Zmin,Zmax])

    # rng = np.random.default_rng()
    # batch_id = rng.integers(0, batches)

    #Run
    rel_steps_avg, rel_hits_avg = run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    annotate_trial(trial, "rel_steps", rel_steps_avg)
    annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"],penalty_below)
