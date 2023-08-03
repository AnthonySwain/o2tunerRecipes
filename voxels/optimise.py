"""
The cut tuning optimisation run
"""

import sys
from os.path import join, abspath
from os import environ

from o2tuner.system import run_command
from o2tuner.utils import annotate_trial
from o2tuner.optimise import optimise
from o2tuner.io import parse_json, dump_json, dump_yaml, parse_yaml, exists_file
from o2tuner.optimise import needs_cwd
from o2tuner.config import resolve_path

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
    Retrieve the average number of original and skipped steps
    """
    search_string = "Original number, skipped, kept, skipped fraction and kept fraction of steps"
    extract_start = len(search_string.split()) #whats going on here? 
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
    return sum(steps_orig) / len(steps_orig), sum(steps_skipped) / len(steps_skipped)


def compute_metrics(hits_path, hits_baseline_path, steps_path, steps_baseline_path, o2_detectors):
    """
    Compute the loss and return steps and hits relative to the baseline
    """
    hits_opt = extract_hits(hits_path, o2_detectors, {n: i for i, n in enumerate(o2_detectors)})
    hits_ref = extract_hits(hits_baseline_path, o2_detectors, {n: i for i, n in enumerate(o2_detectors)})
    rel_hits = [h / r if (h is not None and r is not None and r > 0) else None for h, r in zip(hits_opt, hits_ref)] #wtf is happening here 

    steps = extract_avg_steps(steps_path)
    steps_baseline = extract_avg_steps(steps_baseline_path)

    # baseline steps
    steps_baseline = steps_baseline[0]
    steps_remaining = steps[1]
    rel_steps = 1 - (steps_remaining / steps_baseline)

    return rel_steps, rel_hits


def compute_loss(rel_hits, rel_steps, rel_hits_cutoff, penalty_below):
    """
    Compute the loss and return steps and hits relative to the baseline
    """
    rel_hits_valid = [rh for rh in rel_hits if rh is not None] #Gets rid of hits that just say none, need to look into why some say none...
    
    loss = rel_steps**2 #this isn't correct if we want to maximise rel_steps ? -> Misleading name. 1- steps/stepsbase

    for rvh in rel_hits_valid: #I think i understand the logic but I don't think that this is correct way to do it, 
        if rvh < rel_hits_cutoff:
            loss += (penalty_below * (rel_hits_cutoff - rvh))**2
        else:
            loss += (rel_hits_cutoff - rvh)**2

    return loss / (len(rel_hits_valid) + 1)


def run_on_batch(config):
    # in the reference directory we have the MCStepoLoggerOutput.root file
    reference_dir = resolve_path(config['reference_dir']) #what does resolve path do
    kine_file = join(reference_dir, "o2sim_Kine.root")
    steplogger_file = join(reference_dir, "MCStepLoggerOutput.root")
    # in the baseline directory we have the baseline steps and baseline hits
    baseline_dir = resolve_path(config['baseline_dir'])
    sim_log_baseline = join(baseline_dir, config["o2_sim_log"])
    baseline_hits_file = join(baseline_dir, "hits.dat")

    # replay the simulation
    cmd = f'o2-sim-serial -n {config["events"]} -g extkinO2 --extKinFile {kine_file} -e MCReplay ' \
          f'--configKeyValues="MCReplayParam.allowStopTrack=true;MCReplayParam.stepFilename={steplogger_file};GlobalSimProcs.blackholeVoxelFile={config["hashmap_file"]}"'
    _, sim_log = run_command(cmd, log_file="sim.log")

    # extract the hits using O2 macro and pipe to file
    extract_hits_root = abspath(join(O2_ROOT, "share", "macro", "analyzeHits.C"))
    cmd_extract_hits = f"root -l -b -q {extract_hits_root}"
    _, hit_file = run_command(cmd_extract_hits, log_file="hits.dat")

    # compute the loss and further metrics...
    return compute_metrics(hit_file, baseline_hits_file, sim_log, sim_log_baseline, config["O2DETECTORS"])


def sample_voxels(trial, n_voxels, save_file_line_by_line):
    """
    create a simple text file, each line with 0 or 1
    """
    with open(save_file_line_by_line, "w") as f:
        for nv in range(n_voxels):
            on_or_off = trial.suggest_categorical(f"voxel_{nv}", [0, 1]) #trial.suggest_catergorical <- what does this do?
            f.write(str(on_or_off)) #writing as a string -> is there a better way?:) 


def create_hash_map(macro_path, save_file_line_by_line, nx, ny, nz, save_root_hashmap_file):
    """
    extract list from save_file, construct hash map and save to save_root_hashmap_file

    you would need to run a ROOT macro somewhat like
    """
    
    CreateHashMap = f"root -l -b -q {macro_path}({save_file_line_by_line},{save_root_hashmap_file},{nx},{ny},{nz})"
    _, hashmap_file = run_command(CreateHashMap, log_file="hits.dat")
    


def CreateRadialHashMap(trial, RadialMacroPath, Nx, Ny, Nz, RootHashMapSaveLoc,layers):
    """
    figure out the voxels that belong to layer i
    """
    #Function made in dump hashmaps macro - its currently overloaded and needs testing but I don't think overloaded functions are supported
    # in root macros, -> Just call the macro directly, skip making a text file to pass to it. 
    # this will be sampled from the predefined search space. Here, None is simply a placeholder
    
    #Get the predefined search space from the yaml file/
    layer_i = trial.suggest_categorical("i_layer_xy", layers)
    NumbLayers = len(layers)

    XYmax = 1000 
    minRadius = (XYmax/NumbLayers)*layer_i 

    CreateRadialHashMap = f"root -l -b -q {RadialMacroPath}({RootHashMapSaveLoc},{Nx},{Ny},{Nz},{minRadius})"
    _, hashmap_file = run_command(CreateRadialHashMap, log_file="hits.dat")
    

@needs_cwd #what does this decorator do
def objective(trial, config):
    """
    The central objective function for the optimisation
    """

    # construct voxel hash map
    # 1. how many voxels in x, y, z
    #    could use the same logic as here: https://gitlab.cern.ch/bvolkel/VecGeom/-/blob/master/VecGeom/base/FlatVoxelHashMap.h#L163 to only have one 1D list
    #    e.g. write true/false (or 0, 1) to a file which will then be read by ROOT to make the actual HashMap and store it in a ROOT file
    # 2. which to switch on

    # e.g.
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    save_file_line_by_line = config["voxels_sampled_file"]
    save_root_hashmap_file = config["hashmap_file"]
    sample_voxels(trial, nx * ny * nz, save_file_line_by_line)
    create_hash_map(config["CreateHashMapFromTxtMacroFullPath"], save_file_line_by_line, nx, ny, nz, save_root_hashmap_file)

    # rng = np.random.default_rng()
    # batch_id = rng.integers(0, batches)
    rel_steps_avg, rel_hits_avg = run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    annotate_trial(trial, "rel_steps", rel_steps_avg)
    annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"])#, penalty below!)

@needs_cwd #what does this decorator do
def iterate_layers_xy(trial, config):
    """
    This would be a function to simply iterate layer by layer <- I thought I understood what this meant but I do not anymore. Compute loss?
    Is this slowly bringing in the walls around the center of the detector in a cylindrical sense (/ circular in xy place) to determine
    the best outside cutoff for the detectors... 
    """

    # construct voxel hash map
    # 1. how many voxels in x, y, z
    #    could use the same logic as here: https://gitlab.cern.ch/bvolkel/VecGeom/-/blob/master/VecGeom/base/FlatVoxelHashMap.h#L163 to only have one 1D list
    #    e.g. write true/false (or 0, 1) to a file which will then be read by ROOT to make the actual HashMap and store it in a ROOT file
    # 2. which to switch on

    # e.g.
    nx = config["n_voxels_x"]
    ny = config["n_voxels_y"]
    nz = config["n_voxels_z"]
    save_root_hashmap_file = config["hashmap_file"]

    #I hope this works in a similar way but this may be wrong:( 
    layers = config["search_space"]

    #Edit these so it just calls the cpp file directly and makes the circular layer, change so layer is actually just radius (or converetd into a fraction of the radius ig)
    #     make_voxel_layer(trial, nx, ny, nz, save_file_line_by_line)
    create_hash_map(trial, config["CreateRadialHashMapFullPath"], nx, ny, nz, save_root_hashmap_file,layers)

    # rng = np.random.default_rng()
    # batch_id = rng.integers(0, batches)
    rel_steps_avg, rel_hits_avg = run_on_batch(config)

    # annotate drawn space and metrics to trial so we can re-use it
    annotate_trial(trial, "rel_steps", rel_steps_avg)
    annotate_trial(trial, "rel_hits", rel_hits_avg)
    # annotate with other data if you want

    return compute_loss(rel_hits_avg, rel_steps_avg, config["rel_hits_cutoff"])
