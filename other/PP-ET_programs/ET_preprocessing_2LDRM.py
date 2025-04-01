import xarray as xr
import numpy as np
import glob
import os.path
import pickle
from multiprocessing import Pool
import time
from multiprocessing import Process
from multiprocessing import current_process
from itertools import repeat
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

start = time.time()
print("hello")

# %% INPUTS
# file_list_inputs = "/data/keeling/a/erickkc2/f/SAAG3/00_preprocessing/lists_rawfiles-CPM_by4months/2018_months01to04.txt"
# path_out = "/data/keeling/a/erickkc2/f/SAAG3/01_inputs2LDRM/input_files_2LDRM_by4months__4km/"
# time_factor = 24
# file_list_inputs = "/data/keeling/a/erickkc2/a/preprocessing_DL-DRM/SAAG2/lists_rawfiles-nonCPM_byyears/2018edit.txt"
# path_out = "/data/keeling/a/erickkc2/f/SAAG2/input_files_2LDRM-nonCPM_byyears__24km/"
# time_factor = 8 # 24

var_ET = None    # in kg m-2 s-1
var_LH = "LH" #  # in J s-1 m-2
var_temp = "T2B" # in K
# time_step_raw = 60*60  # 60*60*3 nonCPM # 60*60 CPM  # time step in seconds raw data
only_nonnegativeET = True

# coarse_factor = 1

# n_cores = 18
# name_var_2LDRM = "ET24"



def main():
    # %% ############# READING ARGUMENTS FROM COMMAND LINE (overwriting to default in script) ###################
    #https://machinelearningmastery.com/command-line-arguments-for-your-python-script/
    # https://stackoverflow.com/questions/30487767/check-if-argparse-optional-argument-is-set-or-not
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--file_inputs", help="Text file containing list of files to process", required = True)
    parser.add_argument("-o", "--path_outputs", help="Path where outputs will be saved", required=True)
    parser.add_argument("-np", "--nproc", type = int, help="number of processors for parallelization", default=2)
    parser.add_argument("-ns", "--nseconds", help="Number of secs per time step in raw data",required=True)
    parser.add_argument("-tf", "--time_factor", type = int, help="Time factor for ET",required=True)
    parser.add_argument("-cf", "--coarse_factor", type = int, help="Coarsening factor. Default value is 1 (no coarsening).",default=1)
    parser.add_argument("-ovar", "--output_var", help="Variable name to be used for the outputs",required=True)
    parser.add_argument("-rf", "--rem_file", help=("Sometimes the initial and finalfiles are stored at 00:00:00 to 00:00:00, and "
                                                   "the number of files is not multiple of time_factor (is a multiple + 1)."
                                                   "Enter 1 to remove 1st file in computation or 2 to remove the last file"),default=0)
    # args, leftovers = parser.parse_known_args()
    args = vars(parser.parse_args())



    if args["file_inputs"] is not None:
        file_list_inputs = args["file_inputs"]
    if args["path_outputs"] is not None:
        path_out = args["path_outputs"]
    if args["nproc"] is not None:
        n_cores = args["nproc"]
    if args["nseconds"] is not None:
        time_step_raw = eval(args["nseconds"])
        if (not isinstance(time_step_raw, int)):
            # https://stackoverflow.com/questions/22633544/how-to-throw-error-and-exit-with-a-custom-message-in-python
            raise SystemExit("ERROR: nseconds must be an integer or an expression that can be calculated as an integer")
    if args["time_factor"] is not None:
        time_factor = args["time_factor"]
    if args["coarse_factor"] is not None:
        coarse_factor = args["coarse_factor"]
    if args["output_var"] is not None:
        name_var_2LDRM = args["output_var"]
    if args["rem_file"] is not None:
        rem_file = args["rem_file"]


    # %% ########### reading  list of files

    with open(file_list_inputs) as f:
        files = f.read().splitlines()
    print(f"Working on files read from {file_list_inputs }")
    if (len(files) == 0):
        raise SystemExit("The file is empty. No files found!")
    
    # remove a file to have a number of files which is multiple of time_factor
    if (rem_file == 1):
        files = files[1:]
    if (rem_file == 2):
        files = files[:-1]


    # %% ###########  initial variables
    n_timesteps_out = len(files)//time_factor

    if (n_timesteps_out < n_cores):
        n_cores = n_timesteps_out
    print("Number of time steps in ouputs is:", n_timesteps_out)
    print("Number of cores for processing: ", n_cores)

    # %% ###########  calculating dimensions of coarsened grid
    var_name0 = [x for x in [var_ET, var_LH, var_temp] if x is not None][0]
    xda_0 = xr.open_dataset(files[0])[var_name0]
    if coarse_factor >= 1 :
        xda_0 = xda_0.coarsen(south_north = coarse_factor, west_east = coarse_factor, boundary='trim').mean() 
    dims_xda_0 = xda_0.shape
    del xda_0


    # %% ###########  array to store time aggregated data
    basename_wo_ext = os.path.splitext( os.path.basename(file_list_inputs))[0]
    file_out = os.path.join(path_out, f"{basename_wo_ext}_{name_var_2LDRM}.dat")

    arr_out_var =   np.memmap(file_out, dtype=np.float32, mode='w+', shape = (n_timesteps_out,dims_xda_0[1],dims_xda_0[2]))



    # %% ########## PROCESS
    lst = range(n_timesteps_out)
    subsets = np.array_split(lst, n_cores)

    # for subset in subsets:
    #     for el in subset:
    #         files[[i in range(el*(time_factor), (el+1)*time_factor)]] 

    files_subsets = [[files[i] for i in range(subset[0]*time_factor, (subset[-1]+1)*time_factor)] for subset in subsets]
    # print([len(f) for f in files_subsets])

    with Pool(n_cores) as pool: 
        parallel_output =   pool.starmap(process_data2D, zip(files_subsets, repeat([var_ET,var_LH, var_temp]), repeat(time_factor), repeat(time_step_raw), repeat(coarse_factor), repeat(only_nonnegativeET)), chunksize = 1)
        # parallel_output =   pool.map(process_data, subsets, chunksize = 1)

    start_f = time.time() 
    for i, subset in enumerate(subsets):
        arr_out_var[list(subset),...] = parallel_output[i] 
    del parallel_output
    end_f = time.time()
    cum_time = end_f - start_f
    print("cum_time copyiing: ", cum_time, " secs",end="\n")

        
    # f = open(file_out, "wb")
    # f.write(arr_out_var.astype('float32', order = "C"))  # !!! MODIFIED
    # f.close() 
        
    file_pickle = os.path.join(path_out, f"{basename_wo_ext}_{name_var_2LDRM}.pkl")
    f = open(file_pickle, "wb")
    # pickle.dump(arr_all_pr_d.transpose([0,2,1])astype('float32', order = "C"), file = f)
    pickle.dump(arr_out_var.astype('float32', order = "C"), file = f)
    f.close()

    arr_out_var.flush()
    print(file_out + "  written !")

    end = time.time()
    print("Ellapsed time:", end - start, "secs")


##############################################
#### ------------ FUNCTIONS ------------- ####
##############################################

def process_data2D__core(files, vars_related_ET, time_factor, time_step_raw, coarse_factor, only_nonnegative_ET = True, arr_out_var = None):
    default_temp_celcius = 20
    
    current = current_process()
    division = divmod(len(files), time_factor)
    n_timesteps_out = division[0]
    if (division[1] != 0):
        raise ValueError('Length of files should be a multiple of time factor')
    if (arr_out_var is None):   return_arr = True
    else:                       return_arr = False
        
    # for i_timestep_out in range(ntimesteps_out):
    for i_timestep_out in range(n_timesteps_out):
        print(current.name,", i_timestep_out: ", i_timestep_out + 1, " /",n_timesteps_out,"   ",end="\r")
        # with open(os.path.join(path_out, "log.txt"), 'a') as fa:
        #     print(f"Reading file {str(i_timestep_out)} / {ntimesteps_out}   "  ,end="\n", file = fa)
        for i in range(time_factor):
            i_timestep_input = i_timestep_out * time_factor + i
            # with open(os.path.join(path_out, "log.txt"), 'a') as fa:
                # print(f"Reading file {str(i_timestep_input+1)} / {len(files)}   "  ,end="\n", file = fa)
            # print(f"Reading file {str(i_timestep_input+1)} / {len(files)}   "  ,end="\n")
            file = files[i_timestep_input]

            if (vars_related_ET[0] is not None):  # var_ET
                xda = xr.open_dataset(file)[vars_related_ET[0]]
            else:
                xda = xr.open_dataset(file)[vars_related_ET[1]]  # latent heat 
                if (vars_related_ET[2] is not None):
                    xda_temp = xr.open_dataset(file)[vars_related_ET[2]]
                else:
                    xda_temp = default_temp_celcius + 273.15
                #     J. Bringfelt. Test of a forest evapotranspiration model. Meteorology and
                #     Climatology Reports 52, SMHI, Norrköpping, Sweden, 1986.
                # other way to do that may be https://earthscience.stackexchange.com/questions/20733/fluxnet15-how-to-convert-latent-heat-flux-to-actual-evapotranspiration
                Lv = 4185.5 * (751.78 - 0.5655 * (xda_temp))  #     # Calculate lambda or Lv()
                del xda_temp

                xda = (xda / Lv)  # from latent heat to evap ( in kg m-2 s-1 or mm s-1)
             
            if coarse_factor >= 1 :
                xda = xda.coarsen(south_north = coarse_factor, west_east = coarse_factor, boundary='trim').mean()  

            if (i_timestep_input == 0 and  i_timestep_out == 0):
                dims_xda = xda.shape
                arr_temp_var = np.zeros([time_factor,dims_xda[1],dims_xda[2]], order='C')
                if (arr_out_var is None):
                    arr_out_var = np.zeros([n_timesteps_out,dims_xda[1],dims_xda[2]], order='C')

            arr_temp_var[i,:,:] = xda.values * time_step_raw  # multiplied by time_step_raw to transform from  kg m-2 s-1   to  kg m-2
        
        if (only_nonnegative_ET):
            arr_out_var[i_timestep_out,:,:] = np.maximum(np.sum(arr_temp_var, axis = 0),0)
        else:
            arr_out_var[i_timestep_out,:,:] = np.sum(arr_temp_var, axis = 0)

    
    if (return_arr):
        return arr_out_var
    
def process_data2D(files, var_name, time_factor, time_step_raw, coarse_factor, only_nonnegative_ET = True):
    return process_data2D__core(files, var_name, time_factor, time_step_raw,  coarse_factor, only_nonnegative_ET, None)

def process_data2D_inplace(files, var_name, time_factor, time_step_raw, coarse_factor, only_nonnegative_ET, arr_out_var):
    process_data2D__core(files, var_name, time_factor, time_step_raw, coarse_factor, only_nonnegative_ET, arr_out_var)








if __name__ == "__main__":
    main()
