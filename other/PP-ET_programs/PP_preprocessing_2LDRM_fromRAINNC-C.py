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

# See references:
# https://forum.mmm.ucar.edu/threads/how-to-set-the-prec_acc_dt-properly.5357/
# https://mailman.ucar.edu/pipermail/ncl-talk/2015-May/002521.html

start = time.time()
print("hello")

# %% INPUTS
# file_list_inputs = "/data/keeling/a/erickkc2/a/preprocessing_DL-DRM/SAAG2/lists_rawfiles-nonCPM_byyears__forprecipRAINNC-C/2019.txt"
# path_out = "/data/keeling/a/erickkc2/f/SAAG2/input_files_2LDRM-nonCPM_byyears__24km/"
# time_factor = 8 # 24

# yyyy = 2018

# var_names = {"PP": "PREC_ACC_NC", "PPC": "PREC_ACC_C", "ET": "ACETLSM"}    # for non-convection permitting models you would be to add PREC_ACC_NC and PREC_ACC_C
vars_raincum = ["RAINNC", "RAINC"] # ["PREC_ACC_NC", "PREC_ACC_C"]   for non-convectuve permitting
# bucket_size = 100        # use None in case buckets option was not used in WRF
vars_buckets = ["I_RAINNC","I_RAINC"] # USE None in case there is 


coarse_factor = 1

n_cores = 18

name_var_2LDRM = "PP24"
# names_vars_out = {"PP24", "ET24"}


def main():
    # %% ############# READING ARGUMENTS FROM COMMAND LINE (overwriting to default in script) ###################
    #https://machinelearningmastery.com/command-line-arguments-for-your-python-script/
    # https://stackoverflow.com/questions/30487767/check-if-argparse-optional-argument-is-set-or-not
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--file_inputs", help="Text file containing list of files to process", required = True)
    parser.add_argument("-o", "--path_outputs", help="Path where outputs will be saved", required=True)
    parser.add_argument("-np", "--nproc", type = int, help="number of processors for parallelization", default=2)
    parser.add_argument("-tf", "--time_factor", type = int, help="Time factor for ET",required=True)
    parser.add_argument("-cf", "--coarse_factor", type = int, help="Coarsening factor. Default value is 1 (no coarsening).",default=1)
    parser.add_argument("-ovar", "--output_var", help="Variable name to be used for the outputs",required=True)
    parser.add_argument("-cpm", "--isCPM", type = int,  
                        help= ("If the model is a Convection Permitting Model (CPM), the variable PREC_ACC_C will be ignored"
                               " and only PREC_ACC_NC will be used. If it is a non-CPM model both variables will be used."
                               " Use 1 if it is a CPM model or 0 if it is not."),required=True)
    parser.add_argument("-b", "--is_bucket", type = int, help=("RAINNC and RAINC were saved using the bucket option in WRF?."
                                                    " Use 1 if it is true or 0 for false, Default = 1"),default=1)
    parser.add_argument("-bs", "--size_bucket", type = int, help=("if is_bucket is true, use this variable to set the size of the bucket"
                                                    " Use 1 if it is true or 0 for false, Default = 100"),default=100)
    # args, leftovers = parser.parse_known_args()
    args = vars(parser.parse_args())



    if args["file_inputs"] is not None:
        file_list_inputs = args["file_inputs"]
    if args["path_outputs"] is not None:
        path_out = args["path_outputs"]
    if args["nproc"] is not None:
        n_cores = args["nproc"]
    if args["time_factor"] is not None:
        time_factor = args["time_factor"]
    if args["coarse_factor"] is not None:
        coarse_factor = args["coarse_factor"]
    if args["output_var"] is not None:
        name_var_2LDRM = args["output_var"]
    if args["isCPM"] is not None:
        isCPM = args["isCPM"]
        if (isCPM == 1):
            vars_raincum = ["RAINNC"] 
            vars_buckets = ["I_RAINNC"] 
        elif (isCPM == 0):
            vars_raincum = ["RAINNC", "RAINC"] 
            vars_buckets = ["I_RAINNC","I_RAINC"] 
        else:
            raise SystemExit("ERROR:  isCPM should be 0 or 1")
    if args["is_bucket"] is not None:
        is_bucket = args["is_bucket"]
    
    if args["size_bucket"] is not None:
        bucket_size = args["size_bucket"]
    if (not is_bucket) :
        bucket_size = None   # set to None cause in  process_data2D, bucket_size = None  will be used to identify if RAIN... was saved withoutusing the bucket option of WRF

    # %% ############## getting list of files #########
    # files = glob.glob(path_inputs + pattern_files, recursive=True)
    # files = sorted(files)
    with open(file_list_inputs) as f:
        files = f.read().splitlines()
    print(f"Number of files: {len(files)}")
    print(f"Working on files read from {file_list_inputs }")

    # %% ############ initial variables 3##############
    n_timesteps_out = len(files)//time_factor

    # %% ########### calculating dimensions of coarsened grid #############
    var_name0 = vars_raincum[0]
    xda_0 = xr.open_dataset(files[0])[var_name0]
    if coarse_factor >= 1 :
        xda_0 = xda_0.coarsen(south_north = coarse_factor, west_east = coarse_factor, boundary='trim').mean() 
    dims_xda_0 = xda_0.shape
    del xda_0

    # %% ########### array to store time aggregated data #################
    basename_wo_ext = os.path.splitext( os.path.basename(file_list_inputs))[0]
    file_out = os.path.join(path_out, f"{basename_wo_ext}_{name_var_2LDRM}.dat")

    arr_out_var =   np.memmap(file_out, dtype=np.float32, mode='w+', shape = (n_timesteps_out,dims_xda_0[1],dims_xda_0[2]))







    # %% ######################## PROCESS ###############################
    lst = range(n_timesteps_out)
    subsets = np.array_split(lst, n_cores)

    # for i_subset in range(1,len(subsets)):
    #     subsets[i_subset] = np.concatenate([subsets[i_subset-1][-1:], subsets[i_subset]])

    files_subsets = [[files[1:][i] for i in range(subset[0]*time_factor, (subset[-1]+1)*time_factor)] for subset in subsets]

    files_subsets[0] = files[:1] + files_subsets[0]
    for i_subset in range(1,len(subsets)):
        files_subsets[i_subset] = files_subsets[i_subset-1][-1:] + files_subsets[i_subset] #np.concatenate([subsets[i_subset-1][-1:], subsets[i_subset]])

    with Pool(n_cores) as pool: 
        parallel_output =   pool.starmap(process_data2D, zip(files_subsets, repeat(vars_raincum), repeat(time_factor), repeat(coarse_factor), repeat(bucket_size)), chunksize = 1)
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
        
    print(file_out + "  written !")

    end = time.time()
    print("Ellapsed time:", end - start, "secs")



def process_data2D__core(files, vars_names, time_factor, coarse_factor, bucket_size, arr_out_var = None):
    current = current_process()
    division = divmod(len(files)-1, time_factor)
    n_timesteps_out = division[0]
    if (division[1] != 0):
        raise ValueError('Length of files - 1  should be a multiple of time factor')
    if (arr_out_var is None):   return_arr = True
    else:                       return_arr = False
        
    # for i_timestep_out in range(ntimesteps_out):
    # probably improve to not have to be using all files, but only those separated by the target timestep ouput
    list_prcum = []

    def read_vars(this_file):
        xds = xr.open_dataset(this_file)
        list_xda = [None] * len(vars_names)
        for i_var, var_name in enumerate(vars_names):
            var_bucket = vars_buckets[i_var]
            list_xda[i_var] = xds[var_name]
            if (bucket_size != None):  # correction in case of bucket option
                list_xda[i_var] = list_xda[i_var] + bucket_size * xds[var_bucket]
            if coarse_factor >= 1 :
                list_xda[i_var] = list_xda[i_var].coarsen(south_north = coarse_factor, west_east = coarse_factor, boundary='trim').mean()
        return(list_xda)

    list_xda = read_vars(files[0])
    list_prcum.append(sum([xda.values for xda in list_xda]))

    if (arr_out_var is None):
        dims_xda = list_prcum[0].shape
        arr_out_var = np.zeros([n_timesteps_out,dims_xda[1],dims_xda[2]], order='C')

    
    for i_timestep_out in range(n_timesteps_out):
        print(current.name,", i_timestep_out: ", i_timestep_out + 1, " /",n_timesteps_out,"   ",end="\r")
        
        i_timestep_input_end = i_timestep_out * time_factor + time_factor

        file = files[i_timestep_input_end]

        list_xda = read_vars(file)
        list_prcum.append(sum([xda.values for xda in list_xda]))

        arr_out_var[i_timestep_out,:,:] = list_prcum[1] - list_prcum[0]
        list_prcum.pop(0)
    
    if (return_arr):
        return arr_out_var
        
    
def process_data2D(files, var_name, time_factor, coarse_factor, bucket_size):
    return process_data2D__core(files, var_name, time_factor, coarse_factor, bucket_size, None)

def process_data2D_inplace(files, var_name, time_factor, coarse_factor, bucket_size, arr_out_var):
    process_data2D__core(files, var_name, time_factor, coarse_factor, bucket_size, arr_out_var)


if __name__ == "__main__":
    main()
