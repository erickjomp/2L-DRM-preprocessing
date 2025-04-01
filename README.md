# prepro2LDRM-slabvars 


This repository contains the CLI preprocessing program for the 2L-DRM model. This program takes WRF files (wrfout files) as inputs and writes  binary files that will be the inputs of the 2LDRM model. You obtain binary files for: PW, U, V and PWflux. Other variables are written only for validation purposes. 

The generated binary files have dimensions (nx, ny, nslabs, ntimesteps), with the exception of PWflux which has dimensions (nx, ny, nslabs-1, ntimesteps) since it is calcualted at the interface between slabs. The binary files are stored in Fortran-memory order (see Notes below).

Once installed, please use:

    prepro2LDRM-slabvars --help

for more information about the options and arguments of this preprocessing program.


This program is entirely written in Fortran and is computationally very efficient. To build and compile it, you need to have the following installed:
1. A Fortran compiler that supports the module ieee_arithmetic is required (any compiler that implements the Fortran 2003 standard or higher should be fine). I did it with the compiler included in GNU/GCC 9.3.0. 
2. A compatible NetCDF4 libary.
3. A compatible MPI implementation (MPI is used for parallelization capabilities)
4. Fortran Package Manager (fpm). See intructions for installing it here: https://github.com/fortran-lang/fpm. I installed it using a conda environment.

Then you can clone this repository in a local directory: 

    clone https://github.com/erickjomp/2L-DRM-preprocessing.git

And you can build the 2L-DRM program using:

    fpm build

If this does not work, probably you will have to manually specify some linking flags. For instance this worked for me: `fpm build  --link-flag "-L/sw/netcdf4-4.7.4-gnu-9.3.0/lib`

And install it using:

    fpm install


In addition, if using a Linux system, you may want to add `~/.local/bin/` to your PATH environment. For that, you can edit your `~/.profile` file by adding a line like this:

    export PATH=~/.local/bin/:$PATH

Then you will be able to call the program `prepro2LDRM-slabvars` from any directory. For a list of the required and optional inputs of this program, please use:

    prepro2LDRM-slabvars --help


## Notes
- The data in the generated  binary files is in Fortran-like memory order (read more about it here https://manik.cc/2021/02/25/memory-order.html). As a rule of thumb, if you want an input binary file with dimensions (nx, ny, 2, ndays) and you are going to write it from a Python numpy array, the numpy array must have dimensions (ndays, 2, ny, nx).
- The precipitable water flux at the interslab (PWflux) as calculated by $\rho_{air} w_{wind} q$, where $\rho_{air}$ is the density of the air, $w_{wind}$ is the vertical wind velocity and $q$ is the specific humidity at the interslab pressure level.

## Additional binary files required by 2L-DRM
- Precipitation (PP) and Evapotranspiration can be calculated using the python programs in the folder `other_programs/PP-ET_programs`. To know the arguments required by those programs, please use `python ET_preprocessing_2LDRM.py --help` or  `python PP_preprocessing_2LDRM_fromRAINNC-C.py --help` .
- Other scripts for additional inputs can be found in `other_programs/additional_scripts`
