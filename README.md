# prepro2LDRM-slabvars


This repository contains the CLI preprocessing program for the 2L-DRM model. This program takes WRF files (wrfout files) as inputs and writes the files that will be the inputs of the 2L-DRM model (Dominguez et al., 2020; Claros et al., 2025): for each slab of the atmosphere,

| Variable | Description | Units | Symbol in Claros et al. (2025) |
|---|---|---|---|
| PW | precipitable water of the slab: water vapor contained in a unit area column of the slab | kg m-2 | w |
| U, V | moisture-weighted average of the wind velocity in the slab, along the grid x and y axes | m s-1 | u, v |
| PWflux | vertical moisture flux at the interslab level, accumulated over the output step: upward if positive, downward if negative | kg m-2 | Fu, Fd |

and at the surface:

| Variable | Description | Units | Symbol in Claros et al. (2025) |
|---|---|---|---|
| PP | precipitation accumulated over the output step | kg m-2 | P |
| ET | evapotranspiration accumulated over the output step | kg m-2 | E |

Other variables at the interslab level (W: vertical wind, Q: specific humidity, moistdensity: density of moist air) are written only for validation purposes.

The program runs on Linux and Windows, in parallel with MPI or as a serial program.

Once installed, please use:

    prepro2LDRM-slabvars --help

for more information about the options and arguments of this preprocessing program.


## Outputs

Each output time step averages `--timefactor` consecutive WRF files (e.g. `-tf 3` with hourly files gives 3-hourly outputs); PWflux, PP and ET are accumulated over the step. The output files are named after the list of files given with `--listfiles`, without its extension (`AMZ2_CPM01.txt` gives `AMZ2_CPM01`):

- `--format nc` (default): one NetCDF file, `AMZ2_CPM01.nc`, with all the variables, the time of each output step (start of the step, in hours since the first WRF time), XLAT/XLONG, the pressure of the interslab levels and the options of the run.
- `--format bin`: one raw binary file per variable, `AMZ2_CPM01_PW.dat`, `AMZ2_CPM01_U.dat`, ..., as read by 2L-DRM. They have dimensions (nx, ny, nslabs, ntimesteps), with the exception of PWflux, W, Q and moistdensity which have dimensions (nx, ny, nslabs-1, ntimesteps) since they are calculated at the interface between slabs, and PP and ET, (nx, ny, ntimesteps). The binary files are stored in Fortran-memory order (see Notes below).

With `--split month` (or `year`) a separate set of outputs is written for each month (or year) of the WRF times, e.g. `AMZ2_CPM01_2018-05.nc` or `AMZ2_CPM01_2018-05_PW.dat`. Every output step must be within one month: the list of files should start at 00:00 of a day (see `--remfile`).

The daily PW24 and PWflux24 outputs are no longer written; `--extratimefactor` (`-etf`) is still accepted but ignored.

To keep using existing scripts that feed 2L-DRM with binary files, add `--format bin` to them.

### PP and ET

PP and ET are written by default (they replace the former python programs of `other/PP-ET_programs`, with the same results). `--no_pp_et` skips them, e.g. for WRF files without these variables.
- PP is the accumulated rain at the end of the output step minus at its start: `RAINNC + RAINC + BUCKET_MM x (I_RAINNC + I_RAINC)`, with the bucket size of the `BUCKET_MM` attribute of the WRF files (`RAINC`, the bucket counters or `BUCKET_MM` count as 0 if missing, as in convection-permitting runs or without the bucket option). So it needs the WRF file at the end of the last step (or at the start of the first one): the list of files must have [n x timefactor + 1] files, and `--remfile` must be 2 (the default) or 1. Output step n goes from file (n-1) x timefactor to file n x timefactor of the whole list.
- ET is the sum over the files of the output step of `LH / Lv` times the seconds between files, where `Lv = 4185.5 (751.78 - 0.5655 T2B)` (J kg-1) is the latent heat of vaporization; negative values are set to 0.

The seconds between WRF files (for PWflux and ET) are calculated from the `Times` of the first two files; `--nseconds` is optional and must agree with them.


## Build

This program is entirely written in Fortran and is computationally very efficient. To build and compile it, you need to have the following installed:
1. A Fortran compiler that supports Fortran 2018 `stop ..., quiet=` (needed by FLAP), e.g. gfortran 12 or newer, or ifx.
2. NetCDF-Fortran (with NetCDF-C), built with a compatible compiler: the program uses its `netcdf` module.
3. Optionally, an MPI library (MPI is used for parallelization capabilities): Open MPI or MPICH on Linux, MS-MPI on Windows.
4. Fortran Package Manager (fpm), version 0.13 or newer. See intructions for installing it here: https://github.com/fortran-lang/fpm. It can be installed in a conda environment (`conda install -c conda-forge fpm`).

The other Fortran packages (FLAP and the Fortran standard library, at the same versions as 2L-DRM) are downloaded by fpm.

Then you can clone this repository in a local directory:

    git clone https://github.com/erickjomp/2L-DRM-preprocessing.git

And you can build the program with the `release` profile (optimized), giving the directory of the `netcdf` module, which `nf-config` prints:

    fpm build --profile release --flag "$(nf-config --fflags)"

Without `--profile release`, fpm uses its `debug` profile (no optimization, with run-time checks), which is much slower. The `--flag` is needed when the `netcdf` module is not in the compiler's default search path, as on RHEL-based systems (`/usr/lib64/gfortran/modules`) or conda; it does no harm elsewhere.

If the NetCDF libraries are not found when linking, add their directory, e.g. `--link-flag "-L/path/to/netcdf/lib"`. A specific compiler is chosen with `--compiler`, e.g. `--compiler gfortran-13`.

And install it using (with the same options as in `fpm build`):

    fpm install --profile release --flag "$(nf-config --fflags)"

In addition, if using a Linux system, you may want to add `~/.local/bin/` to your PATH environment. For that, you can edit your `~/.profile` file by adding a line like this:

    export PATH=~/.local/bin/:$PATH

Then you will be able to call the program `prepro2LDRM-slabvars` from any directory.

### Serial build (no MPI)

To build without MPI, remove `USE_MPI` from `macros` in the `[preprocess.cpp]` section of `fpm.toml`, and the line `mpi = "*"` from `[dependencies]`. The program then runs without `mpirun` and processes all the files in one process.

### Windows

The simplest way is [MSYS2](https://www.msys2.org) (UCRT64 environment):

    pacman -S git mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-netcdf-fortran mingw-w64-ucrt-x86_64-fpm
    fpm build --profile release --flag "-I/ucrt64/include"

For the parallel version, also install `mingw-w64-ucrt-x86_64-msmpi` and the MS-MPI runtime (`msmpisetup.exe` from Microsoft); otherwise use the serial build. On Windows:
- WRF file names contain `:`, which Windows does not allow: rename them (e.g. `wrfout_d01_2018-05-31_00_00_00`) and use the new names in the list of files.
- Paths longer than 260 characters need long paths to be enabled in Windows.

The Windows builds are run in the GitHub Actions workflow, but have not been checked on a real machine yet.


## Tests

`test/regression` has a regression test with pytest: it runs the program on a small case cropped from the AMZ2 CPM 4 km run (16 x 14 cells over the Andes, 2 days) and compares the outputs with the original run. It needs Python with numpy, netCDF4 and pytest:

    python -m pip install numpy netCDF4 pytest
    python -m pytest test/regression --exe $(fpm run prepro2LDRM-slabvars --profile release --flag "$(nf-config --fflags)" --runner echo | tail -1)

Add `--serial` for a serial build, and `--mpiexec PATH` if `mpiexec` is not in PATH. After an intended change in the results, rerun with `--update-reference` and commit the new references. The tests also run on GitHub Actions (`.github/workflows/tests.yml`) for each push.

### keeling

The default gfortran (11) cannot compile FLAP; use the gcc-toolset-13 one with the system Open MPI and NetCDF-Fortran:

    module load mpi/openmpi-x86_64
    fpm build --profile release --compiler /opt/rh/gcc-toolset-13/root/usr/bin/gfortran --flag "$(nf-config --fflags)"

On the login node, MPI runs need `export OMPI_MCA_pml=ob1` (the default UCX transport fails there with "Connection refused"); the tests set it. Jobs on the compute nodes do not need it. The tests also let Open MPI start more processes than there are cores.


## Notes
- The data in the generated  binary files is in Fortran-like memory order (read more about it here https://manik.cc/2021/02/25/memory-order.html). As a rule of thumb, when reading the binary files  as arrays in python, you should reshape it to dimensions (ntimes, nslabs, ny, nx), except PWflux which has to be reshaped to (ntimes, nslabs-1, ny, nx). The NetCDF files have the same dimensions in this order: (time, slab, south_north, west_east). Each NetCDF variable has a short `long_name`, a `description` with its definition (as in the `description` attribute of WRF files) and a `comment` about NaN values.

- The precipitable water flux between slabs (PWflux) is calculated by $\rho_{air} w_{wind} q$, where $\rho_{air}$ is the density of the air, $w_{wind}$ is the vertical wind velocity and $q$ is the specific humidity at the interslab level. It is accumulated over each output step (kg m-2).

- Where the ground is above an interslab level, the variables of the slab below it and at that level are NaN.

- U and V are along the axes of the WRF grid (as the WRF U and V), not true east and north.

## Additional binary files required by 2L-DRM
- Precipitation (PP) and evapotranspiration (ET) are now written by this program (see Outputs). The former python programs are kept in `other/PP-ET_programs` for reference.
- Other scripts for additional 2L-DRM inputs can be found in `other/additional_scripts`

## References
- Claros, E., Dominguez, F., & Liu, C. (2025). Moisture sources of precipitation using convection-permitting simulations: A study over South America. Geophysical Research Letters, 52, e2025GL118806. https://doi.org/10.1029/2025GL118806
- Dominguez, F., Hu, H., & Martinez, J. A. (2020). Two-layer Dynamic Recycling Model (2L-DRM): Learning from moisture tracking models of different complexity. Journal of Hydrometeorology, 21(1), 3-16. https://doi.org/10.1175/JHM-D-19-0101.1
