! program that process by chunks, each chunk of files is tempraily agreggated and becomes only one time step in the output
! number of files by chunk is:  time_factor / (number of timesteps by file)
! output variables are precipitable water and U and V (moisture averaged velocities)

program main
    USE netcdf
    ! USE OMP_LIB
    use, intrinsic :: ieee_arithmetic
    use mod_transforming
    use IO
#ifdef USE_MPI
    use mpi
#endif
    use FLAP
    use stdlib_system, only: join_path, is_directory, make_directory_all
    use stdlib_error, only: state_type
    use stdlib_datetime, only: datetime_type, operator(-), total_seconds
    use mod_output, only: output_files, run_info
    use mod_surface, only: read_accum_rain, add_evap_rate

    implicit none
    ! include 'mpif.h'

    integer                                      :: nx, ny, nz, nt, nslabs
    character(len=len_path)                      :: path_outputs
    character(len=:), allocatable                :: basename
    character(len=:),dimension(:), allocatable   :: files_all
    character(len=:),dimension(:), allocatable   :: files_chunk
    character(len=len_path)                      :: file_list_files
    character(len=:), allocatable                :: file
    character(len=:), allocatable                :: prefix_temp  ! path_outputs + start of the temporary file names
    type(state_type)                             :: error_mkdir
    character(len=3)                             :: format              ! nc or bin
    character(len=5)                             :: split               ! none, month or year
    type(output_files)                           :: outputs
    type(run_info)                               :: info
    character(len=:), allocatable                :: prefix, period
    character(len=19)                            :: time0_str, time_start_str, time_end_str
    type(datetime_type)                          :: time0
    real(8)                                      :: time_step           ! hours since time0
    integer                                      :: i_step
    integer                                      :: n_vars_P
    character(len=20),dimension(:), allocatable  :: vars_P
    character(len=20)                            :: var_mixing_ratio, var_Uhat, var_Vhat, var_What
    integer                                      :: ncid, status, varid, i, ncid_out
    real(4), dimension(:,:,:,:), allocatable     :: data_P, data_Uhat, data_Vhat, data_What, data_Uhat_stag,&
                                                     data_Vhat_stag, data_What_stag
    real(4), dimension(:,:,:,:), target, allocatable :: data_Q
    real(4), dimension(:,:,:,:), pointer         :: data_mixingratio                                        
    real(4), dimension(:,:,:,:,:), allocatable   :: data_Pvars
    real(4), dimension(:),allocatable            :: pressure_interslabs   ! in Pa
    real(4), dimension(:,:,:,:), allocatable     :: data_PW, data_U, data_V, data_What_inter, data_Q_inter
    integer                                      :: nx_coarse, ny_coarse, coarse_factor, time_factor
    real(4), dimension(:,:,:,:), allocatable     :: data_U_coarse, data_V_coarse, data_PW_coarse, data_Q_coarse, data_What_coarse
    real(4), dimension(:,:,:,:), allocatable     :: data_U_out, data_V_out, data_PW_out, data_Q_out, data_What_out
    integer                                      :: extra_time_factor_PW  ! deprecated -etf, ignored
    real(4)                                      :: startTime,stopTime
    integer                                      :: i_file, n_files_by_chunk, i_chunk, n_all_chunks, counter_chunk
    character(len=50), dimension(5)              :: vars_out_names
    ! real(4), dimension(:), allocatable           :: o

    integer                                      :: process_Rank, size_Of_Cluster, ierror, n_chunks_max, n_chunks_last_processor,&
                                                    n_chunks_accum, i_proc  ! MPI variables
    integer, dimension(:), allocatable           :: n_chunks  !MPI_variables
    ! integer, dimension(:), allocatable           :: i_chunk_start, i_chunk_end
    integer                                      :: ppos, ncdf_dim1_id, ncdf_dim2_id, ncdf_dim3_id, ncdf_dim4_id, ncdf_varout_id

    ! for getting arguments from command line (only for filename_list_files, path_outputs and coarse_factor for now)
    type(command_line_interface)                 :: cli    ! Command Line Interface (CLI).
    ! character(99)                                :: string ! String value.
    integer                                      :: error_cli  ! Error trapping flag.
    integer                                      :: rem_file
    character(len=:), allocatable                :: removed_file   ! file removed by --remfile ('' if none)

    ! precipitation (PP) and evapotranspiration (ET)
    logical                                      :: with_pp_et
    character(len=20)                            :: vars_rain(2), vars_rain_bucket(2), var_QFX, var_LH, var_temp_ET
    character(len=3)                             :: et_from   ! QFX or LH
    real(4), dimension(:,:), allocatable         :: work_2d_a, work_2d_b, work_2d_coarse, data_PP_out, data_ET_out
    integer, dimension(:,:), allocatable         :: work_2d_int
    real(8), dimension(:,:), allocatable         :: et_sum, rain_start, rain_end
    logical                                      :: has_rain_start

    real(4), dimension(:,:,:,:), allocatable     :: data_pottemp 
    real(4), dimension(:,:,:,:), allocatable     :: data_pottemp_inter, data_temp_inter, data_mixingratio_inter, data_density_inter
    real(4), dimension(:,:,:,:), allocatable     :: data_PWflux, data_PWflux_coarse, data_PWflux_out, data_density_coarse, &
                                                    data_density_out
    character(len=20)                            :: var_PWflux, var_density, var_perturbation_pottemp ! from PWFLUX program
    logical                                      :: cum_per_timestep, save_density   ! from PWFLUX program
    real(4)                                      :: nseconds_in_timestep  ! seconds between WRF files, from their Times
    integer                                      :: nseconds_cli          ! --nseconds, 0 if not given
    real(4)                                      :: base_state_pottemp               ! from PWFLUX program
    integer                                      :: i_interslab                      ! from PWFLUX program
    real(4)                                      :: Rd, kappa, P0                    ! from PWFLUX program


    call cpu_time(startTime)

    !##########################################################################################!
    !!!!!!!########################## INPUTS (modify herre) ##########################!!!!!!!!!!
    !##########################################################################################!

    ! path_inputs= "/data/keeling/a/erickkc2/a/preprocessing_DL-DRM/SAAG2/lists_raw_files_by4months/"
    ! filename_list_files = "avasae.txt" !"2018_months09to12.txt"  ! it is overwritten in case of entering command line argument
    ! filename_list_files  can be taken from command line argument
    ! path_outputs = "/data/keeling/a/erickkc2/f/SAAG2/input_files_2LDRM_by4months__4km/"


    ! nx = 1471
    ! ny = 2027
    ! nz = 60
    nt = 1   ! number of raw timesteps per file
    ! nslabs = 2
    coarse_factor = 1  ! optionally read in CLI
    ! read via CLI
    ! time_factor = 3
    rem_file = 0
    

    n_vars_P = 2
    allocate(vars_P(n_vars_P))
    vars_P(1) = "P"
    vars_P(2) = "PB"
    var_Uhat = "U"
    var_Vhat = "V"
    var_What = "W"
    var_mixing_ratio = "QVAPOR"

    ! PP: accumulated rain (resolved and convective), and bucket counters (with the bucket
    ! size of the BUCKET_MM attribute); missing ones count as 0
    vars_rain = [character(len=20) :: "RAINNC", "RAINC"]
    vars_rain_bucket = [character(len=20) :: "I_RAINNC", "I_RAINC"]
    ! ET: surface moisture flux (kg m-2 s-1), or latent heat flux (W m-2) and temperature
    ! for the latent heat of vaporization (K) with --et_from LH
    var_QFX = "QFX"
    var_LH = "LH"
    var_temp_ET = "T2B"

    ! allocate(pressure_interslabs(nslabs-1))
    ! pressure_interslabs(1) = 70000
    ! if (nslabs > 2) pressure_interslabs(2) = 60000


    vars_out_names = (/"PW","U ","V ","W ","Q "/)


    !############################################################################################################!
    !#######################################  from PWFLUX program ###############################################!
    !############################################################################################################!
    var_perturbation_pottemp = "T"       ! in K
    base_state_pottemp = 300     ! https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/output.html

    var_PWflux = "PWflux"
    cum_per_timestep = .true.
    save_density = .true.
    var_density = "moistdensity"
    !############################################################################################################!
    !##################################### END from PWFLUX program ##############################################!
    !############################################################################################################!




    !##########################################################################################!
    !!!!!!!################################## CLI ####################################!!!!!!!!!!
    !##########################################################################################!

    ! reading arguments from command line if any
    ! https://stackoverflow.com/questions/13843772/how-to-use-command-line-arguments-in-fortran ! no logner used
    ! nor using FLAP

    ! maybe create a module for this
    ! call cli%init(help = "ave")

    call cli%init(description = 'Program that processes WRF outputs and writes the input files for running 2L-DRM: &
                                    &PW, U, V, PWflux, PP (precipitation) and ET (evapotranspiration), in one NetCDF &
                                    &file (default) or as binary files (--format bin). &
                                    &W, Q and moistdensity are written only for validation purposes.', &
                  help = "Usage: mpirun [mpirun_optional_arguments] (no mpirun needed for a serial build) ")
    call cli%add(switch='--listfiles', &
                 switch_ab='-lf',    &
                 help='Text file containing the list of WRF output files to read. &
                    &The number of files should be [n x timefactor + 1], where n is an integer. &
                    &For other options see the --remfile option. ',   &
                 required=.true.,   &
                 act='store',       &
                 error=error_cli)
    if (error_cli/=0) error stop


    call cli%add(switch='--pathout', &
                 switch_ab='-po',    &
                 help='Path where outputs will be saved. Created (with its parent directories) if it does not exist. &
                            &The output files are named after the file of --listfiles without extension: <name>.nc, &
                            &or <name>_<variable>.dat with --format bin; with --split the period is added, &
                            &e.g. <name>_2018-05.nc.',   &
                 required=.true.,   &
                 act='store',       &
                 error=error_cli)
    if (error_cli/=0) error stop


    call cli%add(switch='--timefactor', &
                 switch_ab='-tf',    &
                 help='time factor used for time coarsening. &
                      &Each timefactor files in the input files become one time step in the outputs. ',   &
                 required=.true.,   &
                 act='store',       &
                 error=error_cli)
    if (error_cli/=0) error stop


    call cli%add(switch='--extratimefactor', &
        switch_ab='-etf',    &
        help='Deprecated and ignored: the daily PW24 and PWflux24 outputs are no longer written. &
              &Will be removed in a future version.',   &
        required=.false.,   &
        act='store',       &
        def="0",          &
        error=error_cli)
    if (error_cli/=0) error stop


    call cli%add(switch='--coarse_factor', &
        switch_ab='-cf',    &
        help='Spatial coarsening factor. The default is 1 (no coarsening).',   &
        required=.false.,   &
        act='store',       &
        def="1",          &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--remfile', &
        switch_ab='-rf',    &
        help='By default the program expects the list of files to have [n x timefactor + 1] files, where n is an integer. &
            &(i.e. the first  and last file have time 00:00:00). &
            &Then the program ignores the last file so it process  [n x timefactor] files. &
            &If you prefer to ignore the first file instead, use option 1. &
            &If the number number of files in the list of files is actually [n x timefactor], use option 0.',    &   
        required=.false.,   &
        act='store',       &
        def="2",          &
        error=error_cli)
    if (error_cli/=0) error stop


    call cli%add(switch='--format', &
        switch_ab='-fmt',    &
        help='Format of the outputs. nc: one NetCDF file with all the variables. &
            &bin: one raw binary file per variable (float32, Fortran order), as read by 2L-DRM.',   &
        required=.false.,   &
        act='store',       &
        def="nc",          &
        choices="nc,bin",  &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--split', &
        switch_ab='-sp',    &
        help='Write a separate set of outputs for each month or year (from the WRF Times, the period is &
            &added to the file names). Every output step must be within one period.',   &
        required=.false.,   &
        act='store',       &
        def="none",          &
        choices="none,month,year",  &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--no_pp_et', &
        switch_ab='-noppet',    &
        help='Do not write PP and ET (for WRF files without RAINNC/RAINC, LH or T2B). &
            &PP of an output step is the accumulated rain at its end minus at its start, so it needs &
            &[n x timefactor + 1] files and --remfile 1 or 2: the removed file is used as the end (or start) &
            &of the last (or first) step. ET is the sum over the files of a step of the evapotranspiration &
            &rate (see --et_from).',   &
        required=.false.,   &
        act='store_true',       &
        def='.false.',          &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--et_from', &
        switch_ab='-etfrom',    &
        help='Evapotranspiration rate used for ET. QFX: the surface moisture flux, the water WRF evaporates. &
            &LH: LH / Lv(T2B), as the former python program; 2-3 % more than QFX, since WRF calculates LH &
            &as QFX times a constant latent heat of vaporization.',   &
        required=.false.,   &
        act='store',       &
        def="QFX",          &
        choices="QFX,LH",  &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--nseconds', &
        switch_ab='-nsec',    &
        help='Seconds between consecutive WRF files (used for PWflux and ET). By default it is calculated &
            &from the Times of the first two files; if given, it must agree with them.',   &
        required=.false.,   &
        act='store',       &
        def="0",          &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--nx', &
        switch_ab='-nx',    &
        help='Number of grid cells in x direction',   &
        required=.true.,   &
        act='store',       &
        ! def="1",          &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--ny', &
        switch_ab='-ny',    &
        help='Number of grid cells in y direction',   &
        required=.true.,   &
        act='store',       &
        ! def="1",          &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--nz', &
        switch_ab='-nz',    &
        help='Number of grid cells in the vertical direction',   &
        required=.true.,   &
        act='store',       &
        ! def="1",          &
        error=error_cli)
    if (error_cli/=0) error stop


    call cli%add(switch='--nslabs', &
        switch_ab='-ns',    &
        help='Number of slabs to interpolate to. By default it is 2.',   &
        required=.false.,   &
        act='store',       &
        def="2",          &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%add(switch='--pressurelevels', &
        switch_ab='-pl',    &
        help='Pressure levels to interpolate to (in Pa). There should be [nslanbs - 1] pressure levels.',   &
        required=.true.,   &
        act='store',       &
        ! def="70000",          &
        nargs="+", &
        error=error_cli)
    if (error_cli/=0) error stop

    call cli%get(switch='-lf', val=file_list_files, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-po', val=path_outputs, error=error_cli)
    if (error_cli/=0) error stop
    ! a value that fills the buffer was probably cut
    if (len_trim(file_list_files) == len_path .or. len_trim(path_outputs) == len_path) then
        write(*,'(A,I0,A)') "The paths given in -lf and -po must be shorter than ", len_path, " characters"
        error stop
    end if
    ! join_path adds the separator, so a trailing one is removed (but not from a root
    ! like / or C:\)
    i = len_trim(path_outputs)
    if (i > 1) then
        if (scan(path_outputs(i:i), "/\") == 1 .and. path_outputs(i-1:i-1) /= ":") path_outputs(i:) = " "
    end if
    call cli%get(switch='-tf', val=time_factor, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-etf', val=extra_time_factor_PW, error=error_cli)
    if (error_cli/=0) error stop
    if (cli%is_passed(switch='-etf')) print *, "Warning: -etf/--extratimefactor is deprecated and ignored &
                                               &(PW24 and PWflux24 are no longer written)"
    call cli%get(switch='-cf', val=coarse_factor, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-rf', val=rem_file, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-fmt', val=format, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-sp', val=split, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-noppet', val=with_pp_et, error=error_cli)
    if (error_cli/=0) error stop
    with_pp_et = .not. with_pp_et
    call cli%get(switch='-etfrom', val=et_from, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-nsec', val=nseconds_cli, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-nx', val=nx, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-ny', val=ny, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-nz', val=nz, error=error_cli)
    if (error_cli/=0) error stop
    call cli%get(switch='-ns', val=nslabs, error=error_cli)
    if (error_cli/=0) error stop


    allocate(pressure_interslabs(nslabs-1))
    ! pressure_interslabs(1) = 70000
    ! if (nslabs > 2) pressure_interslabs(2) = 60000

    call cli%get_varying(switch='-pl', val=pressure_interslabs, error=error_cli)
    if (error_cli/=0) error stop
    print *, pressure_interslabs

    !##########################################################################################!
    !!!!!!!####################### PROCESS (dont touch) ##############################!!!!!!!!!!
    !##########################################################################################!

    ! geometry and number of files
    nx_coarse = nx / coarse_factor
    ny_coarse = ny / coarse_factor

    n_files_by_chunk = time_factor/nt
    

    print *, file_list_files
    call read_list_files(file_list_files, files_all, rem_file, removed_file)
    if (with_pp_et .and. len(removed_file) == 0) then
        print *, "PP needs the accumulated rain at the start and end of each output step, so the list of files &
                 &must have [n x timefactor + 1] files and --remfile must be 1 or 2. Use --no_pp_et to skip PP and ET."
        error stop
    end if
    call get_nseconds()
    n_all_chunks = (size(files_all)) / n_files_by_chunk
    ! print *,size(files_all)
    
    if (MOD(size(files_all), n_files_by_chunk) /= 0) then
        print *, "The number of files is not compatible with the required timefactor&
        &. By default he program expects [n x timefactor + 1 ] files and ignores the last file (n is an integer). &
        & If only [n x timefactor] files are being provided use option 0 for --remfile argument."
        error stop
    end if
    ! getting basename for outputs from input filename
    ! ppos = scan(trim(filename_list_files),".", BACK= .true.)
    ! basename = filename_list_files
    ! if ( ppos > 0 )  basename = filename_list_files(1:ppos-1)

    basename = stem(file_list_files)
    print *, "BASEANAME = ", basename
    prefix_temp = join_path(path_outputs, "tempdata_" // basename // "_")

#ifdef USE_MPI
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)
#else
    size_Of_Cluster = 1
    process_Rank = 0
#endif

    ! the output directory (and its parents) is created if it does not exist; only by
    ! process 0, the others wait for it before opening their temporary files
    if (process_Rank == 0 .and. .not. is_directory(trim(path_outputs))) then
        call make_directory_all(trim(path_outputs), error_mkdir)
        if (error_mkdir%error()) then
            write(*,'(A)') "Cannot create output directory '" // trim(path_outputs) // "': " // &
                           error_mkdir%print_msg()
            error stop
        end if
    end if
#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
#endif


    ! n_chunks_max = ceiling(n_all_chunks*1.0/size_Of_Cluster)
    ! print *,"n chunks max : ", n_chunks_max
    ! allocate(i_chunk_start(size_Of_Cluster))
    ! allocate(i_chunk_end(size_Of_Cluster))
    ! i_chunk_start(1) = 1
    ! do i_task = 1,size_Of_Cluster-1
    !     i_chunk_end(i_task) = n_chunks_max * i_task
    !     i_chunk_start(i_task+1) = n_chunks_max * i_task + 1
    ! end do
    ! i_chunk_end(size_Of_Cluster) = n_all_chunks

    ! call MPI_Scatter(i_chunks_all, n_chunks, MPI_INT, i_chunks, n_chunks, MPI_INT, 0, MPI_COMM_WORLD, ierror)


    ! $OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(data_U_out,data_V_out,data_PW_out,data_Q_out,data_What_out,&
    ! $OMP                                       data_U_coarse,data_V_coarse,data_PW_coarse,data_Q_coarse,data_What_coarse,&
    ! $OMP                                       data_P,data_Uhat,data_Vhat,data_What,data_Q,&
    ! $OMP                                       data_U,data_V,data_PW,data_Q_inter,data_What_inter, files_chunk)
    
    ! $OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(files_all,list_files,nx,ny,nz,nt,nslabs,coarse_factor,time_factor,vars_P,&
    ! $OMP                                      var_Uhat,var_Vhat,var_What,var_mixing_ratio,pressure_interslabs,nx_coarse,ny_coarse,&
    ! $OMP                                      n_files_by_chunk, n_chunks, startTime1, stopTime1) &
    ! $OMP                                    PRIVATE(ncid)

    ! do i_chunk = i_chunk_start(process_Rank+1),i_chunk_end(process_Rank+1)

    allocate(character(len=len(files_all)) :: files_chunk(n_files_by_chunk))

    allocate(data_U_coarse(nx_coarse, ny_coarse, nslabs,time_factor))      
    allocate(data_V_coarse(nx_coarse, ny_coarse, nslabs,time_factor)) 
    allocate(data_PW_coarse(nx_coarse, ny_coarse, nslabs,time_factor)) 
    allocate(data_Q_coarse(nx_coarse, ny_coarse, nslabs-1,time_factor)) 
    allocate(data_What_coarse(nx_coarse, ny_coarse, nslabs-1,time_factor))     

    allocate(data_P(nx,ny,nz,nt))
    allocate(data_Uhat(nx,ny,nz,nt))
    allocate(data_Vhat(nx,ny,nz,nt))
    allocate(data_What(nx,ny,nz,nt))
    allocate(data_Q(nx,ny,nz,nt))
    allocate(data_pottemp(nx,ny,nz,nt)) ! from PWFLUX program

    allocate(data_U(nx,ny,nslabs,nt))
    allocate(data_V(nx,ny,nslabs,nt))
    allocate(data_PW(nx,ny,nslabs,nt))
    ! interpolated data (interslabs) 
    allocate(data_Q_inter(nx,ny,nslabs-1,nt))
    allocate(data_What_inter(nx,ny,nslabs-1,nt))

    ! raw data to store
    allocate(data_Pvars(nx,ny,nz,nt,size(vars_P)))
    allocate(data_Uhat_stag(nx+1,ny,nz,nt))
    allocate(data_Vhat_stag(nx,ny+1,nz,nt))
    allocate(data_What_stag(nx,ny,nz+1,nt))

    !  data for output
    allocate(data_U_out(nx_coarse, ny_coarse, nslabs,1))      
    allocate(data_V_out(nx_coarse, ny_coarse, nslabs,1)) 
    allocate(data_PW_out(nx_coarse, ny_coarse, nslabs, 1)) 
    allocate(data_Q_out(nx_coarse, ny_coarse, nslabs-1,1)) 
    allocate(data_What_out(nx_coarse, ny_coarse, nslabs-1,1))  
    !############################################################################################################!
    !#######################################  from PWFLUX program ###############################################!
    !############################################################################################################!
    ! interpolated data (interslabs) 
    allocate(data_pottemp_inter(nx,ny,nslabs-1,nt))      ! from PWFLUX program
    allocate(data_temp_inter(nx,ny,nslabs-1,nt))         ! from PWFLUX program
    allocate(data_mixingratio_inter(nx,ny,nslabs-1,nt))  ! from PWFLUX program
    ! allocate(data_What_inter(nx,ny,nslabs-1,nt))       ! from PWFLUX program
    allocate(data_PWflux(nx,ny,nslabs-1,nt))             ! from PWFLUX program
    allocate(data_density_inter(nx,ny,nslabs-1,nt))      ! from PWFLUX program
    ! for ouput
    allocate(data_PWflux_coarse(nx_coarse, ny_coarse, nslabs-1,time_factor))     
    allocate(data_PWflux_out(nx_coarse, ny_coarse, nslabs-1,1)) 
    if (save_density) then
        allocate(data_density_coarse(nx_coarse, ny_coarse, nslabs-1,time_factor))     
        allocate(data_density_out(nx_coarse, ny_coarse, nslabs-1,1)) 
    end if
    ! PP and ET: only 2-D arrays
    if (with_pp_et) then
        allocate(work_2d_a(nx, ny), work_2d_int(nx, ny))
        ! only when used, otherwise not present in add_evap_rate
        if (et_from == "LH") allocate(work_2d_b(nx, ny))
        if (coarse_factor > 1) allocate(work_2d_coarse(nx_coarse, ny_coarse))
        allocate(et_sum(nx_coarse, ny_coarse), rain_start(nx_coarse, ny_coarse), rain_end(nx_coarse, ny_coarse))
        allocate(data_PP_out(nx_coarse, ny_coarse), data_ET_out(nx_coarse, ny_coarse))
    else
        allocate(data_PP_out(1, 1), data_ET_out(1, 1))   ! not written
    end if
    has_rain_start = .false.
    !############################################################################################################!
    !##################################### END from PWFLUX program ##############################################!
    !############################################################################################################!

    
    counter_chunk = 1
    n_chunks_max = ceiling(n_all_chunks*1.0/size_Of_Cluster)
    n_chunks_last_processor = mod(n_all_chunks-1,n_chunks_max) +1  ! -1 and +1 to account for the case when mod would be 0
    ! n_chunks = n_chunks_max  ! default, in loop it is updated for last processor


    open(11,  file = prefix_temp  // trim(vars_out_names(1)) // &
                    "_proc" // str(process_Rank), & 
                    status='REPLACE', access='stream')
    open(12,  file = prefix_temp  // trim(vars_out_names(2)) // &
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    open(13,  file = prefix_temp  // trim(vars_out_names(3)) // &
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    open(14,  file = prefix_temp  // trim(vars_out_names(4)) // & 
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    open(15,  file = prefix_temp  // trim(vars_out_names(5)) // & 
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    
    !############################################################################################################!
    !#######################################  from PWFLUX program ###############################################!
    !############################################################################################################!
    open(21,  file = prefix_temp  // trim(var_PWflux) // &
    "_proc" // str(process_Rank), & 
    status='REPLACE', access='stream')
    if (save_density) then
        open(22,  file = prefix_temp  // trim(var_density) // &
        "_proc" // str(process_Rank), status='REPLACE', access='stream')
    end if
    if (with_pp_et) then
        open(23, file = prefix_temp // "PP_proc" // str(process_Rank), status='REPLACE', access='stream')
        open(24, file = prefix_temp // "ET_proc" // str(process_Rank), status='REPLACE', access='stream')
    end if
    !############################################################################################################!
    !##################################### END from PWFLUX program ##############################################!
    !############################################################################################################!



    allocate(n_chunks(size_Of_Cluster))
    n_chunks = n_all_chunks / size_Of_Cluster
    do i_proc = 1, mod(n_all_chunks, size_Of_Cluster)
        n_chunks(i_proc) = n_chunks(i_proc) + 1
    end do
    
    i_proc = 0
    n_chunks_accum = n_chunks(i_proc+1)
    do i_chunk = 1, n_all_chunks

        ! REMOVE ! if you dont want to run parallel
        ! if ( ((i_chunk*size_Of_Cluster - 1) / n_all_chunks) == process_Rank) then 

        if (i_chunk > n_chunks_accum) then
            i_proc = i_proc + 1
            n_chunks_accum = n_chunks_accum + n_chunks(i_proc+1)
        end if

        if (i_proc == process_Rank) then 
            write (*, "('Processor ', I0, ' started chunk ', I0, ' (', I0,'/',I0, ' chunks in this processor)  ')") &
                      process_Rank, i_chunk, counter_chunk, n_chunks(i_proc+1)
            counter_chunk = counter_chunk + 1
        else
            cycle
        end if
        ! end REMOVE

        files_chunk = files_all(1+n_files_by_chunk*(i_chunk-1):n_files_by_chunk*i_chunk)
        if (with_pp_et) et_sum = 0


        do i_file = 1,size(files_chunk)
            file = trim(files_chunk(i_file))
            write (*,"('Processor ', I0, ', working in chunk ',I0,', reading file ', A)") process_Rank, i_chunk, trim(file)

            
            call nc_check(nf90_open(path = trim(file), mode = nf90_nowrite, ncid = ncid), trim(file))

            if (size(vars_P) == 1) then
                call nc_check(nf90_inq_varid(ncid,vars_P(1),varid), trim(file) // ": " // trim(vars_P(1)))
                call nc_check(nf90_get_var(ncid,varid,data_P), trim(file))
            else 
                do i = 1,size(vars_P)
                    call nc_check(nf90_inq_varid(ncid,vars_P(i),varid), trim(file) // ": " // trim(vars_P(i)))
                    call nc_check(nf90_get_var(ncid,varid,data_Pvars(:,:,:,:,i)), trim(file))
                end do

            end if

            call nc_check(nf90_inq_varid(ncid,var_Uhat,varid), trim(file) // ": " // trim(var_Uhat))
            call nc_check(nf90_get_var(ncid,varid,data_Uhat_stag), trim(file))

            call nc_check(nf90_inq_varid(ncid,var_Vhat,varid), trim(file) // ": " // trim(var_Vhat))
            call nc_check(nf90_get_var(ncid,varid,data_Vhat_stag), trim(file))

            call nc_check(nf90_inq_varid(ncid,var_What,varid), trim(file) // ": " // trim(var_What))
            call nc_check(nf90_get_var(ncid,varid,data_What_stag), trim(file))

            call nc_check(nf90_inq_varid(ncid,var_mixing_ratio,varid), trim(file) // ": " // trim(var_mixing_ratio))
            ! actually is mixing ratio, but will be transformed to specific humidity Q
            call nc_check(nf90_get_var(ncid,varid,data_Q), trim(file))
            ! data_Q = data_Q / (1+data_Q) ! transforming mixing ratio to Q (specific_humidity)
            data_mixingratio => data_Q

            call nc_check(nf90_inq_varid(ncid,var_perturbation_pottemp,varid), &                ! from PWFLUX program
                          trim(file) // ": " // trim(var_perturbation_pottemp))
            call nc_check(nf90_get_var(ncid,varid,data_pottemp), trim(file))                    ! from PWFLUX program
            data_pottemp = data_pottemp + base_state_pottemp                   ! from PWFLUX program

            if (with_pp_et) call add_evap_rate(ncid, file, et_from, var_QFX, var_LH, var_temp_ET, nx, ny, &
                                               coarse_factor, et_sum, work_2d_a, work_2d_b, work_2d_coarse)

            call nc_check(nf90_close(ncid), trim(file))

            data_P = sum(data_Pvars,dim = 5)
            call destagger(data_Uhat_stag, data_Uhat, nx, ny, nz, nt, 1)
            call destagger(data_Vhat_stag, data_Vhat, nx, ny, nz, nt, 2)
            call destagger(data_What_stag, data_What, nx, ny, nz, nt, 3)

    !############################################################################################################!
    !#######################################  from PWFLUX program ###############################################!
    !############################################################################################################!
            ! interpolating to interslabs
            call DINTERP3DZ(data_mixingratio, data_mixingratio_inter, data_P, pressure_interslabs, nx, ny, nz, nslabs-1, &
                        ieee_value(pressure_interslabs(1), ieee_quiet_nan))
            call DINTERP3DZ(data_What, data_What_inter, data_P, pressure_interslabs, nx, ny, nz, nslabs-1, &
                        ieee_value(pressure_interslabs(1), ieee_quiet_nan))
            call DINTERP3DZ(data_pottemp, data_pottemp_inter, data_P, pressure_interslabs, nx, ny, nz, nslabs-1, &
                        ieee_value(pressure_interslabs(1), ieee_quiet_nan))

            kappa = 0.286   ! Rd/Cp  ! https://en.wikipedia.org/wiki/Potential_temperature  
            P0 = 100000.0 ! Pascales
            ! https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html
            ! data_temp = data_pottemp * (data_P / P0)**kappa

            ! transfort

            Rd = 287.052874
            ! eq 14 from http://www.atmo.arizona.edu/students/courselinks/fall11/atmo551a/ATMO_451a_551a_files/GasLawHydrostatic.pdf
            ! https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.density.html
            ! data_density = 0.622 * data_P * (1 + data_mixingratio) / (Rd * data_temp * (data_mixingratio + 0.622))
            do i_interslab = 1, nslabs-1
                data_temp_inter(:,:,i_interslab,:) = data_pottemp_inter(:,:,i_interslab,:) * &
                                                        (pressure_interslabs(i_interslab) / P0)** kappa
                data_density_inter(:,:,i_interslab,:) = &
                        0.622 * pressure_interslabs(i_interslab) * &
                        (1 + data_mixingratio_inter(:,:,i_interslab,:)) / &
                        (Rd * data_temp_inter(:,:,i_interslab,:) * (data_mixingratio_inter(:,:,i_interslab,:) + 0.622))
            end do

            ! data_PWflux = data_Q * data_density * data_What
            data_PWflux = (data_mixingratio_inter  / (1+data_mixingratio_inter))  * data_density_inter * data_What_inter

            if (coarse_factor > 1) then              
                call coarsen_xy_grid(data_PWflux,  data_PWflux_coarse(:,:,:,1+nt*(i_file-1):nt*i_file),   &
                                    nx, ny, nslabs-1, nt, coarse_factor)    
                if (save_density) then
                    call coarsen_xy_grid(data_density_inter,  data_density_coarse(:,:,:,1+nt*(i_file-1):nt*i_file),   &
                                    nx, ny, nslabs-1, nt, coarse_factor)  
                end if
            else
                data_PWflux_coarse(:,:,:,1+nt*(i_file-1):nt*i_file) = data_PWflux
                if (save_density) then
                    data_density_coarse(:,:,:,1+nt*(i_file-1):nt*i_file) = data_density_inter
                end if
            end if
    !############################################################################################################!
    !##################################### END from PWFLUX program ##############################################!
    !############################################################################################################!




            data_Q = data_mixingratio / (1+data_mixingratio) ! transforming mixing ratio to Q (specific_humidity)
            ! print *, "P(100,50,1,1) = ", data_P(100,50,1,1)

            ! calculate_slabs should be improved for the case when two or more pressure levels fall between two presure values of the data
            call calculate_slabs(data_Uhat, data_Vhat, data_Q, data_P, pressure_interslabs, nx, ny, nz, nslabs, nt, data_PW, &
                                    data_U, data_V, data_Q_inter)
            ! call DINTERP3DZ(data_What, data_What_inter, data_P, pressure_interslabs, nx, ny, nz, nslabs-1, &
            !                         ieee_value(pressure_interslabs(1), ieee_quiet_nan))

            if (coarse_factor > 1) then              
                call coarsen_xy_grid(data_U,          data_U_coarse(:,:,:,1+nt*(i_file-1):nt*i_file),   nx, ny, nslabs,   &
                                                                                                        nt, coarse_factor)
                call coarsen_xy_grid(data_V,          data_V_coarse(:,:,:,1+nt*(i_file-1):nt*i_file),   nx, ny, nslabs,   &
                                                                                                        nt, coarse_factor)
                call coarsen_xy_grid(data_PW,         data_PW_coarse(:,:,:,1+nt*(i_file-1):nt*i_file),  nx, ny, nslabs,   &
                                                                                                        nt, coarse_factor)
                call coarsen_xy_grid(data_Q_inter,    data_Q_coarse(:,:,:,1+nt*(i_file-1):nt*i_file),   nx, ny, nslabs-1, &
                                                                                                        nt, coarse_factor)
                call coarsen_xy_grid(data_What_inter, data_What_coarse(:,:,:,1+nt*(i_file-1):nt*i_file),nx, ny, nslabs-1, &
                                                                                                        nt, coarse_factor)
            else
                data_U_coarse(:,:,:,1+nt*(i_file-1):nt*i_file) = data_U
                data_V_coarse(:,:,:,1+nt*(i_file-1):nt*i_file) = data_V
                data_PW_coarse(:,:,:,1+nt*(i_file-1):nt*i_file) = data_PW
                data_Q_coarse(:,:,:,1+nt*(i_file-1):nt*i_file) = data_Q_inter
                data_What_coarse(:,:,:,1+nt*(i_file-1):nt*i_file) = data_What_inter
            end if

        end do  

        data_U_out(:,:,:,1) = sum(data_U_coarse,dim = 4)/time_factor
        data_V_out(:,:,:,1) = sum(data_V_coarse,dim = 4)/time_factor
        data_PW_out(:,:,:,1) = sum(data_PW_coarse,dim = 4)/time_factor
        data_Q_out(:,:,:,1) = sum(data_Q_coarse,dim = 4)/time_factor
        data_What_out(:,:,:,1) = sum(data_What_coarse,dim = 4)/time_factor

        write(11) data_PW_out 
        write(12) data_U_out 
        write(13) data_V_out 
        write(14) data_What_out 
        write(15) data_Q_out 

    !############################################################################################################!
    !#######################################  from PWFLUX program ###############################################!
    !############################################################################################################!
        data_PWflux_out(:,:,:,1) = sum(data_PWflux_coarse,dim = 4)/time_factor
        write(21) data_PWflux_out 

        if (save_density) then
            data_density_out(:,:,:,1)  = sum(data_density_coarse,dim = 4)/time_factor
            write(22)   data_density_out
        end if

        if (with_pp_et) then
            ! accumulated rain at the start of the step (the end of the previous one) and at its end
            if (.not. has_rain_start) then
                call read_accum_rain(boundary_file(n_files_by_chunk*(i_chunk-1)), vars_rain, vars_rain_bucket, &
                                     nx, ny, coarse_factor, rain_start, work_2d_a, work_2d_int)
                has_rain_start = .true.
            end if
            call read_accum_rain(boundary_file(n_files_by_chunk*i_chunk), vars_rain, vars_rain_bucket, &
                                 nx, ny, coarse_factor, rain_end, work_2d_a, work_2d_int)
            data_PP_out = real(rain_end - rain_start, 4)
            rain_start = rain_end
            data_ET_out = real(max(et_sum * nseconds_in_timestep, 0d0), 4)
            write(23) data_PP_out
            write(24) data_ET_out
        end if
    !############################################################################################################!
    !##################################### END from PWFLUX program ##############################################!
    !############################################################################################################!


    end do
    ! $OMP END PARALLEL DO
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(21)  ! from PWFLUX program
    close(22)  ! from PWFLUX program
    if (with_pp_et) then
        close(23)
        close(24)
        deallocate(work_2d_a, work_2d_int, et_sum, rain_start, rain_end)
        if (allocated(work_2d_b)) deallocate(work_2d_b)
        if (allocated(work_2d_coarse)) deallocate(work_2d_coarse)
    end if

    ! CLOSING raw data to store
    deallocate(data_Pvars)
    deallocate(data_Uhat_stag)
    deallocate(data_Vhat_stag)
    deallocate(data_What_stag)
    

    ! CLOSING transformed variables
    deallocate(data_P)
    deallocate(data_Uhat)
    deallocate(data_Vhat)
    deallocate(data_What)
    deallocate(data_Q)

        
    deallocate(data_U)
    deallocate(data_V)
    deallocate(data_PW)
    deallocate(data_Q_inter)
    deallocate(data_What_inter)  


    deallocate(data_U_coarse)      
    deallocate(data_V_coarse) 
    deallocate(data_PW_coarse) 
    deallocate(data_Q_coarse) 
    deallocate(data_What_coarse)  

    !############################################################################################################!
    !#######################################  from PWFLUX program ###############################################!
    !############################################################################################################!
    deallocate(data_pottemp)
        
    ! ! interpolated data (interslabs)
    deallocate(data_pottemp_inter)
    deallocate(data_temp_inter)
    ! deallocate(data_mixingratio_intZer)
    ! deallocate(data_What_inter)  
    deallocate(data_PWflux)
    deallocate(data_PWflux_coarse)  
    ! deallocate(data_PWflux_out)
    
    if (save_density) then
        deallocate(data_density_coarse)
    end if
    !############################################################################################################!
    !##################################### END from PWFLUX program ##############################################!
    !############################################################################################################!


    deallocate(files_chunk)


#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
#endif
    if (process_Rank==0) then
        ! run description for the NetCDF global attributes, and time origin of the outputs
        info%list_files = trim(file_list_files)
        info%split = trim(split)
        info%time_factor = time_factor
        info%coarse_factor = coarse_factor
        info%rem_file = rem_file
        info%et_from = trim(et_from)
        info%pressure_interslabs = pressure_interslabs
        time_start_str = ""
        if (format == "nc" .or. split /= "none") then
            time0_str = read_wrf_time(files_all(1))
            time0 = parse_datetime(time0_str)
            info%time_units = "hours since " // time0_str(1:10) // " " // time0_str(12:19)
        end if
        if (format == "nc") call read_lat_lon(files_all(1), nx, ny, coarse_factor, info%lat, info%lon)
        if (.not. allocated(data_density_out)) then
            allocate(data_density_out(nx_coarse, ny_coarse, nslabs-1, 1))
            data_density_out = ieee_value(data_density_out(1,1,1,1), ieee_quiet_nan)
        end if

        period = ""
        i_step = 0
        do i_proc = 0,(size_Of_Cluster-1)
            write (*, "('Reading results of processor ',I0, ' and writing to final files')") i_proc
            open(11,  file = prefix_temp // trim(vars_out_names(1)) // "_proc" // & 
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            open(12,  file = prefix_temp // trim(vars_out_names(2)) // "_proc" // &
                    str(i_proc), &
                    status='OLD', action="read", access='stream')
            open(13,  file = prefix_temp // trim(vars_out_names(3)) // "_proc" // &
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            open(14,  file = prefix_temp // trim(vars_out_names(4)) // "_proc" // &
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            open(15,  file = prefix_temp // trim(vars_out_names(5)) // "_proc" // &
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
    !############################################################################################################!
    !#######################################  from PWFLUX program ###############################################!
    !############################################################################################################!
            open(21,  file = prefix_temp // trim(var_PWflux) // "_proc" // & 
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            if (save_density) then
                open(22,  file = prefix_temp // trim(var_density) // "_proc" // & 
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            end if
            if (with_pp_et) then
                open(23, file = prefix_temp // "PP_proc" // str(i_proc), status='OLD', action="read", access='stream')
                open(24, file = prefix_temp // "ET_proc" // str(i_proc), status='OLD', action="read", access='stream')
            end if
    !############################################################################################################!
    !##################################### END from PWFLUX program ##############################################!
    !############################################################################################################!

            ! the processes have consecutive chunks, so reading them in order gives the output steps in order
            do i_chunk = 1,n_chunks(i_proc+1)
                i_step = i_step + 1
                read(11)   data_PW_out
                read(12)   data_U_out
                read(13)   data_V_out
                read(14)   data_What_out
                read(15)   data_Q_out
                read(21)   data_PWflux_out
                if (cum_per_timestep) data_PWflux_out = data_PWflux_out * (time_factor * nseconds_in_timestep)
                if (save_density) read(22) data_density_out
                if (with_pp_et) read(23) data_PP_out
                if (with_pp_et) read(24) data_ET_out

                ! times of the first and last file of the output step: the step is labelled with the first
                ! one and, when splitting, both must be in the same period
                time_step = 0
                if (format == "nc" .or. split /= "none") then
                    time_start_str = read_wrf_time(files_all(1 + n_files_by_chunk*(i_step-1)))
                    time_step = total_seconds(parse_datetime(time_start_str) - time0) / 3600
                    ! equally spaced files, as assumed for PWflux and ET
                    if (abs(time_step - (i_step-1) * time_factor * nseconds_in_timestep / 3600d0) > 1d-6) then
                        write(*,'(A,I0,A)') "Output step ", i_step, " starts at " // time_start_str // &
                            ", but the WRF files are not equally spaced (" // trim(str(nint(nseconds_in_timestep))) // &
                            " s between the first two files)"
                        error stop
                    end if
                end if
                if (split /= "none") then
                    time_end_str = read_wrf_time(files_all(n_files_by_chunk*i_step))
                    if (period_name(time_end_str, split) /= period_name(time_start_str, split)) then
                        write(*,'(A,I0,A)') "Output step ", i_step, " (" // time_start_str // " to " // &
                            time_end_str // ") is not within one " // trim(split) // ", so it cannot be split. &
                            &Every output step must start and end in the same " // trim(split) // &
                            " (check where the list of files starts and the -rf option)."
                        error stop
                    end if
                end if

                ! a new set of output files for each period
                if (.not. outputs%is_open .or. period_name(time_start_str, split) /= period) then
                    call outputs%close()
                    period = period_name(time_start_str, split)
                    if (split == "none") then
                        prefix = join_path(path_outputs, basename)
                    else
                        prefix = join_path(path_outputs, basename // "_" // period)
                    end if
                    write (*, "('Writing ', A)") prefix // trim(merge(".nc   ", "_*.dat", format == "nc"))
                    call outputs%create(format, prefix, nx_coarse, ny_coarse, nslabs, with_pp_et, info)
                end if
                call outputs%write_step(time_step, data_PW_out, data_U_out, data_V_out, data_What_out, &
                                        data_Q_out, data_PWflux_out, data_density_out, data_PP_out, data_ET_out)
            end do

            close(11, status = "delete")
            close(12, status = "delete")
            close(13, status = "delete")
            close(14, status = "delete")
            close(15, status = "delete")

            close(21, status = "delete")      ! from PWFLUX program
            if (save_density) close(22, status = "delete")      ! from PWFLUX program
            if (with_pp_et) then
                close(23, status = "delete")
                close(24, status = "delete")
            end if
        end do

        call outputs%close()
    end if

    deallocate(data_U_out)      
    deallocate(data_V_out) 
    deallocate(data_PW_out) 
    deallocate(data_Q_out) 
    deallocate(data_What_out)

#ifdef USE_MPI
    call MPI_FINALIZE(ierror)
#endif

    

    call cpu_time(stopTime)


    write(*, '(A, F8.4)') 'Elapsed time, hr : ',  (stopTime - startTime)/(60*60)

contains

    ! file j of the list before --remfile removed one (j = 0 .. n x timefactor): the
    ! instants that bound the output steps for PP (step i goes from file (i-1) x timefactor
    ! to file i x timefactor)
    function boundary_file(j) result(f)
        integer, intent(in)           :: j
        character(len=:), allocatable :: f
        integer                       :: i

        i = j + merge(0, 1, rem_file == 1)
        if (i < 1 .or. i > size(files_all)) then
            f = removed_file
        else
            f = trim(files_all(i))
        end if
    end function boundary_file

    ! seconds between WRF files, from the Times of the first two files (with the removed
    ! file if there is only one), checked against --nseconds if given
    subroutine get_nseconds()
        character(len=:), allocatable :: second_file
        real(8)                       :: seconds

        if (size(files_all) > 1) then
            second_file = trim(files_all(2))
        else
            second_file = removed_file
        end if
        if (len(second_file) == 0) then
            if (nseconds_cli <= 0) then
                print *, "Only one WRF file: give the seconds between WRF files with --nseconds"
                error stop
            end if
            nseconds_in_timestep = nseconds_cli
            return
        end if
        seconds = total_seconds(parse_datetime(read_wrf_time(second_file)) - &
                                parse_datetime(read_wrf_time(trim(files_all(1)))))
        seconds = abs(seconds)
        if (seconds <= 0) then
            print *, "The first two WRF files have the same time"
            error stop
        end if
        if (nseconds_cli > 0 .and. abs(seconds - nseconds_cli) > 1d-6) then
            write(*,'(A,I0,A,I0,A)') "--nseconds ", nseconds_cli, " does not agree with the WRF Times (", &
                                     nint(seconds), " s between the first two files)"
            error stop
        end if
        nseconds_in_timestep = real(seconds, 4)
    end subroutine get_nseconds

end program  main


