! program that process by chunks, each chunk of files is tempraily agreggated and becomes only one time step in the output
! number of files by chunk is:  time_factor / (number of timesteps by file)
! output variables are precipitable water and U and V (moisture averaged velocities)

program main
    USE netcdf
    ! USE OMP_LIB
    use, intrinsic :: ieee_arithmetic
    use mod_transforming
    use IO
    use mpi 
    use FLAP
    use filesystem, only: stem      ! see other functions here https://github.com/scivision/fortran-filesystem/blob/main/API.md?plain=1

    implicit none
    ! include 'mpif.h'

    integer                                      :: nx, ny, nz, nt, nslabs
    character(len=150)                           :: path_inputs, path_outputs
    character(len=50)                            :: filename_list_files, basename
    character(len=200),dimension(:), allocatable :: files_all
    character(len=200),dimension(:), allocatable :: files_chunk
    character(len=200)                           :: file_list_files, file, file_output
    integer                                      :: n_vars_P
    character(len=20),dimension(:), allocatable  :: vars_P
    character(len=20)                            :: var_mixing_ratio, var_Uhat, var_Vhat, var_What
    integer                                      :: ncid, status, varid, i, ncid_out
    real(4), dimension(:,:,:,:), allocatable     :: data_P, data_Uhat, data_Vhat, data_What, data_Q, data_Uhat_stag,&
                                                     data_Vhat_stag, data_What_stag
    real(4), dimension(:,:,:,:,:), allocatable   :: data_Pvars
    real(4), dimension(:),allocatable            :: pressure_interslabs   ! in Pa
    real(4), dimension(:,:,:,:), allocatable     :: data_PW, data_U, data_V, data_What_inter, data_Q_inter
    integer                                      :: nx_coarse, ny_coarse, coarse_factor, time_factor
    real(4), dimension(:,:,:,:), allocatable     :: data_U_coarse, data_V_coarse, data_PW_coarse, data_Q_coarse, data_What_coarse
    real(4), dimension(:,:,:,:), allocatable     :: data_U_out, data_V_out, data_PW_out, data_Q_out, data_What_out
    real(4), dimension(:,:,:,:), allocatable     :: data_PW_temp  ! for extra time aggreagtion of PW
    integer                                      :: extra_time_factor_PW, time_factor_extra, i_extra  ! for extra time aggreagtion of PW
    character(len=50)                            :: var_PW_extra  ! for extra time aggreagtion of PW
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
    ! two next ones read via CLI
    ! time_factor = 3
    ! extra_time_factor_PW = 24 ! must be  multiple of time_factor, if same value as time_factor is provided, 
                        ! then no extra aggregation occurs
    rem_file = 0
    

    n_vars_P = 2
    allocate(vars_P(n_vars_P))
    vars_P(1) = "P"
    vars_P(2) = "PB"
    var_Uhat = "U"
    var_Vhat = "V"
    var_What = "W"
    var_mixing_ratio = "QVAPOR"

    ! allocate(pressure_interslabs(nslabs-1))
    ! pressure_interslabs(1) = 70000
    ! if (nslabs > 2) pressure_interslabs(2) = 60000


    vars_out_names = (/"PW","U ","V ","W ","Q "/)

    var_PW_extra = "PW24"  ! only in case  extra_time_factor_PW  is not time_factor




    !##########################################################################################!
    !!!!!!!################################## CLI ####################################!!!!!!!!!!
    !##########################################################################################!

    ! reading arguments from command line if any
    ! https://stackoverflow.com/questions/13843772/how-to-use-command-line-arguments-in-fortran ! no logner used
    ! nor using FLAP

    ! maybe create a module for this
    call cli%init(description = 'Program that processes WRF outputs and writes imput binary files for running 2L-DRM. &
                                    &You obatain binary files for: PW, etc')
    call cli%add(switch='--listfiles', &
                 switch_ab='-lf',    &
                 help='Text file containing the list of WRF output files to read',   &
                 required=.true.,   &
                 act='store',       &
                 error=error_cli)
    if (error_cli/=0) stop


    call cli%add(switch='--pathout', &
                 switch_ab='-po',    &
                 help='Path where outputs will be saved. The outputs files will take the name of the file listfiles &
                            &plus the variable name.',   &
                 required=.true.,   &
                 act='store',       &
                 error=error_cli)
    if (error_cli/=0) stop


    call cli%add(switch='--timefactor', &
                 switch_ab='-tf',    &
                 help='time factor to which the variables will be calculated',   &
                 required=.true.,   &
                 act='store',       &
                 error=error_cli)
    if (error_cli/=0) stop


    call cli%add(switch='--extratimefactor', &
        switch_ab='-etf',    &
        help='extra time factor to which precipitable water (PW) will be calculated',   &
        required=.true.,   &
        act='store',       &
        error=error_cli)
    if (error_cli/=0) stop


    call cli%add(switch='--coarse_factor', &
        switch_ab='-cf',    &
        help='coarsening factor',   &
        required=.false.,   &
        act='store',       &
        def="1",          &
        error=error_cli)
    if (error_cli/=0) stop

    call cli%add(switch='--remfile', &
        switch_ab='-rf',    &
        help='Should be 0, 1 or 2. sometime the number of files is a multiple of timefactor + 1, because the list of files &
                &start with 00:00:00 and ends with 00:00:00. Since the number of files must be a multiple of timefactor, &
                & this options allows to remove the first file (enter 1) or the last file (enter 2) ',   &
        required=.false.,   &
        act='store',       &
        def="0",          &
        error=error_cli)
    if (error_cli/=0) stop


    call cli%add(switch='--nx', &
        switch_ab='-nx',    &
        help='Number of grid cells in x direction',   &
        required=.true.,   &
        act='store',       &
        ! def="1",          &
        error=error_cli)
    if (error_cli/=0) stop

    call cli%add(switch='--ny', &
        switch_ab='-ny',    &
        help='Number of grid cells in y direction',   &
        required=.true.,   &
        act='store',       &
        ! def="1",          &
        error=error_cli)
    if (error_cli/=0) stop

    call cli%add(switch='--nz', &
        switch_ab='-nz',    &
        help='Number of grid cells in the vertical',   &
        required=.true.,   &
        act='store',       &
        ! def="1",          &
        error=error_cli)
    if (error_cli/=0) stop


    call cli%add(switch='--nslabs', &
        switch_ab='-ns',    &
        help='Number of slabs to interpolate to',   &
        required=.false.,   &
        act='store',       &
        def="2",          &
        error=error_cli)
    if (error_cli/=0) stop

    call cli%add(switch='--pressurelevels', &
        switch_ab='-pl',    &
        help='Pressure levels to interpolate to (nslanbs -1 pressure levels)',   &
        required=.false.,   &
        act='store',       &
        def="700",          &
        nargs="+", &
        error=error_cli)
    if (error_cli/=0) stop

    call cli%get(switch='-lf', val=file_list_files, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-po', val=path_outputs, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-tf', val=time_factor, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-etf', val=extra_time_factor_PW, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-cf', val=coarse_factor, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-rf', val=rem_file, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-nx', val=nx, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-ny', val=ny, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-nz', val=nz, error=error_cli)
    if (error_cli/=0) stop
    call cli%get(switch='-ns', val=nslabs, error=error_cli)
    if (error_cli/=0) stop


    allocate(pressure_interslabs(nslabs-1))
    pressure_interslabs(1) = 70000
    if (nslabs > 2) pressure_interslabs(2) = 60000

    call cli%get_varying(switch='-pl', val=pressure_interslabs, error=error_cli)
    if (error_cli/=0) stop
    print *, pressure_interslabs

    !##########################################################################################!
    !!!!!!!####################### PROCESS (dont touch) ##############################!!!!!!!!!!
    !##########################################################################################!

    ! geometry and number of files
    nx_coarse = nx / coarse_factor
    ny_coarse = ny / coarse_factor

    n_files_by_chunk = time_factor/nt
    

    print *, file_list_files
    call read_list_files(file_list_files, files_all, rem_file)
    n_all_chunks = size(files_all) / n_files_by_chunk
    ! print *,size(files_all)
    
    ! STOP

    ! getting basename for outputs from input filename
    ! ppos = scan(trim(filename_list_files),".", BACK= .true.)
    ! basename = filename_list_files
    ! if ( ppos > 0 )  basename = filename_list_files(1:ppos-1)

    basename = stem(file_list_files)
    print *, "BASEANAME = ", basename

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)
    

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

    allocate(files_chunk(n_files_by_chunk))

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

    allocate(data_U(nx,ny,nslabs,nt))
    allocate(data_V(nx,ny,nslabs,nt))
    allocate(data_PW(nx,ny,nslabs,nt))
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

    counter_chunk = 1
    n_chunks_max = ceiling(n_all_chunks*1.0/size_Of_Cluster)
    n_chunks_last_processor = mod(n_all_chunks-1,n_chunks_max) +1  ! -1 and +1 to account for the case when mod would be 0
    ! n_chunks = n_chunks_max  ! default, in loop it is updated for last processor


    open(11,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_"  // trim(vars_out_names(1)) // &
                    "_proc" // str(process_Rank), & 
                    status='REPLACE', access='stream')
    open(12,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_"  // trim(vars_out_names(2)) // &
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    open(13,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_"  // trim(vars_out_names(3)) // &
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    open(14,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_"  // trim(vars_out_names(4)) // & 
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    open(15,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_"  // trim(vars_out_names(5)) // & 
                    "_proc" // str(process_Rank), &
                    status='REPLACE', access='stream')
    ! read(10) PW
    ! close(10)

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


        do i_file = 1,size(files_chunk)
            file = files_chunk(i_file)
            write (*,"('Processor ', I0, ', working in chunk ',I0,', reading file ', A)") process_Rank, i_chunk, trim(file)

            
            status = nf90_open(path = trim(file), mode = nf90_nowrite, ncid = ncid)

            if (size(vars_P) == 1) then
                status = nf90_inq_varid(ncid,vars_P(1),varid)
                status = nf90_get_var(ncid,varid,data_P)
            else 
                do i = 1,size(vars_P)
                    status = nf90_inq_varid(ncid,vars_P(i),varid)
                    status = nf90_get_var(ncid,varid,data_Pvars(:,:,:,:,i))
                end do

            end if

            status = nf90_inq_varid(ncid,var_Uhat,varid)
            status = nf90_get_var(ncid,varid,data_Uhat_stag)

            status = nf90_inq_varid(ncid,var_Vhat,varid)
            status = nf90_get_var(ncid,varid,data_Vhat_stag)

            status = nf90_inq_varid(ncid,var_What,varid)
            status = nf90_get_var(ncid,varid,data_What_stag)

            status = nf90_inq_varid(ncid,var_mixing_ratio,varid)
            status = nf90_get_var(ncid,varid,data_Q)   ! actually is mixing ratio, but will be transformed to specific humidity Q
            ! data_Q = data_Q / (1+data_Q) ! transforming mixing ratio to Q (specific_humidity)

            

            status = nf90_close(ncid)

            data_P = sum(data_Pvars,dim = 5)
            call destagger(data_Uhat_stag, data_Uhat, nx, ny, nz, nt, 1)
            call destagger(data_Vhat_stag, data_Vhat, nx, ny, nz, nt, 2)
            call destagger(data_What_stag, data_What, nx, ny, nz, nt, 3)
            data_Q = data_Q / (1+data_Q) ! transforming mixing ratio to Q (specific_humidity)
            ! print *, "P(100,50,1,1) = ", data_P(100,50,1,1)

            ! calculate_slabs should be improved for the case when two or more pressure levels fall between two presure values of the data
            call calculate_slabs(data_Uhat, data_Vhat, data_Q, data_P, pressure_interslabs, nx, ny, nz, nslabs, nt, data_PW, &
                                    data_U, data_V, data_Q_inter)
            call DINTERP3DZ(data_What, data_What_inter, data_P, pressure_interslabs, nx, ny, nz, nslabs-1, &
                                    ieee_value(pressure_interslabs(1), ieee_quiet_nan))

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

        ! file_output = trim(path_outputs) // trim(basename) // "_" //  trim(vars_out_names(1)) // ".nc"
        ! status = nf90_create(file_output, 0, ncid_out)
        ! status = nf90_def_dim(ncid_out, "x", nx_coarse, ncdf_dim1_id)
        ! status = nf90_def_dim(ncid_out, "y", ny_coarse, ncdf_dim2_id)
        ! status = nf90_def_dim(ncid_out, "slab", nslabs, ncdf_dim3_id)
        ! status = nf90_def_dim(ncid_out, "time", nf90_unlimited, ncdf_dim4_id)
        ! if (sizeof(data_PW(1,1,1,1)) == 4) then
        !     status = nf90_def_dim(ncid_out,vars_out_names(1), NF90_FLOAT,
        !                             (/ncdf_dim1_id, ncdf_dim2_id, ncdf_dim3_id, ncdf_dim4_id/), ncdf_varout_id)
        ! else if (sizeof(data_PW(1,1,1,1)) == 8) then
        !     status = nf90_def_dim(ncid_out,vars_out_names(1), NF90_DOUBLE,&
        !                             (/ncdf_dim1_id, ncdf_dim2_id, ncdf_dim3_id, ncdf_dim4_id/), ncdf_varout_id)
        ! end if  
        ! status = nf90_put_var(ncid_out, ncdf_varout_id, data_PW_out, (/1,1,1,i_chunk/))



    end do
    ! $OMP END PARALLEL DO
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

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
    

    deallocate(files_chunk)


    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    if (process_Rank==0) then
        ! for extra time aggragation PW
        if (mod(extra_time_factor_PW, time_factor) == 0) then
            time_factor_extra = extra_time_factor_PW / time_factor
            if (time_factor_extra > 1) then
                i_extra = 1
                allocate(data_PW_temp(nx_coarse, ny_coarse, nslabs, time_factor_extra))
                open(201,  file = trim(path_outputs) // trim(basename) //"_" // trim(var_PW_extra) //".dat", &
                                status='REPLACE', action="write", access='stream')
            end if
        else
            print *, "extra_time_factor_PW must be a multiple of time_factor"
            stop
        end if
        ! END for extra time aggragation PW

        open(101,file = trim(path_outputs) // trim(basename) //"_" // trim(vars_out_names(1)) //".dat", &
                status = "REPLACE", action = "write", access='stream')
        open(102,file = trim(path_outputs) // trim(basename) //"_" // trim(vars_out_names(2)) //".dat", &
                status = "REPLACE", action = "write", access='stream')
        open(103,file = trim(path_outputs) // trim(basename) //"_" // trim(vars_out_names(3)) //".dat", & 
                status = "REPLACE", action = "write", access='stream')
        open(104,file = trim(path_outputs) // trim(basename) //"_" // trim(vars_out_names(4)) //".dat", & 
                status = "REPLACE", action = "write", access='stream')
        open(105,file = trim(path_outputs) // trim(basename) //"_" // trim(vars_out_names(5)) //".dat", & 
                status = "REPLACE", action = "write", access='stream')

        do i_proc = 0,(size_Of_Cluster-1)
            write (*, "('Reading results of processor ',I0, ' and writing to final files')") i_proc
            open(11,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_" // trim(vars_out_names(1)) // "_proc" // & 
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            open(12,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_" // trim(vars_out_names(2)) // "_proc" // &
                    str(i_proc), &
                    status='OLD', action="read", access='stream')
            open(13,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_" // trim(vars_out_names(3)) // "_proc" // &
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            open(14,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_" // trim(vars_out_names(4)) // "_proc" // &
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')
            open(15,  file = trim(path_outputs) // "tempdata_" //trim(basename) //"_" // trim(vars_out_names(5)) // "_proc" // &
                    str(i_proc), & 
                    status='OLD', action="read", access='stream')


            ! n_chunks = n_chunks_max
            ! if (i_proc == (size_Of_Cluster-1))    n_chunks = n_chunks_last_processor
                
            do i_chunk = 1,n_chunks(i_proc+1)
                read(11)   data_PW_out
                read(12)   data_U_out
                read(13)   data_V_out
                read(14)   data_What_out
                read(15)   data_Q_out

                write(101)   data_PW_out
                write(102)   data_U_out
                write(103)   data_V_out
                write(104)   data_What_out
                write(105)   data_Q_out

                ! for extra time aggragation PW
                if (time_factor_extra > 1) then
                    data_PW_temp(:,:,:,i_extra) = data_PW_out(:,:,:,1)
                    if (i_extra == time_factor_extra) then
                        write(201)  sum(data_PW_temp,dim = 4) / time_factor_extra
                        i_extra = 0
                    end if
                    i_extra = i_extra + 1
                end if
                ! END for extra time aggragation PW  
            end do


            ! vars_out_names = (/"PW","U ","V ","W ","Q "/)
            close(11, status = "delete")
            close(12, status = "delete")
            close(13, status = "delete")
            close(14, status = "delete")
            close(15, status = "delete")
        end do

        close(101)
        close(102)
        close(103)
        close(104)
        close(105)

        if (time_factor_extra > 1) then
            close(201)
        end if
    end if

    deallocate(data_U_out)      
    deallocate(data_V_out) 
    deallocate(data_PW_out) 
    deallocate(data_Q_out) 
    deallocate(data_What_out)   

    call MPI_FINALIZE(ierror)

    

    call cpu_time(stopTime)


    write(*, '(A, F8.4)') 'Elapsed time, hr : ',  (stopTime - startTime)/(60*60)

end program  main


