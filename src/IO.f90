module IO
    implicit none

    ! length of the buffers for paths and file names (the maximum path length in Linux)
    integer, parameter :: len_path = 4096
contains
    ! if not removing initial or last file, rem_file should equal to 0; removed is the
    ! removed file ('' if none)
    subroutine read_list_files(text_file, files, rem_file, removed)
        ! https://fortran-lang.discourse.group/t/allocation-of-array-for-reading-strings-from-a-text-file/5986/4
        ! files is allocated with the length of the longest line, not len_path: every MPI
        ! process keeps the whole list, which can have tens of thousands of files
        integer :: iun, nrecs, istat, max_len
        character(len=*), intent(in)  :: text_file
        character(len=:), intent(out), allocatable :: files(:)
        integer, intent(in)               :: rem_file
        character(len=:), allocatable, intent(out), optional :: removed
        character(len=len_path) :: buffer
        integer          :: i
        character(len=:), allocatable :: files_temp(:)

        open(newunit=iun, file=text_file, form="formatted", status = "old")
        ! count number of records and find the longest one
        nrecs = 0
        max_len = 1
        count_loop: do
            read(iun, '(a)', iostat=istat) buffer
            if (istat /= 0) exit count_loop
            nrecs = nrecs+1
            ! a line that fills the buffer was probably cut
            if (len_trim(buffer) == len_path) then
                write(*,'(A,I0,A,I0,A)') "Line ", nrecs, " of " // trim(text_file) // " has ", len_path, &
                                         " characters or more"
                error stop
            end if
            max_len = max(max_len, len_trim(buffer))
        end do count_loop

        allocate(character(len=max_len) :: files(nrecs))
        rewind(iun)
        read_loop: do i=1,nrecs
            read(iun,'(a)', iostat=istat) buffer
            ! print *,buffer
            if (istat /= 0) then
                print *,' error while reading text file'
                exit read_loop
            end if
            files(i) = trim(buffer)
        end do read_loop
        close(iun)

        ! only in case removing file
        if (present(removed)) removed = ""
        if ((rem_file == 1 .or. rem_file == 2) .and. nrecs > 0) then
            if (rem_file == 1) then
                if (present(removed)) removed = trim(files(1))
                files_temp = files(2:)
            else if (rem_file == 2) then
                if (present(removed)) removed = trim(files(nrecs))
                files_temp = files(:size(files)-1)
            end if
            call move_alloc(files_temp, files)
        end if
    end subroutine

    ! file name without directory and extension; both / and backslash are directory
    ! separators, so that it works with any path on Windows. Not base_name of
    ! stdlib_system: it only splits on the OS separator and returns '' for a name
    ! without directory (no backslash at the end of a comment line: cpp would join
    ! it with the next line)
    function stem(path) result(base)
        character(len=*), intent(in)  :: path
        character(len=:), allocatable :: base
        integer                       :: i_sep, i_dot

        i_sep = scan(trim(path), "/\", back=.true.)
        base = trim(path(i_sep+1:))
        i_dot = index(base, ".", back=.true.)
        if (i_dot > 1) base = base(:i_dot-1)
    end function stem

    ! first time of a wrfout file (variable Times), as "YYYY-MM-DD_hh:mm:ss"
    function read_wrf_time(file) result(time_str)
        use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_nowrite
        character(len=*), intent(in) :: file
        character(len=19)            :: time_str
        integer                      :: ncid, varid

        call nc_check(nf90_open(trim(file), nf90_nowrite, ncid), trim(file))
        call nc_check(nf90_inq_varid(ncid, "Times", varid), trim(file) // ": Times")
        call nc_check(nf90_get_var(ncid, varid, time_str, start=[1, 1], count=[19, 1]), trim(file) // ": Times")
        call nc_check(nf90_close(ncid), trim(file))
    end function read_wrf_time

    ! parses "YYYY-MM-DD_hh:mm:ss" (WRF Times), as in 2L-DRM
    function parse_datetime(str) result(dt)
        use stdlib_datetime, only: datetime_type, datetime
        character(len=*), intent(in) :: str
        type(datetime_type)          :: dt
        character(len=len(str))      :: buffer
        integer                      :: year, month, day, hour, minute, second, ios, i

        buffer = str
        do i = 1, len(buffer)
            if (index("-_:", buffer(i:i)) > 0) buffer(i:i) = " "
        end do
        read(buffer, *, iostat=ios) year, month, day, hour, minute, second
        if (ios /= 0) then
            write(*,*) "Invalid datetime '", trim(str), "'. Expected format: YYYY-MM-DD_hh:mm:ss"
            error stop
        end if
        dt = datetime(year, month, day, hour, minute, second)
    end function

    ! period of a time "YYYY-MM-DD_hh:mm:ss" for --split: "YYYY-MM", "YYYY" or ""
    function period_name(time_str, split) result(period)
        character(len=*), intent(in)  :: time_str, split
        character(len=:), allocatable :: period

        select case (trim(split))
        case ("month")
            period = time_str(1:7)
        case ("year")
            period = time_str(1:4)
        case default
            period = ""
        end select
    end function period_name

    ! XLAT and XLONG of a wrfout file, coarsened like the outputs; lat and lon are
    ! allocated as (nx/coarse_factor, ny/coarse_factor)
    subroutine read_lat_lon(file, nx, ny, coarse_factor, lat, lon)
        use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_nowrite
        use mod_transforming, only: coarsen_xy_grid
        character(len=*), intent(in)         :: file
        integer, intent(in)                  :: nx, ny, coarse_factor
        real(4), allocatable, intent(out)    :: lat(:,:), lon(:,:)
        real(4)                              :: fine(nx, ny, 1, 1), coarse(nx/coarse_factor, ny/coarse_factor, 1, 1)
        character(len=5), parameter          :: names(2) = ["XLAT ", "XLONG"]
        integer                              :: ncid, varid, i

        call nc_check(nf90_open(trim(file), nf90_nowrite, ncid), trim(file))
        do i = 1, 2
            call nc_check(nf90_inq_varid(ncid, trim(names(i)), varid), trim(file) // ": " // trim(names(i)))
            call nc_check(nf90_get_var(ncid, varid, fine, start=[1, 1, 1], count=[nx, ny, 1]), &
                          trim(file) // ": " // trim(names(i)))
            call coarsen_xy_grid(fine, coarse, nx, ny, 1, 1, coarse_factor)
            if (i == 1) lat = coarse(:,:,1,1)
            if (i == 2) lon = coarse(:,:,1,1)
        end do
        call nc_check(nf90_close(ncid), trim(file))
    end subroutine read_lat_lon

    ! stops with the NetCDF error message if status is not nf90_noerr
    subroutine nc_check(status, context)
        use netcdf, only: nf90_noerr, nf90_strerror
        integer, intent(in)          :: status
        character(len=*), intent(in) :: context

        if (status /= nf90_noerr) then
            write(*,'(A)') "NetCDF error (" // context // "): " // trim(nf90_strerror(status))
            error stop
        end if
    end subroutine nc_check

end module IO
