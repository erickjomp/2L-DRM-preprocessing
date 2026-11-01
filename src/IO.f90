module IO
    implicit none

    ! length of the buffers for paths and file names (the maximum path length in Linux)
    integer, parameter :: len_path = 4096
contains
    ! if not removing initial or last file, rem_file should equal to 0
    subroutine read_list_files(text_file, files, rem_file)
        ! https://fortran-lang.discourse.group/t/allocation-of-array-for-reading-strings-from-a-text-file/5986/4
        ! files is allocated with the length of the longest line, not len_path: every MPI
        ! process keeps the whole list, which can have tens of thousands of files
        integer :: iun, nrecs, istat, max_len
        character(len=*), intent(in)  :: text_file
        character(len=:), intent(out), allocatable :: files(:)
        integer, intent(in)               :: rem_file
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
        if (rem_file == 1 .or. rem_file == 2) then
            if (rem_file == 1) then
                files_temp = files(2:)
            else if (rem_file == 2) then
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
