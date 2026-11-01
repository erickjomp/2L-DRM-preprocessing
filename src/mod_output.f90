! writes the output variables, one output step at a time, either as one NetCDF file
! (<prefix>.nc) or as one raw binary file per variable (<prefix>_<variable>.dat, float32
! in Fortran order (nx, ny, levels, time), as read by 2L-DRM)
module mod_output
    use netcdf
    use IO, only: nc_check
    implicit none
    private
    public :: output_files, run_info, nvars, var_names

    integer, parameter :: nvars = 7
    ! order of the arrays given to write_step
    character(len=12), parameter :: var_names(nvars) = &
        [character(len=12) :: "PW", "U", "V", "W", "Q", "PWflux", "moistdensity"]
    logical, parameter :: at_interslab(nvars) = [.false., .false., .false., .true., .true., .true., .true.]
    character(len=8), parameter :: var_units(nvars) = &
        [character(len=8) :: "kg m-2", "m s-1", "m s-1", "m s-1", "kg kg-1", "kg m-2", "kg m-3"]
    ! short labels, and the definitions of the Supporting Information of Claros et al. (2025), Text S2
    character(len=40), parameter :: var_long_names(nvars) = [character(len=40) :: &
        "precipitable water of the slab", &
        "moisture-weighted u wind of the slab", &
        "moisture-weighted v wind of the slab", &
        "vertical wind at the interslab level", &
        "specific humidity at the interslab level", &
        "moisture flux between slabs", &
        "moist air density at the interslab level"]
    character(len=140), parameter :: var_descriptions(nvars) = [character(len=140) :: &
        "water vapor contained in a unit area column of the slab (w in Claros et al., 2025)", &
        "moisture-weighted average of the wind velocity in the slab, along the grid x axis (west_east)", &
        "moisture-weighted average of the wind velocity in the slab, along the grid y axis (south_north)", &
        "vertical wind velocity interpolated to the interslab level", &
        "specific humidity interpolated to the interslab level", &
        "vertical moisture flux at the interslab level, accumulated over the output step: " // &
            "upward (Fu) if positive, downward (Fd) if negative", &
        "density of moist air at the interslab level"]

    ! description of the run, written as global attributes of the NetCDF files
    type :: run_info
        character(len=:), allocatable :: list_files, split, time_units
        integer                       :: time_factor, coarse_factor, rem_file
        real(4), allocatable          :: pressure_interslabs(:)  ! Pa
        real(4), allocatable          :: lat(:,:), lon(:,:)      ! (nx, ny)
    end type run_info

    type :: output_files
        character(len=3)              :: format = "nc"   ! nc or bin
        character(len=:), allocatable :: prefix
        integer                       :: nx, ny, nslabs, nsteps = 0
        integer                       :: ncid = -1, varid_time
        integer                       :: ids(nvars)       ! NetCDF variable ids or binary units
        logical                       :: is_open = .false.
    contains
        procedure :: create
        procedure :: write_step
        procedure :: close => close_files
    end type output_files

contains

    subroutine create(self, format, prefix, nx, ny, nslabs, info)
        class(output_files), intent(inout) :: self
        character(len=*), intent(in)       :: format, prefix
        integer, intent(in)                :: nx, ny, nslabs
        type(run_info), intent(in)         :: info
        integer                            :: i

        self%format = format
        self%prefix = prefix
        self%nx = nx
        self%ny = ny
        self%nslabs = nslabs
        self%nsteps = 0
        if (self%format == "bin") then
            do i = 1, nvars
                open(newunit=self%ids(i), file=self%prefix // "_" // trim(var_names(i)) // ".dat", &
                     status="REPLACE", action="write", access="stream")
            end do
        else
            call create_nc(self, info)
        end if
        self%is_open = .true.
    end subroutine create

    subroutine create_nc(self, info)
        class(output_files), intent(inout) :: self
        type(run_info), intent(in)         :: info
        character(len=*), parameter        :: coords = "XLONG XLAT"
        character(len=:), allocatable      :: file, slab_bounds
        character(len=20)                  :: buffer
        character(len=8)                   :: date
        character(len=10)                  :: clock
        integer :: dim_x, dim_y, dim_slab, dim_interslab, dim_time, varid_slab, varid_interslab, &
                   varid_lat, varid_lon, i, nlev

        file = self%prefix // ".nc"
        call nc_check(nf90_create(file, ior(nf90_clobber, nf90_netcdf4), self%ncid), file)
        call nc_check(nf90_def_dim(self%ncid, "west_east", self%nx, dim_x), file)
        call nc_check(nf90_def_dim(self%ncid, "south_north", self%ny, dim_y), file)
        call nc_check(nf90_def_dim(self%ncid, "slab", self%nslabs, dim_slab), file)
        call nc_check(nf90_def_dim(self%ncid, "interslab", self%nslabs - 1, dim_interslab), file)
        call nc_check(nf90_def_dim(self%ncid, "time", nf90_unlimited, dim_time), file)

        call nc_check(nf90_def_var(self%ncid, "time", nf90_double, [dim_time], self%varid_time), file)
        call put_att(self%ncid, self%varid_time, "long_name", "start of the output step")
        call put_att(self%ncid, self%varid_time, "units", info%time_units)
        call put_att(self%ncid, self%varid_time, "calendar", "standard")
        write(buffer, '(I0)') info%time_factor
        call put_att(self%ncid, self%varid_time, "comment", "each output step averages the " // trim(buffer) // &
                     " input files that start at this time")

        ! slab bounds: surface, interslab levels, model top
        slab_bounds = ""
        do i = 1, self%nslabs
            write(buffer, '(I0)') i
            slab_bounds = slab_bounds // "slab " // trim(buffer) // ": " // level_name(i - 1) // " to " // &
                          level_name(i) // merge("; ", "  ", i < self%nslabs)
        end do
        call nc_check(nf90_def_var(self%ncid, "slab", nf90_int, [dim_slab], varid_slab), file)
        call put_att(self%ncid, varid_slab, "long_name", "slab number")
        call put_att(self%ncid, varid_slab, "description", "slab number, from the surface up (1: lower slab)")
        call put_att(self%ncid, varid_slab, "bounds_description", trim(slab_bounds))
        call nc_check(nf90_def_var(self%ncid, "interslab", nf90_float, [dim_interslab], varid_interslab), file)
        call put_att(self%ncid, varid_interslab, "long_name", "interslab level pressure")
        call put_att(self%ncid, varid_interslab, "description", &
                     "pressure of the inter-slab interface between slab n and slab n+1")
        call put_att(self%ncid, varid_interslab, "units", "Pa")
        call put_att(self%ncid, varid_interslab, "positive", "down")

        call nc_check(nf90_def_var(self%ncid, "XLAT", nf90_float, [dim_x, dim_y], varid_lat), file)
        call put_att(self%ncid, varid_lat, "long_name", "latitude")
        call put_att(self%ncid, varid_lat, "units", "degree_north")
        call nc_check(nf90_def_var(self%ncid, "XLONG", nf90_float, [dim_x, dim_y], varid_lon), file)
        call put_att(self%ncid, varid_lon, "long_name", "longitude")
        call put_att(self%ncid, varid_lon, "units", "degree_east")

        do i = 1, nvars
            nlev = merge(self%nslabs - 1, self%nslabs, at_interslab(i))
            call nc_check(nf90_def_var(self%ncid, trim(var_names(i)), nf90_float, &
                          [dim_x, dim_y, merge(dim_interslab, dim_slab, at_interslab(i)), dim_time], self%ids(i), &
                          chunksizes=[self%nx, self%ny, nlev, 1], shuffle=.true., deflate_level=1), file)
            call put_att(self%ncid, self%ids(i), "long_name", trim(var_long_names(i)))
            call put_att(self%ncid, self%ids(i), "description", trim(var_descriptions(i)))
            call put_att(self%ncid, self%ids(i), "units", trim(var_units(i)))
            call put_att(self%ncid, self%ids(i), "coordinates", coords)
            if (at_interslab(i)) then
                call put_att(self%ncid, self%ids(i), "comment", "NaN where the interslab level is below the ground")
            else
                call put_att(self%ncid, self%ids(i), "comment", "NaN where the whole slab is below the ground")
            end if
        end do

        call put_att(self%ncid, nf90_global, "title", "2L-DRM inputs: slab variables")
        call put_att(self%ncid, nf90_global, "source", &
                     "prepro2LDRM-slabvars, https://github.com/erickjomp/2L-DRM-preprocessing")
        call put_att(self%ncid, nf90_global, "references", &
                     "Claros, E., Dominguez, F., & Liu, C. (2025), Geophysical Research Letters, 52, e2025GL118806, " // &
                     "https://doi.org/10.1029/2025GL118806; Dominguez, F., Hu, H., & Martinez, J. A. (2020), " // &
                     "Journal of Hydrometeorology, 21(1), 3-16, https://doi.org/10.1175/JHM-D-19-0101.1")
        call put_att(self%ncid, nf90_global, "list_files", info%list_files)
        call nc_check(nf90_put_att(self%ncid, nf90_global, "time_factor", info%time_factor), file)
        call nc_check(nf90_put_att(self%ncid, nf90_global, "coarse_factor", info%coarse_factor), file)
        call nc_check(nf90_put_att(self%ncid, nf90_global, "remfile", info%rem_file), file)
        call nc_check(nf90_put_att(self%ncid, nf90_global, "pressure_levels_Pa", info%pressure_interslabs), file)
        call put_att(self%ncid, nf90_global, "split", info%split)
        call date_and_time(date, clock)
        call put_att(self%ncid, nf90_global, "history", "created " // date(1:4) // "-" // date(5:6) // "-" // &
                     date(7:8) // " " // clock(1:2) // ":" // clock(3:4) // ":" // clock(5:6))
        call nc_check(nf90_enddef(self%ncid), file)

        call nc_check(nf90_put_var(self%ncid, varid_slab, [(i, i = 1, self%nslabs)]), file)
        call nc_check(nf90_put_var(self%ncid, varid_interslab, info%pressure_interslabs), file)
        call nc_check(nf90_put_var(self%ncid, varid_lat, info%lat), file)
        call nc_check(nf90_put_var(self%ncid, varid_lon, info%lon), file)

    contains

        ! name of the bottom of slab ilev + 1 (0: surface, nslabs: model top)
        function level_name(ilev) result(name)
            integer, intent(in)           :: ilev
            character(len=:), allocatable :: name
            character(len=20)             :: num

            if (ilev == 0) then
                name = "surface"
            else if (ilev == self%nslabs) then
                name = "model top"
            else
                write(num, '(I0)') nint(info%pressure_interslabs(ilev))
                name = trim(num) // " Pa"
            end if
        end function level_name
    end subroutine create_nc

    subroutine put_att(ncid, varid, name, value)
        integer, intent(in)          :: ncid, varid
        character(len=*), intent(in) :: name, value

        call nc_check(nf90_put_att(ncid, varid, name, value), "attribute " // name)
    end subroutine put_att

    ! appends one output step; the arrays are (nx, ny, nslabs or nslabs-1, 1)
    subroutine write_step(self, time, PW, U, V, W, Q, PWflux, density)
        class(output_files), intent(inout)      :: self
        real(8), intent(in)                     :: time  ! in the units of run_info%time_units
        real(4), dimension(:,:,:,:), intent(in) :: PW, U, V, W, Q, PWflux, density

        self%nsteps = self%nsteps + 1
        if (self%format == "bin") then
            write(self%ids(1)) PW
            write(self%ids(2)) U
            write(self%ids(3)) V
            write(self%ids(4)) W
            write(self%ids(5)) Q
            write(self%ids(6)) PWflux
            write(self%ids(7)) density
        else
            call nc_check(nf90_put_var(self%ncid, self%varid_time, [time], start=[self%nsteps]), self%prefix)
            call put_step(1, PW)
            call put_step(2, U)
            call put_step(3, V)
            call put_step(4, W)
            call put_step(5, Q)
            call put_step(6, PWflux)
            call put_step(7, density)
        end if

    contains

        subroutine put_step(i, data)
            integer, intent(in)                     :: i
            real(4), dimension(:,:,:,:), intent(in) :: data

            call nc_check(nf90_put_var(self%ncid, self%ids(i), data, start=[1, 1, 1, self%nsteps]), &
                          self%prefix // ".nc: " // trim(var_names(i)))
        end subroutine put_step
    end subroutine write_step

    subroutine close_files(self)
        class(output_files), intent(inout) :: self
        integer                            :: i

        if (.not. self%is_open) return
        if (self%format == "bin") then
            do i = 1, nvars
                close(self%ids(i))
            end do
        else
            call nc_check(nf90_close(self%ncid), self%prefix // ".nc")
        end if
        self%is_open = .false.
    end subroutine close_files

end module mod_output
