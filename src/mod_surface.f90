! precipitation (PP) and evapotranspiration (ET) from the 2-D WRF fields, as the former
! python programs PP_preprocessing_2LDRM_fromRAINNC-C.py and ET_preprocessing_2LDRM.py
! - PP of an output step: accumulated rain at its end minus at its start
! - ET of an output step: sum over its files of the evapotranspiration rate times the
!   seconds per file: QFX, or LH / Lv(T) as the python program
! Only 2-D arrays are used, coarsened (block means, in double precision) as soon as they
! are read.
module mod_surface
    use netcdf
    use IO, only: nc_check
    implicit none
    private
    public :: read_accum_rain, add_evap_rate

contains

    ! total accumulated rain of a wrfout file (kg m-2), coarsened: the sum over vars_rain
    ! (e.g. RAINNC, RAINC) of rain + BUCKET_MM * bucket counter. A missing variable (RAINC
    ! in some convection-permitting outputs) counts as 0, as a missing bucket counter or
    ! BUCKET_MM attribute (bucket option not used in WRF). In double precision, since the
    ! accumulated rain can be much larger than the rain of one output step.
    subroutine read_accum_rain(file, vars_rain, vars_bucket, nx, ny, coarse_factor, rain, buffer, ibuffer)
        character(len=*), intent(in)  :: file
        character(len=*), intent(in)  :: vars_rain(:), vars_bucket(:)
        integer, intent(in)           :: nx, ny, coarse_factor
        real(8), intent(out)          :: rain(nx/coarse_factor, ny/coarse_factor)
        real(4), intent(out)          :: buffer(nx, ny)     ! work arrays
        integer, intent(out)          :: ibuffer(nx, ny)
        real(4)                       :: bucket_mm
        integer                       :: ncid, varid, i, j, iv, cf
        logical                       :: found

        cf = coarse_factor
        rain = 0
        found = .false.
        call nc_check(nf90_open(trim(file), nf90_nowrite, ncid), trim(file))
        if (nf90_get_att(ncid, nf90_global, "BUCKET_MM", bucket_mm) /= nf90_noerr) bucket_mm = 0
        do iv = 1, size(vars_rain)
            if (nf90_inq_varid(ncid, trim(vars_rain(iv)), varid) /= nf90_noerr) cycle
            found = .true.
            call nc_check(nf90_get_var(ncid, varid, buffer, start=[1, 1, 1], count=[nx, ny, 1]), &
                          trim(file) // ": " // trim(vars_rain(iv)))
            ibuffer = 0
            if (bucket_mm > 0) then
                if (nf90_inq_varid(ncid, trim(vars_bucket(iv)), varid) == nf90_noerr) &
                    call nc_check(nf90_get_var(ncid, varid, ibuffer, start=[1, 1, 1], count=[nx, ny, 1]), &
                                  trim(file) // ": " // trim(vars_bucket(iv)))
            end if
            do j = 1, ny / cf
                do i = 1, nx / cf
                    rain(i, j) = rain(i, j) + &
                        sum(real(buffer(1+cf*(i-1):cf*i, 1+cf*(j-1):cf*j), 8) + &
                            real(bucket_mm, 8) * ibuffer(1+cf*(i-1):cf*i, 1+cf*(j-1):cf*j)) / (cf * cf)
                end do
            end do
        end do
        call nc_check(nf90_close(ncid), trim(file))
        if (.not. found) then
            write(*,'(A)') trim(file) // ": none of the rain variables found (" // trim(vars_rain(1)) // &
                           ", ...). Use --no_pp_et to skip PP and ET."
            error stop
        end if
    end subroutine read_accum_rain

    ! adds the evapotranspiration rate (kg m-2 s-1) of an open wrfout file, coarsened, to
    ! et_sum, from (et_from):
    ! - QFX: the surface moisture flux, the water WRF evaporates (its LH is QFX times a
    !   constant latent heat of vaporization, 2.5104e6 J kg-1 in Noah-MP, 2.5e6 over water)
    ! - LH: LH / Lv, with Lv (J kg-1) from the temperature var_temp (Bringfelt, 1986, Test
    !   of a forest evapotranspiration model, SMHI), as the former python program;
    !   2-3 % more than QFX at 296-302 K
    ! temp is only needed with et_from LH (an unallocated array can be passed for it: it is
    ! not present)
    subroutine add_evap_rate(ncid, file, et_from, var_qfx, var_lh, var_temp, nx, ny, coarse_factor, et_sum, &
                             rate, temp)
        integer, intent(in)           :: ncid
        character(len=*), intent(in)  :: file, et_from, var_qfx, var_lh, var_temp
        integer, intent(in)           :: nx, ny, coarse_factor
        real(8), intent(inout)        :: et_sum(nx/coarse_factor, ny/coarse_factor)
        real(4), intent(out)          :: rate(nx, ny)             ! work arrays
        real(4), intent(out), optional :: temp(nx, ny)
        integer                       :: i, j, cf

        if (et_from /= "QFX" .and. .not. present(temp)) error stop "add_evap_rate: work array temp missing"
        cf = coarse_factor

        if (et_from == "QFX") then
            call read_2d(var_qfx, rate)
        else
            call read_2d(var_lh, rate)
            call read_2d(var_temp, temp)
            ! LH (W m-2) / Lv (J kg-1), in place
            rate = rate / (4185.5 * (751.78 - 0.5655 * temp))
        end if
        do j = 1, ny / cf
            do i = 1, nx / cf
                et_sum(i, j) = et_sum(i, j) + sum(real(rate(1+cf*(i-1):cf*i, 1+cf*(j-1):cf*j), 8)) / (cf * cf)
            end do
        end do

    contains

        subroutine read_2d(name, data)
            character(len=*), intent(in) :: name
            real(4), intent(out)         :: data(nx, ny)
            integer                      :: varid

            call nc_check(nf90_inq_varid(ncid, trim(name), varid), trim(file) // ": " // trim(name) // &
                          " (use --no_pp_et to skip PP and ET)")
            call nc_check(nf90_get_var(ncid, varid, data, start=[1, 1, 1], count=[nx, ny, 1]), &
                          trim(file) // ": " // trim(name))
        end subroutine read_2d
    end subroutine add_evap_rate

end module mod_surface
