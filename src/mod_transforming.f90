module mod_transforming
    implicit none
contains

    subroutine destagger(data_stag, data_out, nx,ny,nz,nt,dim) 
        real(4), dimension(:,:,:,:),intent(in)        :: data_stag
        real(4), dimension(nx,ny,nz,nt),intent(out)        :: data_out
        integer,intent(in)                            :: nx,ny,nz,nt,dim
        integer                                       :: i,j,k,i_t
        if (dim == 1) then
            data_out = (data_stag(1:nx,:,:,:) + data_stag(2:nx+1,:,:,:)) /2 
            ! do i_t = 1,nt
            !     do k = 1,nz
            !         do j = 1,ny
            !             do i = 1, nx

            !             end do
            !         end do
            !     end do
            ! end do
        else if (dim == 2) then
            data_out = (data_stag(:,1:ny,:,:) + data_stag(:,2:ny+1,:,:)) /2 
        else if (dim == 3) then
            data_out = (data_stag(:,:,1:nz,:) + data_stag(:,:,2:nz+1,:)) /2 
        end if
    end subroutine

    subroutine coarsen_xy_grid(data_in, data_out,nx, ny, nz, nt, coarse_factor)
        use, intrinsic :: ieee_arithmetic
        real(4), dimension(nx,ny,nz,nt), intent(in)                              :: data_in
        real(4), dimension(nx/coarse_factor,ny/coarse_factor,nz,nt), intent(out) :: data_out
        integer, intent(in)                                                      :: nx, ny,nz, nt, coarse_factor
        integer                                                                  :: i,j,k,i_t
        integer                                                                  :: nx_out, ny_out
        
        nx_out = nx / coarse_factor
        ny_out = ny / coarse_factor

        do i_t = 1, nt
            do k = 1, nz
                do j = 1,ny_out
                    do i = 1, nx_out
                        data_out(i,j,k,i_t) = sum(data_in((1+coarse_factor*(i-1)):(coarse_factor*i), &
                                                            (1+coarse_factor*(j-1)):(coarse_factor*j),&
                                                            k,i_t)) / (coarse_factor * coarse_factor)
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine calculate_interslabs(What_wind, T, P, pressure_interslabs, nx, ny, nz, nslabs, n_timesteps, Wwind_interslabs)

        use, intrinsic :: ieee_arithmetic
        ! use iso_fortran_env, only: real32, real64
        implicit none
    
        real(kind=4), dimension(nx,ny,nz,n_timesteps), intent(in)               :: What_wind, T, P
        real(kind=4), dimension(nslabs-1), intent(in)                           :: pressure_interslabs
        integer, intent(in)                                                     :: nx,ny,nz,nslabs,n_timesteps
        real(kind=4), dimension(nx,ny,nslabs-1,n_timesteps), intent(out)        :: Wwind_interslabs
    
        integer                                                                 :: itime
    
        do itime = 1,n_timesteps
            call calculate_interslabs0(What_wind(:,:,:,itime), T(:,:,:,itime), P(:,:,:,itime), pressure_interslabs, &
                                       nx,ny,nz,nslabs, Wwind_interslabs(:,:,:,itime))
        end do
    
    end subroutine
    
    
    
    subroutine calculate_interslabs0(What_wind, T, P, pressure_interslabs, nx, ny, nz, nslabs, Wwind_interslabs)
    
        use, intrinsic :: ieee_arithmetic
        ! use iso_fortran_env, only: real32, real64
    
        implicit none
    
        real(kind=4), dimension(nx,ny,nz), intent(in)               :: What_wind, T, P
        real(kind=4), dimension(nslabs-1), intent(in)               :: pressure_interslabs
        integer, intent(in)                                         :: nx,ny,nz,nslabs
        real(kind=4), dimension(nx,ny,nslabs-1), intent(out)        :: Wwind_interslabs
    
        ! integer                                                     :: i,j,k,i_slab
        integer                                                     :: i_slab
        real(kind=4), dimension(nx,ny,nslabs-1)                     :: What_wind_interslabs, T_interslabs
        real(kind=4)                                                :: What_wind_1
        real(kind=4)                                                :: g,Rg
    
        g = 9.81
        Rg = 287.058
    
        call DINTERP3DZ(What_wind, What_wind_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, &
                        ieee_value(What_wind_1, ieee_quiet_nan))  ! output is U_interslabs
        call DINTERP3DZ(T, T_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, &
                        ieee_value(What_wind_1, ieee_quiet_nan))  ! output is U_interslabs
    
        do i_slab = 1,nslabs-1
            Wwind_interslabs(:,:,i_slab) = What_wind_interslabs(:,:,i_slab) * T_interslabs(:,:,i_slab) * &
                                           Rg / (-pressure_interslabs(i_slab) * g)  
        end do
    
    end subroutine

    ! calculate_slabs should be improved for the case when two or more pressure levels fall between two presure values of the data
    
    ! subroutine calculate_slabs(Uhat, Vhat, Q, P, pressure_interslabs, nx, ny, nz, nslabs, n_timesteps, missingval, PW_slabs, &
    !                            U_slabs, V_slabs)
    subroutine calculate_slabs(Uhat, Vhat, Q, P, pressure_interslabs, nx, ny, nz, nslabs, n_timesteps, PW_slabs, &
        U_slabs, V_slabs, Q_interslabs)
    
        ! use iso_fortran_env, only: real32, real64
    
        implicit none
        integer, intent(in)                                               :: nx, ny, nz, nslabs, n_timesteps
        real(kind=4), dimension(nslabs-1), intent(in)                     :: pressure_interslabs
        real(kind=4), dimension(nx,ny,nz,n_timesteps), intent(in)         :: Uhat, Vhat, Q, P 
        ! real(kind=4), INTENT(IN)                                          :: missingval
        real(kind=4), dimension(nx,ny,nslabs,n_timesteps), intent(out)    :: PW_slabs, U_slabs, V_slabs
        real(kind=4), dimension(nx,ny,nslabs-1,n_timesteps), intent(out)  :: Q_interslabs
        integer                                                           :: itime
        
        do itime = 1,n_timesteps
            ! call calculate_slabs0(Uhat(:,:,:,itime), Vhat(:,:,:,itime), Q(:,:,:,itime), P(:,:,:,itime), pressure_interslabs,&
            !                      nx, ny, nz, nslabs, missingval, PW_slabs(:,:,:,itime), U_slabs(:,:,:,itime), V_slabs(:,:,:,itime))
            call calculate_slabs0(Uhat(:,:,:,itime), Vhat(:,:,:,itime), Q(:,:,:,itime), P(:,:,:,itime), pressure_interslabs,&
                                 nx, ny, nz, nslabs, PW_slabs(:,:,:,itime), U_slabs(:,:,:,itime), V_slabs(:,:,:,itime),&
                                 Q_interslabs(:,:,:,itime))
        end do
    
    end subroutine
    
    
    ! subroutine calculate_slabs0(Uhat, Vhat, Q, P, pressure_interslabs, nx, ny, nz, nslabs, missingval, PW_slabs, U_slabs, V_slabs)
    subroutine calculate_slabs0(Uhat, Vhat, Q, P, pressure_interslabs, nx, ny, nz, nslabs, PW_slabs, U_slabs, V_slabs, Q_interslabs)
        use, intrinsic :: ieee_arithmetic
        ! use iso_fortran_env, only: real32, real64
    
        implicit none
        integer, intent(in)                                   :: nx, ny, nz, nslabs
        real(kind=4), dimension(nslabs-1), intent(in)         :: pressure_interslabs
        real(kind=4), dimension(nx,ny,nz), intent(in)         :: Uhat, Vhat, Q, P 
        ! real(kind=4), INTENT(IN)                              :: missingval
        real(kind=4), dimension(nx,ny,nslabs), intent(out)    :: PW_slabs, U_slabs, V_slabs
    
        real(kind=4), dimension(nx,ny,nslabs-1)               :: Uhat_interslabs, Vhat_interslabs
        real(kind=4), dimension(nx,ny,nslabs-1), intent(out)  :: Q_interslabs
        integer, dimension(nx,ny)                             :: i_slab
        real(kind=4)                                          :: p1,p2,q1,q2,uhat1,uhat2,vhat1,vhat2, p_inter, &
                                                                 q_inter, uhat_inter, vhat_inter, pw_layer, u_layer, v_layer
        logical                                               :: is_in_1slab   
        integer                                               :: i,j,k, i_slab0
        real(kind=4)                                          :: g                                                    
    
    
        call DINTERP3DZ(Uhat, Uhat_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, ieee_value(uhat1, ieee_quiet_nan))  ! output is U_interslabs
        call DINTERP3DZ(Vhat, Vhat_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, ieee_value(vhat1, ieee_quiet_nan))  ! output is V_interslabs
        call DINTERP3DZ(Q, Q_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, ieee_value(q1, ieee_quiet_nan))  ! output is Q_interslabs
        ! call DINTERP3DZ(Uhat, Uhat_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, missingval)  ! output is U_interslabs
        ! call DINTERP3DZ(Vhat, Vhat_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, missingval)  ! output is V_interslabs
        ! call DINTERP3DZ(Q, Q_interslabs, P, pressure_interslabs, nx, ny, nz, nslabs-1, missingval)  ! output is Q_interslabs
    
        ! use wrf_constants, only: g
        g = 9.81
        PW_slabs = 0
        U_slabs = 0
        V_slabs = 0
        i_slab = 1
    
        do i_slab0 = 1, nslabs-1
            p_inter = pressure_interslabs(i_slab0)
            do j = 1,ny
                do i = 1,nx
                    if (p_inter > P(i,j,1)) then
                        ! PW_slabs(i,j,i_slab0) = missingval !ieee_value(PW_slabs(i,j,i_slab0), ieee_quiet_nan)
                        ! U_slabs(i,j,i_slab0) = missingval !ieee_value(U_slabs(i,j,i_slab0), ieee_quiet_nan)
                        ! V_slabs(i,j,i_slab0) = missingval !ieee_value(V_slabs(i,j,i_slab0), ieee_quiet_nan)
                        PW_slabs(i,j,i_slab0) = ieee_value(pw_layer, ieee_quiet_nan)
                        U_slabs(i,j,i_slab0) = ieee_value(U_layer, ieee_quiet_nan)
                        V_slabs(i,j,i_slab0) = ieee_value(V_layer, ieee_quiet_nan)
                        i_slab(i,j) = i_slab0 + 1
                    end if
                end do
            end do
        end do
    
        do k = 1,nz-1
            do j = 1,ny
                do i = 1,nx
                    p1 = P(i,j,k)
                    p2 = P(i,j,k+1)
                    
                    q1 = Q(i,j,k)
                    q2 = Q(i,j,k+1)
    
                    uhat1 = Uhat(i,j,k)
                    uhat2 = Uhat(i,j,k+1)
                    vhat1 = Vhat(i,j,k)
                    vhat2 = Vhat(i,j,k+1)
    
                    is_in_1slab = .true.
                    if (i_slab(i,j) < nslabs) then
                        p_inter = pressure_interslabs(i_slab(i,j))
                        if (p_inter > p2) then
                            is_in_1slab = .false.
                        end if
                    end if
    
                    if (is_in_1slab) then
                        pw_layer = (q1 + q2)/2 * (p1 - p2) / g   ! q  * dp / g
                        u_layer =  (uhat1*q1 + uhat2*q2)/2 * (p1 - p2) / g 
                        v_layer =  (vhat1*q1 + vhat2*q2)/2 * (p1 - p2) / g 
                        PW_slabs(i,j,i_slab(i,j)) = PW_slabs(i,j,i_slab(i,j)) + pw_layer 
                        U_slabs(i,j,i_slab(i,j)) = U_slabs(i,j,i_slab(i,j)) + u_layer 
                        V_slabs(i,j,i_slab(i,j)) = V_slabs(i,j,i_slab(i,j)) + v_layer 
                    else           
                        q_inter = Q_interslabs(i,j,i_slab(i,j))
                        uhat_inter = Uhat_interslabs(i,j,i_slab(i,j))
                        vhat_inter = Vhat_interslabs(i,j,i_slab(i,j))
                        ! first layer 
                        pw_layer = (q1 + q_inter)/2 * (p1 - p_inter) / g   ! q  * dp / g   (layer 1)
                        u_layer =  (uhat1*q1 + uhat_inter*q_inter)/2 * (p1 - p_inter) / g 
                        v_layer =  (vhat1*q1 + vhat_inter*q_inter)/2 * (p1 - p_inter) / g 
                        PW_slabs(i,j,i_slab(i,j)) = PW_slabs(i,j,i_slab(i,j)) + pw_layer 
                        U_slabs(i,j,i_slab(i,j)) = U_slabs(i,j,i_slab(i,j)) + u_layer 
                        V_slabs(i,j,i_slab(i,j)) = V_slabs(i,j,i_slab(i,j)) + v_layer 
    
                        ! second layer (in new slab)     ! this sequence is repeated 3 times, probalby make a subroutine
                        i_slab(i,j) = i_slab(i,j) + 1
                        pw_layer = (q_inter + q2)/2 * (p_inter - p2) / g   ! q  * dp / g (layer 2)
                        u_layer =  (uhat_inter*q_inter + uhat2*q2)/2 * (p_inter - p2) / g 
                        v_layer =  (vhat_inter*q_inter + vhat2*q2)/2 * (p_inter - p2) / g 
                        PW_slabs(i,j,i_slab(i,j)) = PW_slabs(i,j,i_slab(i,j)) + pw_layer 
                        U_slabs(i,j,i_slab(i,j)) = U_slabs(i,j,i_slab(i,j)) + u_layer 
                        V_slabs(i,j,i_slab(i,j)) = V_slabs(i,j,i_slab(i,j)) + v_layer 
                    end if
    
                end do
            end do
        end do
    
        do i_slab0 = 1,nslabs
            do j = 1,ny
                do i = 1,nx
                    U_slabs(i,j,i_slab0) = U_slabs(i,j,i_slab0) / PW_slabs(i,j,i_slab0)
                    V_slabs(i,j,i_slab0) = V_slabs(i,j,i_slab0) / PW_slabs(i,j,i_slab0)
                end do
            end do
        end do
    
    end subroutine
    
    ! code taken from wrf python package
    ! NCLFORTSTART
    SUBROUTINE DINTERP3DZ(data3d, out2d, zdata, levels, nx, ny, nz, nlev, missingval)
        IMPLICIT NONE
        ! use iso_fortran_env, only: real32, real64
    
        !f2py threadsafe
        !f2py intent(in,out) :: out2d
    
        INTEGER, INTENT(IN) :: nx, ny, nz, nlev
        real(kind=4), DIMENSION(nx,ny,nz), INTENT(IN) ::  data3d
        real(kind=4), DIMENSION(nx,ny,nlev), INTENT(OUT) :: out2d
        real(kind=4), DIMENSION(nx,ny,nz), INTENT(IN) :: zdata
        real(kind=4), DIMENSION(nlev), INTENT(IN) :: levels
        real(kind=4), INTENT(IN) :: missingval
    
    ! NCLEND
    
        INTEGER :: i,j,kp,ip,im,lev
        LOGICAL :: dointerp
        real(kind=4) :: w1,w2,desiredloc
    
        ! does vertical coordinate increase or decrease with increasing k?
        ! set offset appropriately
    
        ip = 0
        im = 1
        IF (zdata(1,1,1) .GT. zdata(1,1,nz)) THEN
            ip = 1
            im = 0
        END IF
    
        !!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,lev,kp,dointerp,w1,w2,desiredloc) &
        !!$OMP FIRSTPRIVATE(ip,im) SCHEDULE(runtime)
        DO lev = 1,nlev
            DO i = 1,nx
                DO j = 1,ny
                    ! Initialize to missing.  Was initially hard-coded to -999999.
                    out2d(i,j,lev) = missingval
                    dointerp = .FALSE.
                    kp = nz
                    desiredloc = levels(lev)
    
                    DO WHILE ((.NOT. dointerp) .AND. (kp >= 2))
                        ! CHANGED zdata(i,j,kp-ip) >= desiredloc   make pull request
                        ! (zdata(i,j,kp-ip) >= desiredloc) insted of (zdata(i,j,kp-ip) > desiredloc)
                        ! if not does not interplate when there is a pressure level == desired_pressure
                        IF (((zdata(i,j,kp-im) < desiredloc) .AND. (zdata(i,j,kp-ip) >= desiredloc))) THEN
                            w2 = (desiredloc - zdata(i,j,kp-im))/(zdata(i,j,kp-ip) - zdata(i,j,kp-im))
                            w1 = 1.0 - w2 !1.D0 - w2
                            out2d(i,j,lev) = w1*data3d(i,j,kp-im) + w2*data3d(i,j,kp-ip)
                            dointerp = .TRUE.
                        END IF
                        kp = kp - 1
                    END DO
                END DO
            END DO
        END DO
        !!$OMP END PARALLEL DO
    
        RETURN
    
    END SUBROUTINE DINTERP3DZ


    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str
end module mod_transforming