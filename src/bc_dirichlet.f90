module bc_dirichlet_mod
  use udf_type_mod
  use parameters_constant_mod
  use print_msg_mod
  implicit none
  character(18) :: filename(5)


  private :: map_bc_1d_uprofile     
  public  :: initialise_fbcx_given_profile 

  private :: initialise_fbcx_given_const
  private :: initialise_fbcy_given_const
  private :: initialise_fbcz_given_const
  public  :: initialise_fbc_flow_given   ! applied once only, for bc of constant velocity
  public  :: initialise_fbc_thermo_given ! applied once only, for bc of constant temperature
  public  :: enforce_velo_from_fbc

  public :: extract_dirichlet_fbcx
  public :: extract_dirichlet_fbcy
  public :: extract_dirichlet_fbcz
  
contains
  !==========================================================================================================
  !==========================================================================================================
  subroutine extract_dirichlet_fbcx(fbc, var, dtmp)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), intent(out) :: fbc(4,           dtmp%xsz(2), dtmp%xsz(3))
    real(WP), intent(in)  :: var(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3))
    real(WP) :: tmp1, tmp2
    integer  :: j, k
    integer  :: nx, ny, nz

    if(dtmp%xsz(1) /= dtmp%xen(1)) call Print_error_msg("Error. This is not x-pencil.")
    nx = dtmp%xsz(1); ny = dtmp%xsz(2); nz = dtmp%xsz(3)
    !$acc parallel loop collapse(2) default(present) private(tmp1, tmp2)
    do k=1, nz; do j=1, ny
      tmp1 = var(1,  j, k)
      tmp2 = var(nx, j, k)
      fbc(1, j, k) = tmp1
      fbc(2, j, k) = tmp2
      fbc(3, j, k) = tmp1! not used.
      fbc(4, j, k) = tmp2! not used.
    end do; end do
    !$acc end parallel loop

    return
  end subroutine 
  !==========================================================================================================
  subroutine extract_dirichlet_fbcy(fbc, var, dtmp, dm, is_reversed)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    type(t_domain), intent(in) :: dm
    real(WP), intent(out) :: fbc(dtmp%ysz(1), 4,           dtmp%ysz(3))
    real(WP), intent(in)  :: var(dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3))
    logical, optional, intent(in) :: is_reversed

    real(WP), dimension( dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3) ) :: var_zpencil, var_zpencil1 
    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) ) :: var_ypencil1

    integer :: i, j, k
    integer :: nx, ny, nz
    real(WP) :: sign
    !------------------------------------------------------------------------------------------------------
    ! Check if the input data is in y-pencil format
    !------------------------------------------------------------------------------------------------------
    if(dtmp%ysz(2) /= dtmp%yen(2)) call Print_error_msg("Error. This is not y-pencil.")
    !------------------------------------------------------------------------------------------------------
    ! Extract Dirichlet boundary conditions in the y-direction
    !------------------------------------------------------------------------------------------------------
    !$acc data create(var_zpencil, var_zpencil1, var_ypencil1)
    nx = dtmp%ysz(1); ny = dtmp%ysz(2); nz = dtmp%ysz(3)
    !$acc parallel loop collapse(2) default(present)
    do k=1, nz; do i=1, nx
      fbc(i, 1, k) = var(i, 1,  k) ! Lower boundary
      fbc(i, 2, k) = var(i, ny, k) ! Upper boundary
      fbc(i, 4, k) = TWO * var(i, ny, k) - var(i, ny-1, k)! Upper boundary
    end do; end do
    !$acc end parallel loop
    !------------------------------------------------------------------------------------------------------
    ! Handle special treatment of the lower boundary for pipe geometry (ICASE_PIPE)
    ! this part is the same as axis_mirroring
    !------------------------------------------------------------------------------------------------------
    if(dm%icase == ICASE_PIPE) then
      sign = ONE
      if(present(is_reversed)) then
        if(is_reversed) sign = -ONE
      end if
      call transpose_y_to_z(var, var_zpencil, dtmp)
      nx = dtmp%zsz(1); ny = dtmp%zsz(2); nz = dtmp%zsz(3)
      !$acc parallel loop collapse(3) default(present)
      do k = 1, nz; do j = 1, ny; do i = 1, nx
        var_zpencil1(i, j, k) = sign * var_zpencil(i, j, dm%knc_sym(k))
      end do; end do; end do
      !$acc end parallel loop
      call transpose_z_to_y(var_zpencil1, var_ypencil1, dtmp)
      !$acc kernels default(present)
      fbc(:, 3, :) = var_ypencil1(:, 2, :)
      !$acc end kernels
    else
      !$acc kernels default(present)
      fbc(:, 3, :) = TWO * var(:, 1, :) - var(:, 2, :)
      !$acc end kernels
    end if
    !$acc end data

    return
  end subroutine 
  !==========================================================================================================
  subroutine extract_dirichlet_fbcz(fbc, var, dtmp)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), intent(in)  :: var(dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3))
    real(WP), intent(out) :: fbc(dtmp%zsz(1), dtmp%zsz(2), 4          )
    real(WP) :: tmp1, tmp2
    integer  :: i ,j
    integer  :: nx, ny, nz

    if(dtmp%zsz(3) /= dtmp%zen(3)) call Print_error_msg("Error. This is not z-pencil.")
    nx = dtmp%zsz(1); ny = dtmp%zsz(2); nz = dtmp%zsz(3)
    !$acc parallel loop collapse(2) default(present) private(tmp1, tmp2)
    do j=1, ny; do i=1, nx
      tmp1 = var(i, j, 1)
      tmp2 = var(i, j, nz)
      fbc(i, j, 1) = tmp1
      fbc(i, j, 2) = tmp2
      fbc(i, j, 3) = tmp1
      fbc(i, j, 4) = tmp2
    end do; end do
    !$acc end parallel loop

    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine  map_bc_1d_uprofile(filename, n, y, u)
    use io_files_mod
    implicit none
    character(*), intent(in) :: filename
    integer,  intent(in)  :: n
    real(WP), intent(in)  :: y(n)
    real(WP), intent(out) :: u(n)

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    character(len = 80) :: str
    
    real(WP) :: rtmp
    integer :: i, nn
    integer :: pf_unit, j
    real(WP), allocatable :: uprofile(:)
    real(WP), allocatable :: yy(:)

    !----------------------------------------------------------------------------------------------------------
    ! to read given u-velocity profile, dimensionless, x/H, U/Umean
    !----------------------------------------------------------------------------------------------------------
    open ( newunit = inputUnit,     &
          file    = trim(filename), &
          status  = 'old',         &
          action  = 'read',        &
          iostat  = ioerr,         &
          iomsg   = iotxt)
    if(ioerr /= 0) then
      str = 'Problem opening: '//trim(filename)
      call Print_error_msg(trim(str))
    end if

    nn = 0
    read(inputUnit, *, iostat = ioerr) str

    do
      read(inputUnit, *, iostat = ioerr) rtmp, rtmp
      if(ioerr /= 0) exit
      nn = nn + 1
    end do
    rewind(inputUnit)
    !----------------------------------------------------------------------------------------------------------
    ! to read given u-velocity profile, dimensionless, x/H, U/Umean
    !----------------------------------------------------------------------------------------------------------
    allocate ( uprofile (nn) )
    allocate ( yy (nn) )

    read(inputUnit, *, iostat = ioerr) str
    do i = 1, nn
      read(inputUnit, *, iostat = ioerr) yy(i), uprofile(i)
    end do
    close(inputUnit)

    call profile_interpolation(nn, yy, uprofile, n, y, u)

    if(nrank == 0) then
      open ( newunit = pf_unit,     &
              file    = trim(dir_chkp)//'/check_given_ux_profile.dat', &
              status  = 'replace',         &
              action  = 'write')
      write(pf_unit, '(A)') "#j, y, u - original"
      do j = 1, nn
        write(pf_unit, '(1I3.1, 5ES15.7)') j, yy(j), uprofile(j)
      end do
      write(pf_unit, '(A)') "#j, y, u - interpolation"
      do j = 1, n
        write(pf_unit, '(1I3.1, 5ES15.7)') j, y(j), u(j)
      end do
      close(pf_unit)
    end if

    deallocate(uprofile)
    deallocate(yy)

    return
  end subroutine

  !==========================================================================================================
  !==========================================================================================================
  subroutine initialise_fbcx_given_profile(fbcx, var1y, jst, str)
    use io_files_mod
    implicit none
    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: var1y(:)
    integer,  intent(in)    :: jst
    character(2), intent(in)   :: str

    integer :: k, j, jj
    integer :: ny, nz
    integer :: pf_unit

    ny = size(fbcx, 2)
    nz = size(fbcx, 3)
    !$acc parallel loop collapse(2) private(jj) default(present)
    do k = 1, nz
      do j = 1, ny
        jj = jst + j - 1
        fbcx(1, j, k) = var1y(jj)
        fbcx(3, j, k) = fbcx(1, j, k)
      end do
    end do
    !$acc end parallel loop

    !$acc update self(fbcx)
    open ( newunit = pf_unit,     &
            file    = trim(dir_chkp)//'/check_given_'//trim(str)//'_profile.dat', &
            position= 'append',         &
            action  = 'write')
    write(pf_unit, '(A)') "#fbcx"
    do j = 1, size(fbcx, 2)
      write(pf_unit, '(1I3.1, 1ES15.7)') j, fbcx(1, j, 1)
    end do
    close(pf_unit)

    return
  end subroutine
  !==========================================================================================================
  subroutine initialise_fbcx_given_const(fbcx, fbcx_const)
    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: fbcx_const(2)

    integer :: k, j, n
    integer :: ny, nz

    ny = size(fbcx, 2)
    nz = size(fbcx, 3)
    !$acc parallel loop collapse(3) present(fbcx)
    do k = 1, nz
      do j = 1, ny
        do n = 1, 2
          fbcx(n,   j, k) = fbcx_const(n)
          fbcx(n+2, j, k) = fbcx(n, j, k)
        end do
      end do
    end do
    !$acc end parallel loop

    return
  end subroutine
  !==========================================================================================================
  subroutine initialise_fbcy_given_const(fbcy, fbcy_const, ri)
    real(WP), intent(inout) :: fbcy(:, :, :)
    real(WP), intent(in)    :: fbcy_const(2)
    real(WP), intent(in), optional :: ri(:)

    integer :: k, i, n
    integer :: nx, nz
    real(WP) :: ri_new(2)

    if(present(ri)) then
      ri_new(1) = ri(1)
      ri_new(2) = ri(size(ri))
    else
      ri_new = ONE
    end if

    nx = size(fbcy, 1)
    nz = size(fbcy, 3)
    !$acc parallel loop collapse(3) present(fbcy)
    do k = 1, nz
      do i = 1, nx
        do n = 1, 2
          if(ri_new(n) < (MAXP * HALF)) then
            fbcy(i, n,   k) =  fbcy_const(n) * ri_new(n)
          end if
          fbcy(i, n+2, k) =  fbcy(i, n, k)
        end do
      end do
    end do
    !$acc end parallel loop

    return
  end subroutine
  !==========================================================================================================
  subroutine initialise_fbcz_given_const(fbcz, fbcz_const, ri, jst)
    real(WP), intent(inout) :: fbcz(:, :, :)
    real(WP), intent(in)    :: fbcz_const(2)
    integer,  intent(in), optional :: jst
    real(WP), intent(in), optional :: ri(:)

    integer :: i, j, n, jj
    integer :: nx, ny

    jj = 0
    nx = size(fbcz, 1)
    ny = size(fbcz, 2)
    !$acc parallel loop collapse(3) private(jj) present(fbcz)
    do j = 1, ny
      do i = 1, nx
        do n = 1, 2
          if(present(jst) .and. present(ri)) then
            jj = jst + j - 1
            fbcz(i, j, n  ) =  fbcz_const(n) / ri(jj)
          else
            fbcz(i, j, n  ) =  fbcz_const(n)
          end if
          fbcz(i, j, n+2) =  fbcz(i, j, n)
        end do
      end do
    end do
    !$acc end parallel loop

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine initialise_fbc_flow_given (dm) ! apply once only
    type(t_domain), intent(inout)   :: dm

    real(WP) :: var1y(1:dm%np(2))
    
    integer :: nx, ny, nz, n
    integer :: i, j, k
!==========================================================================================================
! to build up bc with constant values
! -3-1-||||-2-4
! for constant bc, 3=1= geometric bc, side 1;
!                  2=4= geometric bc, side 2 
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! x-bc in x-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call initialise_fbcx_given_const(dm%fbcx_qx, dm%fbcx_const(:, 1))
    call initialise_fbcx_given_const(dm%fbcx_qy, dm%fbcx_const(:, 2))
    call initialise_fbcx_given_const(dm%fbcx_qz, dm%fbcx_const(:, 3))
    call initialise_fbcx_given_const(dm%fbcx_pr, dm%fbcx_const(:, 4))
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call initialise_fbcy_given_const(dm%fbcy_qx, dm%fbcy_const(:, 1))
    call initialise_fbcy_given_const(dm%fbcy_qy, dm%fbcy_const(:, 2))
    call initialise_fbcy_given_const(dm%fbcy_qz, dm%fbcy_const(:, 3)) ! geo_bc, rpi, not rci
    call initialise_fbcy_given_const(dm%fbcy_pr, dm%fbcy_const(:, 4))
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qyr = qy/r = uy
!----------------------------------------------------------------------------------------------------------
      call initialise_fbcy_given_const(dm%fbcy_qyr, dm%fbcy_const(:, 2), dm%rpi)
      call initialise_fbcy_given_const(dm%fbcy_qzr, dm%fbcy_const(:, 3), dm%rci)
    end if
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call initialise_fbcz_given_const(dm%fbcz_qx, dm%fbcz_const(:, 1))
    call initialise_fbcz_given_const(dm%fbcz_qy, dm%fbcz_const(:, 2))
    call initialise_fbcz_given_const(dm%fbcz_qz, dm%fbcz_const(:, 3))
    call initialise_fbcz_given_const(dm%fbcz_pr, dm%fbcz_const(:, 4))
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qyr = qy/r = uy
!----------------------------------------------------------------------------------------------------------
      call initialise_fbcz_given_const(dm%fbcz_qyr, dm%fbcz_const(:, 2), dm%rpi, dm%dcpc%zst(2))
      call initialise_fbcz_given_const(dm%fbcz_qzr, dm%fbcz_const(:, 3), dm%rci, dm%dccp%zst(2))
    end if

!==========================================================================================================
! to build up bc for var(x_const, y, z)
!==========================================================================================================
    !$acc data create(var1y)
!----------------------------------------------------------------------------------------------------------
! x-bc1, qx(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 1) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(1) = trim('PF1D_U1Y.DAT') !(undim)
      ! call map_bc_1d_uprofile( filename(1), ny, dm%yc, var1y(1:ny) )
      ! call initialise_fbcx_given_profile(dm%fbcx_qx, var1y, dm%dpcc%xst(2), 'qx')
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qy(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 2) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%np(2)
      filename(2) = trim('PF1D_V1Y.DAT') !(undim)
      call map_bc_1d_uprofile( filename(2), ny, dm%yp, var1y(1:ny) )
      if(dm%icoordinate == ICARTESIAN) var1y(1:ny) =  var1y(1:ny) * dm%rp(1:ny)
      !$acc update device(var1y)
      call initialise_fbcx_given_profile(dm%fbcx_qy, var1y, dm%dcpc%xst(2), 'qy')
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qz(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 3) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(3) = trim('PF1D_W1Y.DAT') !(undim)
      call map_bc_1d_uprofile( filename(3), ny, dm%yc, var1y(1:ny) )
      if(dm%icoordinate == ICARTESIAN) var1y(1:ny) =  var1y(1:ny) * dm%rc(1:ny)
      !$acc update device(var1y)
      call initialise_fbcx_given_profile(dm%fbcx_qz, var1y, dm%dccp%xst(2), 'qz')
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, pr(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 4) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(4) = trim('PF1D_P1Y.DAT') !(undim)
      call map_bc_1d_uprofile( filename(4), ny, dm%yc, var1y(1:ny) )
      !$acc update device(var1y)
      call initialise_fbcx_given_profile(dm%fbcx_pr, var1y, dm%dccc%xst(2), 'pr')
    end if
    !$acc end data

    if(dm%is_thermo) then
      do n = 1, 2
!----------------------------------------------------------------------------------------------------------
! x-bc in x-pencil, gx, gy, gz
!----------------------------------------------------------------------------------------------------------
        ny = dm%dccc%xsz(2)
        nz = dm%dccc%xsz(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1, nz; do j=1, ny
          dm%fbcx_gx(n, j, k) = dm%fbcx_qx(n, j, k) * dm%fbcx_ftp(n, 1, 1)%d
          dm%fbcx_gy(n, j, k) = dm%fbcx_qy(n, j, k) * dm%fbcx_ftp(n, 1, 1)%d
          dm%fbcx_gz(n, j, k) = dm%fbcx_qz(n, j, k) * dm%fbcx_ftp(n, 1, 1)%d
        end do; end do
        !$acc end parallel loop
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, gx, gy, gz, qyr = qy/r = uy
!----------------------------------------------------------------------------------------------------------
        nx = dm%dccc%ysz(1)
        nz = dm%dccc%ysz(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1, nz; do i=1, nx
          dm%fbcy_gx(i, n, k) = dm%fbcy_qx(i, n, k) * dm%fbcy_ftp(1, n, 1)%d
          dm%fbcy_gy(i, n, k) = dm%fbcy_qy(i, n, k) * dm%fbcy_ftp(1, n, 1)%d
          dm%fbcy_gz(i, n, k) = dm%fbcy_qz(i, n, k) * dm%fbcy_ftp(1, n, 1)%d
        end do; end do
        !$acc end parallel loop
        !if(dm%icoordinate == ICYLINDRICAL) then
          !dm%fbcy_gyr(:, n, :) = dm%fbcy_qyr(:, n, :) * dm%fbcy_ftp(1, n, 1)%d
          !dm%fbcy_gzr(:, n, :) = dm%fbcy_qzr(:, n, :) * dm%fbcy_ftp(1, n, 1)%d
        !end if
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, gx, gy, gz, gyr = gy/r
!----------------------------------------------------------------------------------------------------------
        nx = dm%dccc%zsz(1)
        ny = dm%dccc%zsz(2)
        !$acc parallel loop collapse(2) default(present)
        do j=1, ny; do i=1, nx
          dm%fbcz_gx(i, j, n) = dm%fbcz_qx(i, j, n) * dm%fbcz_ftp(1, 1, n)%d
          dm%fbcz_gy(i, j, n) = dm%fbcz_qy(i, j, n) * dm%fbcz_ftp(1, 1, n)%d
          dm%fbcz_gz(i, j, n) = dm%fbcz_qz(i, j, n) * dm%fbcz_ftp(1, 1, n)%d
        end do; end do
        !$acc end parallel loop
        !if(dm%icoordinate == ICYLINDRICAL) then
          !dm%fbcz_gyr(:, :, n) = dm%fbcz_qyr(:, :, n) * dm%fbcz_ftp(1, 1, n)%d
          !dm%fbcz_gzr(:, :, n) = dm%fbcz_qzr(:, :, n) * dm%fbcz_ftp(1, 1, n)%d
        !end if

      end do
    end if

    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine initialise_fbc_thermo_given(tm, dm) ! call this after scaling the fbc_ftp values
    use thermo_info_mod
    use decomp_2d
    type(t_domain), intent(inout) :: dm
    type(t_thermo), intent(in)    :: tm

    real(WP) :: var1y(1:dm%np(2))
    real(WP), allocatable :: ac4c_ypencil(:, :, :), ac4c_xpencil(:, :, :)
    integer :: nxbf, nx, ny, nz, n
    integer :: i, j, k
    type(DECOMP_INFO) :: dtmp
!----------------------------------------------------------------------------------------------------------
! to build up bc with constant values
! -3-1-||||-2-4
! for constant bc, 3=1= geometric bc, side 1;
!                  2=4= geometric bc, side 2 
!----------------------------------------------------------------------------------------------------------
! x-bc1, pr(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    !$acc data create(var1y)
    if(dm%ibcx_nominal(1, 5) == IBC_PROFILE1D) then !
      var1y = ZERO
      ny = dm%nc(2)
      filename(4) = trim('pf1d_T1y_undim.dat') !(undim)
      call map_bc_1d_uprofile( filename(4), ny, dm%yc, var1y(1:ny) )
      !$acc update device(var1y)
      call initialise_fbcx_given_profile(dm%fbcx_ftp(:,:,:)%t, var1y, dm%dccc%xst(2),'Ty')
      call ftp_refresh_thermal_properties_from_T_undim_3Dftp(dm%fbcx_ftp)
    else 
      if(nrank == 0 .and. dm%ibcx_nominal(1, 1) == IBC_DATABASE) &
      call Print_warning_msg("The thermal field's inlet temperature is the same as Tini given.")
    end if
    !$acc end data

    nxbf = 0
    if((dm%inlet_tbuffer_len - dm%h(1)) > MINP) then
      if(dm%inlet_tbuffer_len > dm%lxx) then 
        call Print_warning_msg("The inlet thermal buffer layer exceeds the domain length and has been reduced to 1/10 of the domain length.")
        dm%inlet_tbuffer_len = dm%lxx / TEN
      end if
      call decomp_info_init(dm%nc(1), 4,  dm%nc(3), dtmp) 
      nxbf = floor(dm%inlet_tbuffer_len * dm%h1r(1))
    end if

    do n = 1, 2
!----------------------------------------------------------------------------------------------------------
! x-bc in x-pencil, ftp, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
      if(dm%ibcx_nominal(n, 5) == IBC_DIRICHLET) then
        !$acc kernels default(present)
        dm%fbcx_ftp(n, :, :)%t   = dm%fbcx_const(n, 5)
        !$acc end kernels
        call ftp_refresh_thermal_properties_from_T_undim_3Dftp(dm%fbcx_ftp(n:n, :, :))
      else if (dm%ibcx_nominal(n, 5) == IBC_NEUMANN) then
        !NOTE: this currently works, but may be subtle on GPU
        !      if necessary, switch to the member-by-member operations
        !      similar issues have been observed in Update_thermal_properties
        !$acc kernels default(present)
        dm%fbcx_qw(n, :, :)  = dm%fbcx_const(n, 5)
        dm%fbcx_ftp(n, :, :) = tm%ftp_ini
        !$acc end kernels

        ! TODO: performance comparison with the above may needed
!        ny = dm%dpcc%xsz(2); nz = dm%dpcc%xsz(3)
!        !$acc parallel loop collapse(2) default(present)
!        do k = 1, nz; do j = 1, ny
!          dm%fbcx_qw(n, j, k) = dm%fbcx_const(n, 5)
!          dm%fbcx_ftp(n, j, k)%t = tm%ftp_ini%t
!          dm%fbcx_ftp(n, j, k)%d = tm%ftp_ini%d
!          dm%fbcx_ftp(n, j, k)%m = tm%ftp_ini%m
!          dm%fbcx_ftp(n, j, k)%k = tm%ftp_ini%k
!          dm%fbcx_ftp(n, j, k)%h = tm%ftp_ini%h
!          dm%fbcx_ftp(n, j, k)%rhoh = tm%ftp_ini%rhoh
!          dm%fbcx_ftp(n, j, k)%cp = tm%ftp_ini%cp
!          dm%fbcx_ftp(n, j, k)%b = tm%ftp_ini%b
!          dm%fbcx_ftp(n, j, k)%alpha = tm%ftp_ini%alpha
!          dm%fbcx_ftp(n, j, k)%Pr = tm%ftp_ini%Pr
!        end do; end do
!        !$acc end parallel loop
      else
        !$acc kernels default(present)
        dm%fbcx_ftp(n, :, :) = tm%ftp_ini
        !$acc end kernels

!        ny = dm%dpcc%xsz(2); nz = dm%dpcc%xsz(3)
!        !$acc parallel loop collapse(2) default(present)
!        do k = 1, nz; do j = 1, ny
!          dm%fbcx_ftp(n, j, k)%t = tm%ftp_ini%t
!          dm%fbcx_ftp(n, j, k)%d = tm%ftp_ini%d
!          dm%fbcx_ftp(n, j, k)%m = tm%ftp_ini%m
!          dm%fbcx_ftp(n, j, k)%k = tm%ftp_ini%k
!          dm%fbcx_ftp(n, j, k)%h = tm%ftp_ini%h
!          dm%fbcx_ftp(n, j, k)%rhoh = tm%ftp_ini%rhoh
!          dm%fbcx_ftp(n, j, k)%cp = tm%ftp_ini%cp
!          dm%fbcx_ftp(n, j, k)%b = tm%ftp_ini%b
!          dm%fbcx_ftp(n, j, k)%alpha = tm%ftp_ini%alpha
!          dm%fbcx_ftp(n, j, k)%Pr = tm%ftp_ini%Pr
!        end do; end do
!        !$acc end parallel loop
      end if
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, ftp, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
      if(dm%ibcy_nominal(n, 5) == IBC_DIRICHLET) then
        !$acc kernels default(present)
        dm%fbcy_ftp(:, n, :)%t   = dm%fbcy_const(n, 5)
        !$acc end kernels
        call ftp_refresh_thermal_properties_from_T_undim_3Dftp(dm%fbcy_ftp(:, n:n, :))
      else if (dm%ibcy_nominal(n, 5) == IBC_NEUMANN) then
        !$acc kernels default(present)
        dm%fbcy_qw(:, n, :) = dm%fbcy_const(n, 5)
        dm%fbcy_ftp(:, n, :) = tm%ftp_ini
        !$acc end kernels
      else
        !$acc kernels default(present)
        dm%fbcy_ftp(:, n, :) = tm%ftp_ini
        !$acc end kernels
      end if

      ! a patch for inlet buffer layer
      if((dm%inlet_tbuffer_len - dm%h(1)) > MINP) then
        allocate ( ac4c_xpencil(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)) )
        allocate ( ac4c_ypencil(dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3)) )
        !$acc data create(ac4c_xpencil, ac4c_ypencil)
        !$acc kernels default(present)
        ac4c_ypencil(:,:,:) = dm%fbcy_const(n, 5)
        !$acc end kernels
        call transpose_y_to_x(ac4c_ypencil, ac4c_xpencil, dtmp)
        if(dm%ibcy_nominal(n, 5) == IBC_DIRICHLET) then
          !$acc kernels default(present)
          ac4c_xpencil(1:nxbf, :, :) = ONE
          !$acc end kernels
        else if (dm%ibcy_nominal(n, 5) == IBC_NEUMANN) then
          !$acc kernels default(present)
          ac4c_xpencil(1:nxbf, :, :) = ZERO
          !$acc end kernels
        else
          !$acc kernels default(present)
          ac4c_xpencil(1:nxbf, :, :) = dm%fbcy_const(n, 5)
          !$acc end kernels
        end if
        call transpose_x_to_y(ac4c_xpencil, ac4c_ypencil, dtmp)
        if(dm%ibcy_nominal(n, 5) == IBC_DIRICHLET) then
          !$acc kernels default(present)
          dm%fbcy_ftp(:, n, :)%t = ac4c_ypencil(:, n, :)
          !$acc end kernels
          call ftp_refresh_thermal_properties_from_T_undim_3Dftp(dm%fbcy_ftp(:, n:n, :))
        else if (dm%ibcy_nominal(n, 5) == IBC_NEUMANN) then
          !$acc kernels default(present)
          dm%fbcy_qw(:, n, :) = ac4c_ypencil(:, n, :)
          !$acc end kernels
        end if
        !$acc end data
        deallocate(ac4c_ypencil, ac4c_xpencil)
      end if
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
      if( dm%ibcz_nominal(n, 5) == IBC_DIRICHLET ) then
        !$acc kernels default(present)
        dm%fbcz_ftp(:, :, n)%t   = dm%fbcz_const(n, 5)
        !$acc end kernels
        call ftp_refresh_thermal_properties_from_T_undim_3Dftp(dm%fbcz_ftp(:, :, n:n))
      else if (dm%ibcz_nominal(n, 5) == IBC_NEUMANN) then
        !$acc kernels default(present)
        dm%fbcz_qw(:, :, n) = dm%fbcz_const(n, 5)
        dm%fbcz_ftp(:, :, n) = tm%ftp_ini
        !$acc end kernels
      else
        !$acc kernels default(present)
        dm%fbcz_ftp(:, :, n) = tm%ftp_ini
        !$acc end kernels
      end if
    end do
    ! TODO: This currently works, but can be arranged into the above loops for better performance
    !       Since only doing once, it can remain unchanged
    !$acc kernels default(present)
    dm%fbcx_ftp(3, :, :) = dm%fbcx_ftp(1, :, :)
    dm%fbcx_ftp(4, :, :) = dm%fbcx_ftp(2, :, :)
    dm%fbcy_ftp(:, 3, :) = dm%fbcy_ftp(:, 1, :)
    dm%fbcy_ftp(:, 4, :) = dm%fbcy_ftp(:, 2, :)
    dm%fbcz_ftp(:, :, 3) = dm%fbcz_ftp(:, :, 1)
    dm%fbcz_ftp(:, :, 4) = dm%fbcz_ftp(:, :, 2)
    !$acc end kernels

    return
  end subroutine
!==========================================================================================================
  subroutine enforce_velo_from_fbc(dm, ux, uy, uz, fbcx0, fbcy0, fbcz0)
    use udf_type_mod
    use parameters_constant_mod
    use print_msg_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (inout) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (inout) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (inout) :: uz
    real(WP), dimension(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)), optional, intent (in) :: fbcx0
    real(WP), dimension(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)), optional, intent (in) :: fbcy0
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), 4), optional, intent (in) :: fbcz0

    real(WP), dimension(dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: acpc_ypencil
    real(WP), dimension(dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3)) :: accp_ypencil
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3)) :: accp_zpencil
    real(WP), dimension(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) :: fbcx
    real(WP), dimension(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)) :: fbcy
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), 4) :: fbcz

    integer :: ibcx_qx1, ibcx_qx2, ibcy_qy1, ibcy_qy2, ibcz_qz1, ibcz_qz2
    integer, dimension(3) ::  npccx, ncpcy, nccpz
    integer :: i, j, k, nx, ny ,nz

    npccx = dm%dpcc%xsz
    ncpcy = dm%dcpc%ysz
    nccpz = dm%dccp%zsz

    !$acc enter data create(acpc_ypencil, accp_ypencil, accp_zpencil, fbcx, fbcy, fbcz)

    !$acc kernels default(present)
    if(.not. present(fbcx0)) then
      fbcx(:,:,:) = ZERO
    else
      fbcx(:,:,:) = fbcx0(:,:,:)
    end if
    !$acc end kernels

    !$acc kernels default(present)
    if(.not. present(fbcy0)) then
      fbcy(:,:,:) = ZERO
    else
      fbcy(:,:,:) = fbcy0(:,:,:)
    end if
    !$acc end kernels

    !$acc kernels default(present)
    if(.not. present(fbcz0)) then
      fbcz(:,:,:) = ZERO
    else
      fbcz(:,:,:) = fbcz0(:,:,:)
    end if
    !$acc end kernels

    ! -mx_rhs-
    ibcx_qx1 = dm%ibcx_qx(1)
    ibcx_qx2 = dm%ibcx_qx(2)
    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(2) default(present)
    do k=1,nz; do j=1,ny
      if(ibcx_qx1 == IBC_DIRICHLET) then
        ux(1, j, k) = fbcx(1, j, k)
        fbcx(3, j, k) = TWO*fbcx(1, j, k) - ux(2, j, k)
      end if
      if(ibcx_qx2 == IBC_DIRICHLET) then
        ux(nx, j, k) = fbcx(2, j, k)
        fbcx(4, j, k) = TWO*fbcx(2, j, k) - ux(nx-1, j, k)
      end if
    end do; end do
    !$acc end parallel loop

    !-my_rhs-
    ibcy_qy1 = dm%ibcy_qy(1)
    ibcy_qy2 = dm%ibcy_qy(2)
    if(ibcy_qy1 == IBC_DIRICHLET .or. &
       ibcy_qy2 == IBC_DIRICHLET) then
      call transpose_x_to_y(uy, acpc_ypencil, dm%dcpc)
      nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
      !$acc parallel loop collapse(2) default(present)
      do k=1,nz; do i=1,nx
        if(ibcy_qy1 == IBC_DIRICHLET) then
          acpc_ypencil(i, 1, k) = fbcy(i, 1, k)
          fbcy(i, 3, k) = TWO*fbcy(i, 1, k) - acpc_ypencil(i, 2, k)
        end if
        if(ibcy_qy2 == IBC_DIRICHLET) then
          acpc_ypencil(i, ny, k) = fbcy(i, 2, k)
          fbcy(i, 4, k) = TWO*fbcy(i, 2, k) - acpc_ypencil(i, ny-1, k)
        end if
      end do; end do
      !$acc end parallel loop
      call transpose_y_to_x(acpc_ypencil, uy, dm%dcpc)
    end if
    if(dm%icase == ICASE_PIPE) then ! due to qr=ur * r, when r=0, qr=0
      call transpose_x_to_y(uy, acpc_ypencil, dm%dcpc)
      !$acc kernels default(present)
      acpc_ypencil(:, 1, :) = ZERO
      !$acc end kernels
      call transpose_y_to_x(acpc_ypencil, uy, dm%dcpc)
    end if

    !-mz_rhs-
    ibcz_qz1 = dm%ibcz_qz(1)
    ibcz_qz2 = dm%ibcz_qz(2)
    if(ibcz_qz1  == IBC_DIRICHLET .or. &
       ibcz_qz2  == IBC_DIRICHLET) then
      call transpose_x_to_y(uz, accp_ypencil, dm%dccp)
      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
      nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
      !$acc parallel loop collapse(2) default(present)
      do j=1,ny; do i=1,nx
        if(ibcz_qz1 == IBC_DIRICHLET) then
          accp_zpencil(i, j, 1) = fbcz(i, j, 1)
          fbcz(i, j, 3) = TWO * fbcz(i, j, 1) - accp_zpencil(i, j, 2)
        end if
        if(ibcz_qz2 == IBC_DIRICHLET) then 
          accp_zpencil(i, j, nz) = fbcz(i, j, 2)
          fbcz(i, j, 4) = TWO * fbcz(i, j, 2) - accp_zpencil(i, j, nz-1)
        end if
      end do; end do
      !$acc end parallel loop
      call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccp)
      call transpose_y_to_x(accp_ypencil, uz, dm%dccp)
    end if

    !$acc exit data delete(acpc_ypencil, accp_ypencil, accp_zpencil, fbcx, fbcy, fbcz)

    return
  end subroutine

end module
