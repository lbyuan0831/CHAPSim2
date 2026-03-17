module eq_energy_mod
  use operations
  use decomp_2d
  use wrt_debug_field_mod
  implicit none

  private :: Compute_energy_rhs
  private :: Calculate_energy_fractional_step
  public  :: Calculate_drhodt
  public  :: Update_thermal_properties
  public  :: Solve_energy_eq
contains
!==========================================================================================================
  subroutine Calculate_drhodt(fl, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    use find_max_min_ave_mod
    use solver_tools_mod
    implicit none
    type(t_domain), intent(in) :: dm
    integer, intent(in) :: isub
    type(t_flow), intent(inout) :: fl
    !real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(in)  :: dens, densm1
    !real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(out) :: drhodt
    !real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div
    !real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div0
    integer, dimension(3) ::  ncccx
    integer :: i, j, k
    integer :: nx, ny, nz
    real(WP) :: maxdrhodt!, d1_bulk, d0_bulk
    logical :: do_projection

    ! Default value
    fl%drhodt = ZERO
    if( .not. dm%is_thermo) return
    if (.not. is_RK_proj(isub)) return
    ! thermal field is half time step ahead of the velocity field
    ! -----*-----$-----*-----$-----*-----$-----*-----$-----*-----$-----*
    !           d_(i-1)     d_i   u_i    d_(i+1)
    ! select case (dm%iTimeScheme)
    !   case (ITIME_EULER) ! 1st order
    !     fl%drhodt = fl%dDens - fl%dDens0
    !     fl%drhodt = fl%drhodt / dm%dt
    !   case (ITIME_AB2) ! 2nd order
    !     fl%drhodt = fl%dDens - fl%dDensm2
    !     fl%drhodt = fl%drhodt / (TWO * dm%dt)
    !   case (ITIME_RK3, ITIME_RK3_CN) ! 3rd order
    !     fl%drhodt = -fl%dDensm2 + SIX * fl%dDens0 - THREE * fl%dDens
    !     fl%drhodt = fl%drhodt / (dm%dt * SIX)
    !   case default
    !     fl%drhodt = fl%dDens - fl%dDens0
    !     fl%drhodt = fl%drhodt / dm%dt
    ! end select
    !------------------------------------------------------------
    ! Compute drho/dt
    !------------------------------------------------------------
    !----------------------------------------------------------
    ! Finite-difference form:
    !   drho/dt = (rho - rho0) / dt
    !----------------------------------------------------------
    ncccx = dm%dccc%xsz
    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%drhodt(i,j,k) = fl%dDens(i,j,k) - fl%dDens0(i,j,k)
      ! Apply time scaling (common to both formulations)
      fl%drhodt(i,j,k) = fl%drhodt(i,j,k) / dm%dt
    end do; end do; end do
    !$acc end parallel loop
    !------------------------------------------------------------
    ! damp drho/dt near in/out b.c.
    !------------------------------------------------------------
    if(is_damping_drhodt) call damping_drhodt(fl%drhodt, dm)
    
    return
  end subroutine Calculate_drhodt
!==========================================================================================================
  subroutine Update_thermal_properties(dens, visc, tm, fl, dm)
    use parameters_constant_mod
    use udf_type_mod
    use operations
    use thermo_info_mod
    use cylindrical_rn_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(inout) :: dens, visc

    real(WP), pointer, dimension(:,:,:) :: dh_pcc, d_pcc
    real(WP), pointer, dimension(:,:,:) :: dh_ypencil, d_ypencil
    real(WP), pointer, dimension(:,:,:) :: dh_cpc_ypencil, d_cpc_ypencil
    real(WP), pointer, dimension(:,:,:) :: dh_zpencil, d_zpencil
    real(WP), pointer, dimension(:,:,:) :: dh_ccp_zpencil, d_ccp_zpencil
    real(WP), pointer, dimension(:,:,:) :: fbcy_c4c
    integer, dimension(3) :: ncccx, ncccy, ncccz, npccx, ncpcy, nccpz
    integer :: i, j, k
    integer :: nx, ny, nz
    type(t_fluidThermoProperty) :: ftp

    ncccx = dm%dccc%xsz
    ncccy = dm%dccc%ysz
    ncccz = dm%dccc%zsz
    npccx = dm%dpcc%xsz
    ncpcy = dm%dcpc%ysz
    nccpz = dm%dccp%zsz
!----------------------------------------------------------------------------------------------------------
!   main field
!----------------------------------------------------------------------------------------------------------
    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) private(ftp) default(present)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ftp%rhoh = tm%rhoh(i, j, k)
          ftp%d = dens(i, j, k)
          call ftp_refresh_thermal_properties_from_DH(ftp)
          tm%hEnth(i, j, k) = ftp%h
          tm%tTemp(i, j, k) = ftp%T
          tm%kCond(i, j, k) = ftp%k
          dens(i, j, k) = ftp%d
          visc(i, j, k) = ftp%m
        end do
      end do
    end do
    !$acc end parallel loop
!----------------------------------------------------------------------------------------------------------
!  BC - x
!----------------------------------------------------------------------------------------------------------
  if( dm%ibcx_Tm(1) == IBC_NEUMANN .or. &
      dm%ibcx_Tm(2) == IBC_NEUMANN) then
    dh_pcc(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    d_pcc (1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
    call Get_x_midp_C2P_3D(tm%rhoh, dh_pcc, dm, dm%iAccuracy, dm%ibcx_ftp) ! exterpolation, check
    call Get_x_midp_C2P_3D(dens,d_pcc, dm, dm%iAccuracy, dm%ibcx_ftp)
    if(dm%ibcx_Tm(1) == IBC_NEUMANN .and. &
       dm%dpcc%xst(1) == 1) then
      ny = size(dm%fbcx_ftp, 2)
      nz = size(dm%fbcx_ftp, 3)
      !$acc parallel loop collapse(2) private(ftp) default(present)
      do k = 1, nz
        do j = 1, ny
          ftp%rhoh = dh_pcc(1, j, k)
          ftp%d = d_pcc(1, j, k)
          call ftp_refresh_thermal_properties_from_DH(ftp)
          !NOTE: This cause issues on GPU - switched to the member-by-member operations
!          dm%fbcx_ftp(1, j, k) = ftp
!          dm%fbcx_ftp(3, j, k) = ftp
          dm%fbcx_ftp(1, j, k)%t = ftp%t
          dm%fbcx_ftp(1, j, k)%d = ftp%d
          dm%fbcx_ftp(1, j, k)%m = ftp%m
          dm%fbcx_ftp(1, j, k)%k = ftp%k
          dm%fbcx_ftp(1, j, k)%h = ftp%h
          dm%fbcx_ftp(1, j, k)%rhoh = ftp%rhoh
          dm%fbcx_ftp(1, j, k)%cp = ftp%cp
          dm%fbcx_ftp(1, j, k)%b = ftp%b
          dm%fbcx_ftp(1, j, k)%alpha = ftp%alpha
          dm%fbcx_ftp(1, j, k)%Pr = ftp%Pr
          dm%fbcx_ftp(3, j, k)%t = ftp%t
          dm%fbcx_ftp(3, j, k)%d = ftp%d
          dm%fbcx_ftp(3, j, k)%m = ftp%m
          dm%fbcx_ftp(3, j, k)%k = ftp%k
          dm%fbcx_ftp(3, j, k)%h = ftp%h
          dm%fbcx_ftp(3, j, k)%rhoh = ftp%rhoh
          dm%fbcx_ftp(3, j, k)%cp = ftp%cp
          dm%fbcx_ftp(3, j, k)%b = ftp%b
          dm%fbcx_ftp(3, j, k)%alpha = ftp%alpha
          dm%fbcx_ftp(3, j, k)%Pr = ftp%Pr
        end do
      end do
      !$acc end parallel loop
    end if

    if(dm%ibcx_Tm(2) == IBC_NEUMANN .and. &
       dm%dpcc%xen(1) == dm%np(1)) then
      ny = size(dm%fbcx_ftp, 2)
      nz = size(dm%fbcx_ftp, 3)
      !$acc parallel loop collapse(2) private(ftp) default(present)
      do k = 1, nz
        do j = 1, ny
          ftp%rhoh = dh_pcc(dm%np(1), j, k)
          ftp%d = d_pcc(dm%np(1), j, k)
          call ftp_refresh_thermal_properties_from_DH(ftp)
          !NOTE: This cause issues on GPU - switched to the member-by-member operations
!          dm%fbcx_ftp(2, j, k) = ftp
!          dm%fbcx_ftp(4, j, k) = ftp
          dm%fbcx_ftp(2, j, k)%t = ftp%t
          dm%fbcx_ftp(2, j, k)%d = ftp%d
          dm%fbcx_ftp(2, j, k)%m = ftp%m
          dm%fbcx_ftp(2, j, k)%k = ftp%k
          dm%fbcx_ftp(2, j, k)%h = ftp%h
          dm%fbcx_ftp(2, j, k)%rhoh = ftp%rhoh
          dm%fbcx_ftp(2, j, k)%cp = ftp%cp
          dm%fbcx_ftp(2, j, k)%b = ftp%b
          dm%fbcx_ftp(2, j, k)%alpha = ftp%alpha
          dm%fbcx_ftp(2, j, k)%Pr = ftp%Pr
          dm%fbcx_ftp(4, j, k)%t = ftp%t
          dm%fbcx_ftp(4, j, k)%d = ftp%d
          dm%fbcx_ftp(4, j, k)%m = ftp%m
          dm%fbcx_ftp(4, j, k)%k = ftp%k
          dm%fbcx_ftp(4, j, k)%h = ftp%h
          dm%fbcx_ftp(4, j, k)%rhoh = ftp%rhoh
          dm%fbcx_ftp(4, j, k)%cp = ftp%cp
          dm%fbcx_ftp(4, j, k)%b = ftp%b
          dm%fbcx_ftp(4, j, k)%alpha = ftp%alpha
          dm%fbcx_ftp(4, j, k)%Pr = ftp%Pr
        end do
      end do
      !$acc end parallel loop
    end if
  end if
!----------------------------------------------------------------------------------------------------------
!  BC - y
!----------------------------------------------------------------------------------------------------------
  if( dm%ibcy_Tm(1) == IBC_NEUMANN .or. &
      dm%ibcy_Tm(2) == IBC_NEUMANN) then
    dh_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))     => fl%wk1
    dh_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
    d_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))      => fl%wk3
    d_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))  => fl%wk4
    call transpose_x_to_y(tm%rhoh, dh_ypencil, dm%dccc)
    call Get_y_midp_C2P_3D(dh_ypencil, dh_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp) ! exterpolation, check
    call axis_estimating_radial_xpx(dh_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(1))
    call transpose_x_to_y(dens, d_ypencil, dm%dccc)
    call Get_y_midp_C2P_3D(d_ypencil, d_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp) ! exterpolation, check
    call axis_estimating_radial_xpx(d_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(1))
    
    if(dm%ibcy_Tm(1) == IBC_NEUMANN .and. &
       dm%dcpc%yst(2) == 1) then
      nx = size(dm%fbcy_ftp, 1)
      nz = size(dm%fbcy_ftp, 3)
      !$acc parallel loop collapse(2) private(ftp) default(present)
      do k = 1, nz
        do i = 1, nx
          ftp%rhoh = dh_cpc_ypencil(i, 1, k)
          ftp%d = d_cpc_ypencil(i, 1, k)
          call ftp_refresh_thermal_properties_from_DH(ftp)
!          dm%fbcy_ftp(i, 1, k) = ftp
!          dm%fbcy_ftp(i, 3, k) = ftp
          dm%fbcy_ftp(i, 1, k)%t = ftp%t
          dm%fbcy_ftp(i, 1, k)%d = ftp%d
          dm%fbcy_ftp(i, 1, k)%m = ftp%m
          dm%fbcy_ftp(i, 1, k)%k = ftp%k
          dm%fbcy_ftp(i, 1, k)%h = ftp%h
          dm%fbcy_ftp(i, 1, k)%rhoh = ftp%rhoh
          dm%fbcy_ftp(i, 1, k)%cp = ftp%cp
          dm%fbcy_ftp(i, 1, k)%b = ftp%b
          dm%fbcy_ftp(i, 1, k)%alpha = ftp%alpha
          dm%fbcy_ftp(i, 1, k)%Pr = ftp%Pr
          dm%fbcy_ftp(i, 3, k)%t = ftp%t
          dm%fbcy_ftp(i, 3, k)%d = ftp%d
          dm%fbcy_ftp(i, 3, k)%m = ftp%m
          dm%fbcy_ftp(i, 3, k)%k = ftp%k
          dm%fbcy_ftp(i, 3, k)%h = ftp%h
          dm%fbcy_ftp(i, 3, k)%rhoh = ftp%rhoh
          dm%fbcy_ftp(i, 3, k)%cp = ftp%cp
          dm%fbcy_ftp(i, 3, k)%b = ftp%b
          dm%fbcy_ftp(i, 3, k)%alpha = ftp%alpha
          dm%fbcy_ftp(i, 3, k)%Pr = ftp%Pr
        end do
      end do
      !$acc end parallel loop
    end if
    if(dm%ibcy_Tm(2) == IBC_NEUMANN .and. &
       dm%dcpc%yen(2) == dm%np(2)) then
      nx = size(dm%fbcy_ftp, 1)
      nz = size(dm%fbcy_ftp, 3)
      !$acc parallel loop collapse(2) private(ftp) default(present)
      do k = 1, nz
        do i = 1, nx
          ftp%rhoh = dh_cpc_ypencil(i, dm%np(2), k)
          ftp%d = d_cpc_ypencil(i, dm%np(2), k)
          call ftp_refresh_thermal_properties_from_DH(ftp)
!          dm%fbcy_ftp(i, 2, k) = ftp
!          dm%fbcy_ftp(i, 4, k) = ftp
          dm%fbcy_ftp(i, 2, k)%t = ftp%t
          dm%fbcy_ftp(i, 2, k)%d = ftp%d
          dm%fbcy_ftp(i, 2, k)%m = ftp%m
          dm%fbcy_ftp(i, 2, k)%k = ftp%k
          dm%fbcy_ftp(i, 2, k)%h = ftp%h
          dm%fbcy_ftp(i, 2, k)%rhoh = ftp%rhoh
          dm%fbcy_ftp(i, 2, k)%cp = ftp%cp
          dm%fbcy_ftp(i, 2, k)%b = ftp%b
          dm%fbcy_ftp(i, 2, k)%alpha = ftp%alpha
          dm%fbcy_ftp(i, 2, k)%Pr = ftp%Pr
          dm%fbcy_ftp(i, 4, k)%t = ftp%t
          dm%fbcy_ftp(i, 4, k)%d = ftp%d
          dm%fbcy_ftp(i, 4, k)%m = ftp%m
          dm%fbcy_ftp(i, 4, k)%k = ftp%k
          dm%fbcy_ftp(i, 4, k)%h = ftp%h
          dm%fbcy_ftp(i, 4, k)%rhoh = ftp%rhoh
          dm%fbcy_ftp(i, 4, k)%cp = ftp%cp
          dm%fbcy_ftp(i, 4, k)%b = ftp%b
          dm%fbcy_ftp(i, 4, k)%alpha = ftp%alpha
          dm%fbcy_ftp(i, 4, k)%Pr = ftp%Pr
        end do
      end do
      !$acc end parallel loop
    end if
  end if
!----------------------------------------------------------------------------------------------------------
!  BC - z
!----------------------------------------------------------------------------------------------------------
  if( dm%ibcz_Tm(1) == IBC_NEUMANN .or. &
      dm%ibcz_Tm(2) == IBC_NEUMANN) then
    dh_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))     => fl%wk1
    dh_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))     => fl%wk2
    dh_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk3
    call transpose_x_to_y(tm%rhoh, dh_ypencil, dm%dccc)
    call transpose_y_to_z(dh_ypencil, dh_zpencil, dm%dccc)
    call Get_z_midp_C2P_3D(dh_zpencil, dh_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp) ! exterpolation, check

    d_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))     => fl%wk1
    d_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))     => fl%wk2
    d_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk4
    call transpose_x_to_y(dens, d_ypencil, dm%dccc)
    call transpose_y_to_z(d_ypencil, d_zpencil, dm%dccc)
    call Get_z_midp_C2P_3D(d_zpencil, d_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp) ! exterpolation, check
    
    if(dm%ibcz_Tm(1) == IBC_NEUMANN .and. &
       dm%dccp%zst(1) == 1) then
      nx = size(dm%fbcz_ftp, 1)
      ny = size(dm%fbcz_ftp, 2)
      !$acc parallel loop collapse(2) private(ftp) default(present)
      do j = 1, ny
        do i = 1, nx
          ftp%rhoh = dh_ccp_zpencil(i, j, 1)
          ftp%d = d_ccp_zpencil(i, j, 1)
          call ftp_refresh_thermal_properties_from_DH(ftp)
!          dm%fbcz_ftp(i, j, 1) = ftp
!          dm%fbcz_ftp(i, j, 3) = ftp
          dm%fbcz_ftp(i, j, 1)%t = ftp%t
          dm%fbcz_ftp(i, j, 1)%d = ftp%d
          dm%fbcz_ftp(i, j, 1)%m = ftp%m
          dm%fbcz_ftp(i, j, 1)%k = ftp%k
          dm%fbcz_ftp(i, j, 1)%h = ftp%h
          dm%fbcz_ftp(i, j, 1)%rhoh = ftp%rhoh
          dm%fbcz_ftp(i, j, 1)%cp = ftp%cp
          dm%fbcz_ftp(i, j, 1)%b = ftp%b
          dm%fbcz_ftp(i, j, 1)%alpha = ftp%alpha
          dm%fbcz_ftp(i, j, 1)%Pr = ftp%Pr
          dm%fbcz_ftp(i, j, 3)%t = ftp%t
          dm%fbcz_ftp(i, j, 3)%d = ftp%d
          dm%fbcz_ftp(i, j, 3)%m = ftp%m
          dm%fbcz_ftp(i, j, 3)%k = ftp%k
          dm%fbcz_ftp(i, j, 3)%h = ftp%h
          dm%fbcz_ftp(i, j, 3)%rhoh = ftp%rhoh
          dm%fbcz_ftp(i, j, 3)%cp = ftp%cp
          dm%fbcz_ftp(i, j, 3)%b = ftp%b
          dm%fbcz_ftp(i, j, 3)%alpha = ftp%alpha
          dm%fbcz_ftp(i, j, 3)%Pr = ftp%Pr
        end do
      end do
      !$acc end parallel loop
    end if
    if(dm%ibcz_Tm(2) == IBC_NEUMANN .and. &
       dm%dccp%zen(1) == dm%np(3)) then
      nx = size(dm%fbcz_ftp, 1)
      ny = size(dm%fbcz_ftp, 2)
      !$acc parallel loop collapse(2) private(ftp) default(present)
      do j = 1, ny
        do i = 1, nx
          ftp%rhoh = dh_ccp_zpencil(i, j, dm%np(3))
          ftp%d = d_ccp_zpencil(i, j, dm%np(3))
          call ftp_refresh_thermal_properties_from_DH(ftp)
!          dm%fbcz_ftp(i, j, 2) = ftp
!          dm%fbcz_ftp(i, j, 4) = ftp
          dm%fbcz_ftp(i, j, 2)%t = ftp%t
          dm%fbcz_ftp(i, j, 2)%d = ftp%d
          dm%fbcz_ftp(i, j, 2)%m = ftp%m
          dm%fbcz_ftp(i, j, 2)%k = ftp%k
          dm%fbcz_ftp(i, j, 2)%h = ftp%h
          dm%fbcz_ftp(i, j, 2)%rhoh = ftp%rhoh
          dm%fbcz_ftp(i, j, 2)%cp = ftp%cp
          dm%fbcz_ftp(i, j, 2)%b = ftp%b
          dm%fbcz_ftp(i, j, 2)%alpha = ftp%alpha
          dm%fbcz_ftp(i, j, 2)%Pr = ftp%Pr
          dm%fbcz_ftp(i, j, 4)%t = ftp%t
          dm%fbcz_ftp(i, j, 4)%d = ftp%d
          dm%fbcz_ftp(i, j, 4)%m = ftp%m
          dm%fbcz_ftp(i, j, 4)%k = ftp%k
          dm%fbcz_ftp(i, j, 4)%h = ftp%h
          dm%fbcz_ftp(i, j, 4)%rhoh = ftp%rhoh
          dm%fbcz_ftp(i, j, 4)%cp = ftp%cp
          dm%fbcz_ftp(i, j, 4)%b = ftp%b
          dm%fbcz_ftp(i, j, 4)%alpha = ftp%alpha
          dm%fbcz_ftp(i, j, 4)%Pr = ftp%Pr
        end do
      end do
      !$acc end parallel loop
    end if
  end if

  return
  end subroutine Update_thermal_properties
!==========================================================================================================
  subroutine Calculate_energy_fractional_step(rhs0, rhs1, dtmp, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(inout) :: rhs0, rhs1
    integer,  intent(in) :: isub
    
    real(WP) :: rhs_explicit_current, rhs_explicit_last, rhs_total
    integer :: i, j, k
    integer :: nx, ny, nz

    nx = dtmp%xsz(1)
    ny = dtmp%xsz(2)
    nz = dtmp%xsz(3)
   !$acc  parallel loop collapse(3) default(present)  &
   !$acc& private(rhs_explicit_current, rhs_explicit_last, rhs_total)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

      ! add explicit terms : convection+viscous rhs
          rhs_explicit_current = rhs1(i, j, k) ! not (*dt)
          rhs_explicit_last    = rhs0(i, j, k) ! not (*dt)
          rhs_total = dm%tGamma(isub) * rhs_explicit_current + &
                      dm%tZeta (isub) * rhs_explicit_last
          rhs0(i, j, k) = rhs_explicit_current
      ! times the time step 
          rhs1(i, j, k) = dm%dt * rhs_total ! * dt
        end do
      end do
    end do
    !$acc end parallel loop

    return
  end subroutine
!==========================================================================================================
  subroutine Compute_energy_rhs(gx, gy, gz, tm, fl, dm, isub)
    use operations
    use udf_type_mod
    use thermo_info_mod
    use boundary_conditions_mod
    use wrt_debug_field_mod
    use cylindrical_rn_mod
    use wrt_debug_field_mod
    implicit none
    ! arguments
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in) :: isub
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ), intent(in) :: gx
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ), intent(in) :: gy
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ), intent(in) :: gz

    ! local variables
    real(WP), pointer, dimension(:,:,:) ::  accc_xpencil,           &
                                            apcc_xpencil,           &
                                            accc_ypencil,           &
                                            accp_ypencil,           &
                                            acpc_ypencil,           &
                                            accc_zpencil,           &
                                            accp_zpencil,           &

                                            gz_ccp_zpencil,         &
                                            gy_cpc_ypencil,         &

                                            hEnth_pcc_xpencil,      &
                                            hEnth_cpc_ypencil,      &
                                            hEnth_ccp_zpencil,      &

                                            Ttemp_ccc_ypencil,      &
                                            Ttemp_ccc_zpencil,      &

                                            kCond_pcc_xpencil,      &
                                            kCond_cpc_ypencil,      &
                                            kCond_ccp_zpencil,      &
                                            kCond_ccc_zpencil,      &

                                            fbcx_4cc,               &
                                            fbcy_c4c,               &
                                            fbcz_cc4

    integer, dimension(3) ::  ncccx, ncccy, ncccz, &
                              npccx, npccy, npccz, &
                              ncpcx, ncpcy, ncpcz, &
                              nccpx, nccpy, nccpz
    integer  :: i, j, k
    integer  :: nx, ny, nz
    integer  :: mbc(1:2, 1:3)

    ncccx = dm%dccc%xsz
    ncccy = dm%dccc%ysz
    ncccz = dm%dccc%zsz
    npccx = dm%dpcc%xsz
    npccy = dm%dpcc%ysz
    npccz = dm%dpcc%zsz
    ncpcx = dm%dcpc%xsz
    ncpcy = dm%dcpc%ysz
    ncpcz = dm%dcpc%zsz
    nccpx = dm%dccp%xsz
    nccpy = dm%dccp%ysz
    nccpz = dm%dccp%zsz

    !$acc kernels default(present)
    tm%ene_rhs = ZERO
    !$acc end kernels
!!------------------------------------------------------------------------------
!!  conv-x-e, x-pencil : d (gx * h_pcc) / dx 
!!------------------------------------------------------------------------------
    fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc1
    ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,4
      fbcx_4cc(i,j,k) = dm%fbcx_ftp(i,j,k)%h
    end do; end do; end do
    !$acc end parallel loop
    hEnth_pcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    call Get_x_midp_C2P_3D(tm%hEnth, hEnth_pcc_xpencil, dm, dm%iAccuracy, dm%ibcx_ftp(:), fbcx_4cc)
    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3))      => fl%wk2
    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      apcc_xpencil(i,j,k) = -gx(i,j,k) * hEnth_pcc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4cc, apcc_xpencil, dm%dpcc)
    else
      ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if
    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk1
    call Get_x_1der_P2C_3D(apcc_xpencil, accc_xpencil, dm, dm%iAccuracy, ebcx_conv, fbcx_4cc)

    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      tm%ene_rhs(i,j,k) = tm%ene_rhs(i,j,k) + accc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  conv-y-e, y-pencil : d (gy * h_cpc) / dy  * (1/r)
!!------------------------------------------------------------------------------
    gy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk1
    call transpose_x_to_y(gy, gy_cpc_ypencil, dm%dcpc)

    fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
    nx = ncpcy(1); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,4; do i=1,nx
      fbcy_c4c(i,j,k) = dm%fbcy_ftp(i,j,k)%h
    end do; end do; end do
    !$acc end parallel loop
    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))      => fl%wk2
    hEnth_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
    call transpose_x_to_y (tm%hEnth, accc_ypencil, dm%dccc)
    call Get_y_midp_C2P_3D(accc_ypencil, hEnth_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp(:), fbcy_c4c)
    call axis_estimating_radial_xpx(hEnth_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(1))

    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))      => fl%wk2
    nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      acpc_ypencil(i,j,k) = - gy_cpc_ypencil(i,j,k) * hEnth_cpc_ypencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    if(is_fbcy_velo_required) then
      call extract_dirichlet_fbcy(fbcy_c4c, acpc_ypencil, dm%dcpc, dm, is_reversed = .true.)
    else
      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if
    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))     => fl%wk1
    call Get_y_1der_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, ebcy_conv, fbcy_c4c)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accc_ypencil, dm%dccc, dm%rci, 1, IPENCIL(2))

    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3))     => fl%wk2
    call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      tm%ene_rhs(i,j,k) = tm%ene_rhs(i,j,k) + accc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!! conv-z-e, z-pencil : d (gz * h_ccp) / dz   * (1/r)
!!------------------------------------------------------------------------------
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))   => fl%wk1
    gz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk2
    call transpose_x_to_y(gz,        accp_ypencil,     dm%dccp)
    call transpose_y_to_z(accp_ypencil, gz_ccp_zpencil,   dm%dccp)

    fbcz_cc4(1:nccpz(1),1:nccpz(2),1:4) => fl%wkbc1
    nx = nccpz(1); ny = nccpz(2)
    !$acc parallel loop collapse(3) default(present)
    do k=1,4; do j=1,ny; do i=1,nx
      fbcz_cc4(i,j,k) = dm%fbcz_ftp(i,j,k)%h
    end do; end do; end do
    !$acc end parallel loop

    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))   => fl%wk1
    accc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))   => fl%wk3
    call transpose_x_to_y (tm%hEnth, accc_ypencil, dm%dccc)
    call transpose_y_to_z (accc_ypencil, accc_zpencil, dm%dccc)
    hEnth_ccp_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
    call Get_z_midp_C2P_3D(accc_zpencil, hEnth_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp(:), fbcz_cc4)

    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3))   => fl%wk3
    nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accp_zpencil(i,j,k) = - gz_ccp_zpencil(i,j,k) * hEnth_ccp_zpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    if(is_fbcz_velo_required) then
      call extract_dirichlet_fbcz(fbcz_cc4, accp_zpencil, dm%dccp)
    else
      nx = nccpz(1); ny = nccpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_cc4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
    call Get_z_1der_P2C_3D( accp_zpencil, accc_zpencil, dm, dm%iAccuracy, ebcz_conv, fbcz_cc4)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accc_zpencil, dm%dccc, dm%rci, 1, IPENCIL(3))

    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk3
    call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
    call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      tm%ene_rhs(i,j,k) = tm%ene_rhs(i,j,k) + accc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!! diff-x-e, d ( k_pcc * d (T) / dx ) dx
!!------------------------------------------------------------------------------
    fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc1
    call get_fbcx_iTh(dm%ibcx_Tm, dm, fbcx_4cc)
    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    call Get_x_1der_C2P_3D(tm%tTemp, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_Tm, fbcx_4cc)

    ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,4
      fbcx_4cc(i,j,k) = dm%fbcx_ftp(i,j,k)%k
    end do; end do; end do
    !$acc end parallel loop
    kCond_pcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
    call Get_x_midp_C2P_3D(tm%kCond, kCond_pcc_xpencil, dm, dm%iAccuracy, dm%ibcx_ftp(:), fbcx_4cc)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      apcc_xpencil(i,j,k) = apcc_xpencil(i,j,k) * kCond_pcc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4cc, apcc_xpencil, dm%dpcc)
    else
      ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if
    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk2
    call Get_x_1der_P2C_3D(apcc_xpencil, accc_xpencil, dm, dm%iAccuracy, ebcx_difu, fbcx_4cc)

    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      tm%ene_rhs(i,j,k) = tm%ene_rhs(i,j,k) + accc_xpencil(i,j,k) * tm%rPrRen
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!! diff-y-e, d ( r * k_cpc * d (T) / dy ) dy * 1/r
!!------------------------------------------------------------------------------
    fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
    call get_fbcy_iTh(dm%ibcy_Tm, dm, fbcy_c4c)
    Ttemp_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk1
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))      => fl%wk2
    call transpose_x_to_y (tm%Ttemp, Ttemp_ccc_ypencil, dm%dccc)
    call Get_y_1der_C2P_3D(Ttemp_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_Tm, fbcy_c4c)

    nx = ncpcy(1); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,4; do i=1,nx
      fbcy_c4c(i,j,k) = dm%fbcy_ftp(i,j,k)%k
    end do; end do; end do
    !$acc end parallel loop
    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))      => fl%wk1
    kCond_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
    call transpose_x_to_y (tm%kCond, accc_ypencil, dm%dccc)
    call Get_y_midp_C2P_3D(accc_ypencil,  kCond_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp(:), fbcy_c4c)
    call axis_estimating_radial_xpx(kCond_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(1))

    nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      acpc_ypencil(i,j,k) = acpc_ypencil(i,j,k) * kCond_cpc_ypencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(acpc_ypencil, dm%dcpc, dm%rp, 1, IPENCIL(2))
    if(is_fbcy_velo_required) then
      call extract_dirichlet_fbcy(fbcy_c4c, acpc_ypencil, dm%dcpc, dm, is_reversed = .true.)
    else
      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if
    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk1
    call Get_y_1der_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, ebcy_difu, fbcy_c4c)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accc_ypencil, dm%dccc, dm%rci, 1, IPENCIL(2))

    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk2
    call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      tm%ene_rhs(i,j,k) = tm%ene_rhs(i,j,k) + accc_xpencil(i,j,k) * tm%rPrRen
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!! diff-z-e, d (1/r* k_ccp * d (T) / dz ) / dz * 1/r
!!------------------------------------------------------------------------------
    fbcz_cc4(1:nccpz(1),1:nccpz(2),1:4) => fl%wkbc1
    call get_fbcz_iTh(dm%ibcz_Tm, dm, fbcz_cc4)
    Ttemp_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
    Ttemp_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
    call transpose_x_to_y (tm%Ttemp, Ttemp_ccc_ypencil, dm%dccc)
    call transpose_y_to_z (Ttemp_ccc_ypencil, Ttemp_ccc_zpencil, dm%dccc)
    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3))      => fl%wk2
    call Get_z_1der_C2P_3D(Ttemp_ccc_zpencil, accp_zpencil, dm, dm%iAccuracy, dm%ibcz_Tm, fbcz_cc4)

    nx = nccpz(1); ny = nccpz(2)
    !$acc parallel loop collapse(3) default(present)
    do k=1,4; do j=1,ny; do i=1,nx
      fbcz_cc4(i,j,k) = dm%fbcz_ftp(i,j,k)%k
    end do; end do; end do
    !$acc end parallel loop
    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))      => fl%wk1
    kCond_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk3
    call transpose_x_to_y (tm%kCond, accc_ypencil, dm%dccc)
    call transpose_y_to_z (accc_ypencil, kCond_ccc_zpencil, dm%dccc)
    kCond_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk1
    call Get_z_midp_C2P_3D(kCond_ccc_zpencil, kCond_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp(:), fbcz_cc4)

    nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accp_zpencil(i,j,k) = accp_zpencil(i,j,k) * kCond_ccp_zpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accp_zpencil, dm%dccp, dm%rci, 1, IPENCIL(3))
    if(is_fbcz_velo_required) then
      call extract_dirichlet_fbcz(fbcz_cc4, accp_zpencil, dm%dccp)
    else
      nx = nccpz(1); ny = nccpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_cc4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
    call Get_z_1der_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, ebcz_difu, fbcz_cc4)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accc_zpencil, dm%dccc, dm%rci, 1, IPENCIL(3))

    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk3
    call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
    call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      tm%ene_rhs(i,j,k) = tm%ene_rhs(i,j,k) + accc_xpencil(i,j,k) * tm%rPrRen
    end do; end do; end do
    !$acc end parallel loop

    call Calculate_energy_fractional_step(tm%ene_rhs0, tm%ene_rhs, dm%dccc, dm, isub)

    return
  end subroutine Compute_energy_rhs

!==========================================================================================================
!==========================================================================================================
  subroutine Solve_energy_eq(fl, tm, dm, isub)
    use udf_type_mod
    use thermo_info_mod 
    use solver_tools_mod
    use boundary_conditions_mod
    use bc_convective_outlet_mod
    use convert_primary_conservative_mod
    implicit none
    ! arguments
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in)    :: isub

    ! backup density every RK stage
    !$acc kernels default(present)
    fl%dDens0(:,:,:) = fl%dDens(:,:,:)
    !$acc end kernels
    ! compute b.c. info from convective b.c. if specified.
    call update_convective_outlet_thermo(tm, dm, isub)
    ! calculate rhs of energy equation
    call Compute_energy_rhs(fl%gx, fl%gy, fl%gz, tm, fl, dm, isub)
    !  update rho * h
    !$acc kernels default(present)
    tm%rhoh(:,:,:) = tm%rhoh(:,:,:) + tm%ene_rhs(:,:,:)
    !$acc end kernels
    !  update other properties from rho * h for domain + b.c.
    call Update_thermal_properties(fl%dDens, fl%mVisc, tm, fl, dm)
    if (dm%icase == ICASE_PIPE) call update_fbcy_cc_thermo_halo(tm, dm)

    return
  end subroutine

end module eq_energy_mod
