module eq_energy_mod
  use operations
  use decomp_2d
  implicit none

  private :: Compute_energy_rhs
  private :: Calculate_energy_fractional_step
  public  :: Solve_energy_eq
contains
!===============================================================================
  subroutine Calculate_energy_fractional_step(rhs0, rhs1, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    real(WP), dimension(:, :, :), intent(inout) :: rhs0, rhs1
    integer,                   intent(in   ) :: isub
    type(t_domain),               intent(in   ) :: dm

    real(WP), dimension(size(rhs1, 1), size(rhs1, 2), size(rhs1, 3)) :: rhs_dummy


  ! add explicit terms
    rhs_dummy(:, :, :) = rhs1(:, :, :)
    rhs1(:, :, :) = dm%tGamma(isub) * rhs1(:, :, :) + &
                    dm%tZeta (isub) * rhs0(:, :, :)
    rhs0(:, :, :) = rhs_dummy(:, :, :)

  ! times the time step 
    rhs1(:, :, :) = dm%dt * rhs1(:, :, :)

    return
  end subroutine
!===============================================================================
  subroutine Compute_energy_rhs(fl, tm, dm, isub)
    use operations
    use udf_type_mod
    use thermo_info_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in) :: isub    

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: gy_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: gz_zpencil 

    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dppc%xsz(3) ) :: hEnth_xpcc
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dppc%xsz(3) ) :: kCond_xpcc
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: hEnth_ycpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: hEnth_zccp_zpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: Ttemp_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: ene_rhs_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: kCond_ycpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: kCond_zccp_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: Ttemp_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: kCond_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: ene_rhs_zpencil
    
    real(WP) :: fbcx(2), fbcy(2), fbcz(2)
    integer  :: ibcx(2), ibcy(2), ibcz(2)
    integer  :: i
!===============================================================================
!   preparation
!===============================================================================
    call transpose_x_to_y(fl%gy,        gy_ypencil,   dm%dcpc)   ! for d(g_y h)/dy
    call transpose_x_to_y(fl%gz,        accp_ypencil, dm%dccp)   ! intermediate, accp_ypencil = gz_ypencil
    call transpose_y_to_z(accp_ypencil, gz_zpencil,   dm%dccp)   ! for d(g_z h)/dz
!-------------------------------------------------------------------------------
!    h --> h_xpcc
!      --> h_ypencil --> h_ycpc_ypencil
!                    --> h_zpencil --> h_zccp_zpencil
!-------------------------------------------------------------------------------
    do i = 1, 2
      fbcx(i) = tm%tpbcx(i)%h
      fbcy(i) = tm%tpbcy(i)%h
      fbcz(i) = tm%tpbcz(i)%h
    end do
    call Get_x_midp_C2P_3D(tm%hEnth,     hEnth_xpcc,         dm, dm%ibcx(:, 5), fbcx(:)) ! for d(g_x h_pcc))/dy
    call transpose_x_to_y (tm%hEnth,     accc_ypencil, dm%dccc)                     !intermediate, accc_ypencil = hEnth_ypencil
    call Get_y_midp_C2P_3D(accc_ypencil, hEnth_ycpc_ypencil, dm, dm%ibcy(:, 5), fbcy(:))! for d(g_y h_cpc)/dy
    call transpose_y_to_z (accc_ypencil, accc_zpencil, dm%dccc) !intermediate, accc_zpencil = hEnth_zpencil
    call Get_z_midp_C2P_3D(accc_zpencil, hEnth_zccp_zpencil, dm, dm%ibcz(:, 5), fbcz(:)) ! for d(g_z h_ccp)/dz
!-------------------------------------------------------------------------------
!    k --> k_xpcc
!      --> k_ypencil --> k_ycpc_ypencil
!                    --> k_zpencil --> k_zccp_zpencil              
!-------------------------------------------------------------------------------
    do i = 1, 2
      fbcx(i) = tm%tpbcx(i)%k
      fbcy(i) = tm%tpbcy(i)%k
      fbcz(i) = tm%tpbcz(i)%k
    end do
    call Get_x_midp_C2P_3D(tm%kCond,      kCond_xpcc,     dm, dm%ibcx(:, 5), fbcx(:) ) ! for d(k_pcc * (dT/dx) )/dx
    call transpose_x_to_y (tm%kCond,      accc_ypencil, dm%dccc)  ! for k d2(T)/dy^2
    call Get_y_midp_C2P_3D(accc_ypencil,  kCond_ycpc_ypencil, dm, dm%ibcy(:, 5), fbcy(:))
    call transpose_x_to_y (accc_ypencil,  kCond_zpencil, dm%dccc) 
    call Get_z_midp_C2P_3D(kCond_zpencil, kCond_zccp_zpencil, dm, dm%ibcz(:, 5), fbcz(:))
!-------------------------------------------------------------------------------
!    T --> T_ypencil --> T_zpencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y (tm%Ttemp,      Ttemp_ypencil, dm%dccc)   ! for k d2(T)/dy^2
    call transpose_x_to_y (Ttemp_ypencil, Ttemp_zpencil, dm%dccc)   ! for k d2(T)/dz^2
!===============================================================================
! the RHS of energy equation
! x-pencil : the RHS terms of energy (derivative) operating in the x direction
!===============================================================================
!-------------------------------------------------------------------------------
! x-pencil : d (gx * h_pcc) / dx 
!-------------------------------------------------------------------------------
    tm%ene_rhs = ZERO
    call Get_x_1st_derivative_P2C_3D( - fl%gx * hEnth_xpcc, accc, dm, dm%ibcx(:, 1) ) ! accc = -d(gx * h)/dx
    tm%ene_rhs = tm%ene_rhs + accc
!-------------------------------------------------------------------------------
! x-pencil : d (T) / dx 
!-------------------------------------------------------------------------------
    do i = 1, 2
      if (dm%ibcx(i, 5) == IBC_NEUMANN) then
        ibcx(i) = IBC_INTERIOR
        fbcx(i) = ZERO
      else
        ibcx(i) = dm%ibcx(i, 5)
        fbcx(i) = dm%fbcx(i, 5)
      end if
    end do
    call Get_x_1st_derivative_C2P_3D(tm%tTemp, apcc, dm, ibcx(:), fbcx(:) )
!-------------------------------------------------------------------------------
! x-pencil : k_pcc * d (T) / dx 
!-------------------------------------------------------------------------------
    apcc = apcc * kCond_xpcc
    if (dm%ibcx(1, 5) == IBC_NEUMANN) then
      apcc(1, :, :) = dm%fbcx(1, 5)
    end if
    if (dm%ibcx(2, 5) == IBC_NEUMANN) then
      apcc(dm%dpcc%xen(1), :, :) = dm%fbcx(2, 5)
    end if
!-------------------------------------------------------------------------------
! x-pencil : d ( k_pcc * d (T) / dx ) dx
!-------------------------------------------------------------------------------
    call Get_x_1st_derivative_P2C_3D(apcc, accc, dm, dm%ibcx(:, 5) )

    tm%ene_rhs = tm%ene_rhs + accc
!===============================================================================
! the RHS of energy equation
! y-pencil : the RHS terms of energy (derivative) operating in the y direction
!===============================================================================
!-------------------------------------------------------------------------------
! y-pencil : d (gy * h_cpc) / dy 
!-------------------------------------------------------------------------------
    ene_rhs_ypencil = ZERO
    call Get_y_1st_derivative_P2C_3D( - gy_ypencil * hEnth_ycpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 2) )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil
!-------------------------------------------------------------------------------
! y-pencil : d (T) / dy
!-------------------------------------------------------------------------------
    do i = 1, 2
      if (dm%ibcy(i, 5) == IBC_NEUMANN) then
        ibcy(i) = IBC_INTERIOR
        fbcy(i) = ZERO
      else
        ibcy(i) = dm%ibcy(i, 5)
        fbcy(i) = dm%fbcy(i, 5)
      end if
    end do
    call Get_y_1st_derivative_C2P_3D(tTemp_ypencil, acpc_ypencil, dm, ibcy(:), fbcy(:) )
!-------------------------------------------------------------------------------
! y-pencil : k_cpc * d (T) / dy 
!-------------------------------------------------------------------------------
    acpc_ypencil = acpc_ypencil * kCond_ycpc_ypencil
    if (dm%ibcy(1, 5) == IBC_NEUMANN) then
      acpc_ypencil(:, 1, :) = dm%fbcy(1, 5)
    end if
    if (dm%ibcy(2, 5) == IBC_NEUMANN) then
      acpc_ypencil(:, dm%dcpc%yen(2), :) = dm%fbcy(2, 5)
    end if
!-------------------------------------------------------------------------------
! y-pencil : d ( k_cpc * d (T) / dy ) dy
!-------------------------------------------------------------------------------
    call Get_y_1st_derivative_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 5) )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil

    call transpose_y_to_x(ene_rhs_ypencil, accc, dm%dccc)
    tm%ene_rhs = tm%ene_rhs + accc
!===============================================================================
! the RHS of energy equation
! z-pencil : the RHS terms of energy (derivative) operating in the z direction
!===============================================================================
!-------------------------------------------------------------------------------
! z-pencil : d (gz * h_ccp) / dz 
!-------------------------------------------------------------------------------
    ene_rhs_zpencil = ZERO
    call Get_z_1st_derivative_P2C_3D( - gz_zpencil * hEnth_zccp_zpencil, accc_zpencil, dm, dm%ibcz(:, 3) )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil
!-------------------------------------------------------------------------------
! z-pencil : d (T) / dz
!-------------------------------------------------------------------------------
    do i = 1, 2
      if (dm%ibcz(i, 5) == IBC_NEUMANN) then
        ibcz(i) = IBC_INTERIOR
        fbcz(i) = ZERO
      else
        ibcz(i) = dm%ibcz(i, 5)
        fbcz(i) = dm%fbcz(i, 5)
      end if
    end do
    call Get_z_1st_derivative_C2P_3D(tTemp_zpencil, accp_zpencil, dm, ibcz(:), fbcz(:) )
!-------------------------------------------------------------------------------
! z-pencil : k_ccp * d (T) / dz 
!-------------------------------------------------------------------------------
    accp_zpencil = accp_zpencil * kCond_zccp_zpencil
    if (dm%ibcz(1, 5) == IBC_NEUMANN) then
      accp_zpencil(:, 1, :) = dm%fbcz(1, 5)
    end if
    if (dm%ibcz(2, 5) == IBC_NEUMANN) then
      accp_zpencil(:, :, dm%dccp%zen(3)) = dm%fbcz(2, 5)
    end if
!-------------------------------------------------------------------------------
! z-pencil : d ( k_ccp * d (T) / dz ) / dz
!-------------------------------------------------------------------------------
    call Get_z_1st_derivative_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%ibcz(:, 5) )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil

    call transpose_z_to_y(ene_rhs_zpencil, ene_rhs_ypencil, dm%dccc)
    call transpose_y_to_x(ene_rhs_ypencil, accc,            dm%dccc)
    tm%ene_rhs = tm%ene_rhs + accc

!===============================================================================
! time approaching
!===============================================================================
    call Calculate_energy_fractional_step(tm%ene_rhs0, tm%ene_rhs, dm, isub)

    return
  end subroutine Compute_energy_rhs
!===============================================================================
  subroutine Update_thermal_properties_from_dh(fl, tm, dm)
    use udf_type_mod
    use thermo_info_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    integer :: i, j, k
    type(t_thermoProperty) :: tp
!-------------------------------------------------------------------------------
!   x-pencil
!-------------------------------------------------------------------------------
    do k = dm%dccc%xst(3), dm%dccc%xen(3)
      do j = dm%dccc%xst(2), dm%dccc%xen(2)
        do i = dm%dccc%xst(1), dm%dccc%xen(1)
          tp%dh = tm%dh(i, j, k)
          call tp%Refresh_thermal_properties_from_DH()
          tm%hEnth(i, j, k) = tp%h
          tm%tTemp(i, j, k) = tp%T
          tm%kCond(i, j, k) = tp%k
          fl%dDens(i, j, k) = tp%d
          fl%mVisc(i, j, k) = tp%m
        end do
      end do
    end do

  return
  end subroutine Update_thermal_properties_from_dh
!===============================================================================
  subroutine Solve_energy_eq(fl, tm, dm, isub)
    use udf_type_mod
    use thermo_info_mod 
    implicit none
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in)    :: isub

!-------------------------------------------------------------------------------
!   calculate rhs of energy equation
!-------------------------------------------------------------------------------
    call Compute_energy_rhs(fl, tm, dm, isub)
!-------------------------------------------------------------------------------
!   update rho * h
!-------------------------------------------------------------------------------
    tm%dh = tm%dh + tm%ene_rhs
!-------------------------------------------------------------------------------
!   update other properties from rho * h
!-------------------------------------------------------------------------------
    call Update_thermal_properties_from_dh(fl, tm, dm)
!-------------------------------------------------------------------------------
!   No Need to apply b.c.
!-------------------------------------------------------------------------------
  return
  end subroutine

end module eq_energy_mod
