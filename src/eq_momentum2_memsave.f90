module eq_momentum_mod
  use operations
  use precision_mod
  use decomp_2d
  use print_msg_mod
  use wrt_debug_field_mod
  implicit none

  private :: Calculate_momentum_fractional_step
  private :: Compute_momentum_rhs
  private :: Correct_massflux
  private :: solve_pressure_poisson
  private :: gravity_decomposition_to_rz
  !private :: solve_poisson_x2z
  
  public  :: Solve_momentum_eq

contains
!==========================================================================================================
  subroutine gravity_decomposition_to_rz(dens, iforce, fgravity, fr_cpc_ypencil, ft_ccp_zpencil, &
                                         wk1, wk2, wkbc1, dm)
    use math_mod
    use udf_type_mod
    use cylindrical_rn_mod
    implicit none
    type(t_domain), intent(in) :: dm
    INTEGER, intent(in) :: iforce
    real(WP), intent(in) :: fgravity
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(in)  :: dens
    real(WP), dimension(dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3)), intent(out) :: fr_cpc_ypencil
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3)), intent(out) :: ft_ccp_zpencil
    real(WP), pointer, dimension(:,:,:), contiguous :: wk1, wk2, wkbc1

    real(WP), dimension(dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3)) :: fr_cpp_zpencil, &
                                                                           ft_cpp_zpencil
    integer, dimension(3) :: ncccy, ncpcy, ncpcz, nccpy, ncppy, ncppz
    real(WP), pointer, dimension(:,:,:) :: accc_ypencil, acpc_ypencil, &
                                           accp_ypencil, acpp_ypencil, &
                                           acpc_zpencil, acpp_zpencil, &
                                           fbcy_c4c
    integer  :: i, k, j
    integer  :: nx, ny, nz
    real(WP) :: theta

    if(dm%icoordinate /= ICYLINDRICAL) return

    ncccy = dm%dccc%ysz
    ncpcy = dm%dcpc%ysz
    ncpcz = dm%dcpc%zsz
    nccpy = dm%dccp%ysz
    ncppy = dm%dcpp%ysz
    ncppz = dm%dcpp%zsz

!-------------------------------------------------------------------------------
!   density interpolation from ccc to cpp
!-------------------------------------------------------------------------------
    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => wk1
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
    call transpose_x_to_y(dens, accc_ypencil, dm%dccc)

    fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc1
    nx = ncpcy(1); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1, nz; do j=1, 4; do i=1, nx
      fbcy_c4c(i, j, k) = dm%fbcy_ftp(i, j, k)%d
    end do; end do; end do
    !$acc end parallel loop

    call Get_y_midp_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)

#ifdef DEBUG_STEPS  
    if(dm%icase == ICASE_PIPE) then
      write(*,*) 'density', acpc_ypencil(4, 1, 4),  acpc_ypencil(4, 1, dm%knc_sym(4)) , &
                  acpc_ypencil(4, 1, 4)-acpc_ypencil(4, 1, dm%knc_sym(4))
    end if
#endif

    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
    call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)
    acpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => wk2
    call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp)
!-------------------------------------------------------------------------------
!   force: calculate in z-pencil
!-------------------------------------------------------------------------------
    !$acc data create(fr_cpp_zpencil, ft_cpp_zpencil)
    if(iforce == 2 .or. iforce == -2) then
      nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
      !$acc parallel loop collapse(3) default(present) private(theta)
      do k=1, nz; do j=1, ny; do i=1, nx
        theta = dm%h(3) * REAL(k, WP)
        fr_cpp_zpencil(i, j, k) =  acpp_zpencil(i, j, k) * fgravity * sin_wp(theta)
        ft_cpp_zpencil(i, j, k) = -acpp_zpencil(i, j, k) * fgravity * cos_wp(theta)
      end do; end do; end do
      !$acc end parallel loop
    else if(iforce == 3 .or. iforce == -3) then
      nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
      !$acc parallel loop collapse(3) default(present) private(theta)
      do k=1, nz; do j=1, ny; do i=1, nx
        theta = dm%h(3) * REAL(k, WP)
        fr_cpp_zpencil(i, j, k) =  acpp_zpencil(i, j, k) * fgravity * cos_wp(theta)
        ft_cpp_zpencil(i, j, k) =  acpp_zpencil(i, j, k) * fgravity * sin_wp(theta)
      end do; end do; end do
      !$acc end parallel loop
    else
      call Print_error_msg('gravity direction is not correct for cylindrical coordinates')
    end if
!-------------------------------------------------------------------------------
!   back to momentum eq.
!-------------------------------------------------------------------------------
    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
    call Get_z_midp_P2C_3D(fr_cpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp)
    call transpose_z_to_y(acpc_zpencil, fr_cpc_ypencil, dm%dcpc)
    call multiple_cylindrical_rn(fr_cpc_ypencil, dm%dcpc, dm%rp, 1, IPENCIL(2))
    acpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => wk2
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk1
    call transpose_z_to_y(ft_cpp_zpencil, acpp_ypencil, dm%dcpp)
    call Get_y_midp_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp)
    call transpose_y_to_z(accp_ypencil, ft_ccp_zpencil, dm%dccp)
    !$acc end data

    return

  end subroutine
!==========================================================================================================
!==========================================================================================================
!> \brief To calcuate the convection and diffusion terms in rhs of momentum eq.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  rhs0          the last iteration rhs
!> \param[inout]  rhs1          the current iteration rhs
!> \param[in]     rhs1_semi     the semi-implicit term
!> \param[in]     isub          the RK iteration to get correct Coefficient 
!_______________________________________________________________________________
  subroutine Calculate_momentum_fractional_step(rhs0, rhs1, rhs1_pfc, dtmp, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    type(t_domain), intent(in) :: dm
    real(WP), dimension(:, :, :), intent(inout) :: rhs1_pfc
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
          rhs_explicit_current = rhs1(i, j, k) ! not (* dt)
          rhs_explicit_last    = rhs0(i, j, k) ! not (* dt)
          rhs_total = dm%tGamma(isub) * rhs_explicit_current + &
                      dm%tZeta (isub) * rhs_explicit_last
          rhs0(i, j, k) = rhs_explicit_current ! not (* dt)
      ! add pressure gradient
          rhs_total = rhs_total + &
                      dm%tAlpha(isub) * rhs1_pfc(i, j, k)

      ! times the time step 
          rhs1(i, j, k) = dm%dt * rhs_total ! * dt 
        end do
      end do
    end do
    !$acc end parallel loop

    return
  end subroutine
!==========================================================================================================
!> \brief To calcuate all rhs of momentum eq.
!==========================================================================================================
  subroutine Compute_momentum_rhs(fl, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    use operations
    use solver_tools_mod
    use typeconvert_mod
    use boundary_conditions_mod
    use wrt_debug_field_mod
    use find_max_min_ave_mod
    use cylindrical_rn_mod
    use math_mod
    implicit none
    ! arguments
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(inout) :: dm
    integer,     intent(in ) :: isub

    real(WP), pointer, dimension(:,:,:) ::  mx_rhs_ypencil,     &
                                            mx_rhs_zpencil,     &
                                            mx_rhs_pfc_xpencil, &
                                            my_rhs_ypencil,     &
                                            my_rhs_zpencil,     &
                                            my_rhs_pfc_xpencil, &
                                            my_rhs_pfc_ypencil, &
                                            mz_rhs_ypencil,     &
                                            mz_rhs_zpencil,     &
                                            mz_rhs_pfc_xpencil, &
                                            mz_rhs_pfc_ypencil, &
                                            mz_rhs_pfc_zpencil, &

                                            pres_ypencil,       &
                                            pres_zpencil,       &
                                            div_ccc_xpencil,    &
                                            div_ccc_ypencil,    &
                                            div_ccc_zpencil,    &

                                            qxix_ccc_xpencil,   &
                                            qxiy_ppc_ypencil,   &
                                            qxiy_ppc_xpencil,   &
                                            qxiz_pcp_zpencil,   &
                                            qxiz_pcp_xpencil,   &
                                            gxix_ccc_xpencil,   &
                                            gxiy_ppc_xpencil,   &
                                            gxiy_ppc_ypencil,   &
                                            gxiz_pcp_xpencil,   &
                                            qxdx_pcc_xpencil,   &
                                            qxdx_cpc_ypencil,   &
                                            qxdx_ccp_zpencil,   &
                                            qxdx_ccc_xpencil,   &
                                            qxdy_ppc_xpencil,   &
                                            qxdy_ppc_ypencil,   &
                                            qxdz_pcp_xpencil,   &
                                            qxdz_pcp_zpencil,   &

                                            qyix_ppc_xpencil,   &
                                            qyix_ppc_ypencil,   &
                                            qyiy_ccc_ypencil,   &
                                            qyiz_cpp_ypencil,   &
                                            qyiz_cpp_zpencil,   &
                                            gyix_ppc_ypencil,   &
                                            gyiy_ccc_ypencil,   &
                                            gyiz_cpp_ypencil,   &
                                            qyr_ypencil,        &
                                            qyr2_ypencil,       &
                                            qyriy_ccc_zpencil,  &
                                            qyriy_ccc_ypencil,  &
                                            qyriz_cpp_ypencil,  &
                                            qyriz_cpp_zpencil,  &
                                            qydx_ppc_xpencil,   &
                                            qydx_ppc_ypencil,   &
                                            qydz_cpp_ypencil,   &
                                            qydz_cpp_zpencil,   &
                                            qydy_ccc_ypencil,   &
                                            qydy_pcc_xpencil,   &
                                            qydy_cpc_ypencil,   &
                                            qydy_ccp_zpencil,   &
                                            qyrdy_ccc_ypencil,  &
                                            qyr2dz_cpp_zpencil, &
                                            qyrdz_cpp_ypencil,  &
                                            qyr2dz_cpp_ypencil, &

                                            qz_zpencil,         &
                                            qzix_pcp_xpencil,   &
                                            qzix_pcp_zpencil,   &
                                            qziy_cpp_ypencil,   &
                                            qziy_cpp_zpencil,   &
                                            qziz_ccc_zpencil,   &
                                            gzix_pcp_zpencil,   &
                                            gziy_cpp_zpencil,   &
                                            gziy_cpp_ypencil,   &
                                            gziz_ccc_zpencil,   &
                                            qziz_ccc_ypencil,   &
                                            qzriy_cpp_ypencil,  &
                                            qzriy_cpp_zpencil,  &
                                            gziz_ccc_ypencil,   &
                                            qzdx_pcp_xpencil,   &
                                            qzdx_pcp_zpencil,   &
                                            qzdy_cpp_ypencil,   &
                                            qzdy_cpp_zpencil,   &
                                            qzdz_ccc_ypencil,   &
                                            qzdz_ccc_zpencil,   &
                                            qzdz_pcc_xpencil,   &
                                            qzdz_cpc_ypencil,   &
                                            qzdz_ccp_zpencil,   &
                                            qzrdz_cpc_ypencil,  &

                                            mu_ccc_xpencil,     &
                                            mu_ccc_ypencil,     &
                                            mu_ccc_zpencil,     &
                                            muiz_ccp_zpencil,   &
                                            muixy_ppc_xpencil,  &
                                            muixy_ppc_ypencil,  &
                                            muixz_pcp_xpencil,  &
                                            muixz_pcp_zpencil,  &
                                            muiyz_cpp_ypencil,  &
                                            muiyz_cpp_zpencil,  &
                                            dDens_ypencil,      &
                                            dDens_zpencil,      &

                                            apcc_xpencil,       &
                                            apcc_xpencil1,      &
                                            apcc_xpencil2,      &
                                            apcc_xpencil3,      &
                                            accc_xpencil,       &
                                            appc_xpencil,       &
                                            appc_xpencil1,      &
                                            apcp_xpencil,       &
                                            appc_ypencil,       &
                                            accc_ypencil,       &
                                            accc_ypencil1,      &
                                            acpp_ypencil,       &
                                            acpp_ypencil1,      &
                                            accp_ypencil,       &
                                            acpc_ypencil,       &
                                            acpc_ypencil1,      &
                                            apcc_ypencil,       &
                                            apcp_ypencil,       &
                                            apcp_zpencil,       &
                                            acpp_zpencil,       &
                                            acpp_zpencil1,      &
                                            accc_zpencil,       &
                                            accc_zpencil1,      &
                                            accp_zpencil,       &
                                            accp_zpencil1,      &
                                            accp_zpencil2,      &
                                            accp_zpencil3,      &
                                            acpc_zpencil,       &
                                            acpc_zpencil1,      &
                                            apcc_zpencil,       &
                                            accp_xpencil,       &
                                            acpc_xpencil,       &

                                            fbcx_4cc,           &
                                            fbcx_4cc1,          &
                                            fbcx_div_4cc,       &
                                            fbcx_mu_4cc,        &
                                            fbcy_c4c,           &
                                            fbcy_c4c1,          &
                                            fbcy_div_c4c,       &
                                            fbcy_div_c4c1,      &
                                            fbcy_mu_c4c,        &
                                            fbcz_cc4,           &
                                            fbcz_cc41,          &
                                            fbcz_cc42,          &
                                            fbcz_div_cc4,       &
                                            fbcz_mu_cc4,        &
                                            fbcx_4pc,           &
                                            fbcx_4cp,           &
                                            fbcz_pc4,           &
                                            fbcz_cp4,           &
                                            fbcy_p4c,           &
                                            fbcy_c4p,           &
                                            fbcy_c4p1

    integer, dimension(3) ::  ncccx, ncccy, ncccz, &
                              npccx, npccy, npccz, &
                              ncpcx, ncpcy, ncpcz, &
                              nccpx, nccpy, nccpz, &
                              nppcx, nppcy, nppcz, &
                              npcpx, npcpy, npcpz, &
                              ncppx, ncppy, ncppz

    integer  :: idir
    integer  :: i, j, k
    integer  :: nx, ny, nz
    integer  :: ibcy(2), ibcz(2)
    real(WP) :: rhsx_bulk, rhsz_bulk, tmp

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
    nppcx = dm%dppc%xsz
    nppcy = dm%dppc%ysz
    nppcz = dm%dppc%zsz
    npcpx = dm%dpcp%xsz
    npcpy = dm%dpcp%ysz
    npcpz = dm%dpcp%zsz
    ncppx = dm%dcpp%xsz
    ncppy = dm%dcpp%ysz
    ncppz = dm%dcpp%zsz
!!==============================================================================
!!  test halo
!!==============================================================================
!    call nvtxStartRange("halo")
!    qxdz_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk1
!    call Get_z_1der_C2P_3D_halo(fl%qx, qxdz_pcp_xpencil, dm, dm%iAccuracy, dm%ibcz_qx)
!    call nvtxEndRange

!    !$acc update self(qxdz_pcp_xpencil)
!    if(nrank==1) then
!      write(*,*) 'halo: qxdz_pcp_xpencil @ MPI rank', nrank
!      write(*,*) qxdz_pcp_xpencil(1:4, 1, :)
!    end if

!    call nvtxStartRange("transpose")
!    qxdz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk1
!    call compute_qxdz_pcp_zpencil(qxdz_pcp_zpencil, fl%wk2, fl%wk3)
!    qxdz_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk2
!    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))     => fl%wk3
!    call transpose_z_to_y(qxdz_pcp_zpencil, apcp_ypencil, dm%dpcp)
!    call transpose_y_to_x(apcp_ypencil, qxdz_pcp_xpencil, dm%dpcp)
!    call nvtxEndRange

!    !$acc update self(qxdz_pcp_xpencil)
!    if(nrank==1) then
!      write(*,*) 'transpose: qxdz_pcp_xpencil @ MPI rank', nrank
!      write(*,*) qxdz_pcp_xpencil(1:4, 1, :)
!    end if

!!==============================================================================
!!  x-momentum equation
!!==============================================================================
    idir = 1

    !$acc kernels default(present)
    fl%mx_rhs          = ZERO
    !$acc end kernels

!!------------------------------------------------------------------------------
!!  X-mom convection term 1/3 at (i', j, k): 
!!  conv-x-m1 = - d(gxix * qxix)/dx
!!------------------------------------------------------------------------------
    fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc1
    if (is_fbcx_velo_required) then
      if (.not. dm%is_thermo) then
        ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,4
          fbcx_4cc(i,j,k) = -dm%fbcx_qx(i,j,k) * dm%fbcx_qx(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      else
        ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,4
            fbcx_4cc(i,j,k) = -dm%fbcx_gx(i,j,k) * dm%fbcx_qx(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    else
      ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3))     => fl%wk1
    qxix_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk2
    call Get_x_midp_P2C_3D(fl%qx, qxix_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)

    if (.not. dm%is_thermo) then
      nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        accc_xpencil(i,j,k) = -qxix_ccc_xpencil(i,j,k) * qxix_ccc_xpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      gxix_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk3
      call Get_x_midp_P2C_3D(fl%gx, gxix_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_gx)
      nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        accc_xpencil(i,j,k) = -qxix_ccc_xpencil(i,j,k) * gxix_ccc_xpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
    call Get_x_1der_C2P_3D(accc_xpencil, apcc_xpencil, dm, dm%iAccuracy, mbcx_cov1, fbcx_4cc)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mx_rhs(i,j,k) = fl%mx_rhs(i,j,k) + apcc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  X-mom convection term 2/3 at (i', j, k):
!!  conv-y-m1 = - 1/r  * d(gyix * qxiy)/dy
!!------------------------------------------------------------------------------
    qxiy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk1
    call compute_qxiy_ppc_ypencil(qxiy_ppc_ypencil, fl%wk2)

    qyix_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk2
    qyix_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk3
    call Get_x_midp_C2P_3D(fl%qy, qyix_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
    if(.not. dm%is_thermo) then
      call transpose_x_to_y (qyix_ppc_xpencil, qyix_ppc_ypencil, dm%dppc)
    end if

    appc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk2

    if (.not. dm%is_thermo) then
      nx = nppcy(1); ny = nppcy(2); nz = nppcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        appc_ypencil(i,j,k) = -qyix_ppc_ypencil(i,j,k) * qxiy_ppc_ypencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      appc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3))     => fl%wk4
      gyix_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk3
      call Get_x_midp_C2P_3D(fl%gy, appc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_gy)
      call transpose_x_to_y (appc_xpencil, gyix_ppc_ypencil, dm%dppc)
      nx = nppcy(1); ny = nppcy(2); nz = nppcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        appc_ypencil(i,j,k) = -gyix_ppc_ypencil(i,j,k) * qxiy_ppc_ypencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    fbcy_p4c(1:npccy(1),1:4,1:npccy(3)) => fl%wkbc1
    if(is_fbcy_velo_required) then
      call extract_dirichlet_fbcy(fbcy_p4c, appc_ypencil, dm%dppc, dm, is_reversed = .true.)
    else
      nx = npccy(1); nz = npccy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_p4c(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk1
    call Get_y_1der_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, mbcy_cov1, fbcy_p4c)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))

    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mx_rhs(i,j,k) =  fl%mx_rhs(i,j,k) + apcc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  X-mom convection term 3/3 at (i', j, k):
!!  conv-z-m1 = - 1/r * d(gzix * qxiz)/dz
!!------------------------------------------------------------------------------
    qzix_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk1
    qzix_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk2
    call Get_x_midp_C2P_3D(fl%qz, qzix_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))     => fl%wk3
    call transpose_x_to_y(qzix_pcp_xpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_z(apcp_ypencil, qzix_pcp_zpencil, dm%dpcp)

    qxiz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk2
    call compute_qxiz_pcp_zpencil(qxiz_pcp_zpencil, fl%wk3, fl%wk4)

    apcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3))     => fl%wk3
    if (.not. dm%is_thermo) then
      nx = npcpz(1); ny = npcpz(2); nz = npcpz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        apcp_zpencil(i,j,k) = -qzix_pcp_zpencil(i,j,k) * qxiz_pcp_zpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      gzix_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk1
      apcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3))     => fl%wk4
      apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))     => fl%wk5
      call Get_x_midp_C2P_3D(fl%gz, apcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_gz)
      call transpose_x_to_y(apcp_xpencil, apcp_ypencil, dm%dpcp)
      call transpose_y_to_z(apcp_ypencil, gzix_pcp_zpencil, dm%dpcp)
      nx = npcpz(1); ny = npcpz(2); nz = npcpz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        apcp_zpencil(i,j,k) = -gzix_pcp_zpencil(i,j,k) * qxiz_pcp_zpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    fbcz_pc4(1:npcpz(1),1:npcpz(2),1:4) => fl%wkbc1
    if(is_fbcz_velo_required) then
      call extract_dirichlet_fbcz(fbcz_pc4, apcp_zpencil, dm%dpcp)
    else
      nx = npcpz(1); ny = npcpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_pc4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    apcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => fl%wk1
    call Get_z_1der_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, mbcz_cov1, fbcz_pc4)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 1, IPENCIL(3))

    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk3
    call transpose_z_to_y (apcc_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mx_rhs(i,j,k) =  fl%mx_rhs(i,j,k) + apcc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  X-mom diffusion term 1/3 at (i', j, k)
!!  diff-x-m1 = d[ 2 * mu * (qxdx - 1/3 * div)]/dx
!!------------------------------------------------------------------------------
    fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc1
    if(is_fbcx_velo_required) then
!!    compute qxdx_pcc_xpencil
      qxdx_pcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
      call Get_x_1der_P2P_3D(fl%qx, qxdx_pcc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
      if(dm%icoordinate == ICYLINDRICAL) &
      call multiple_cylindrical_rn(qxdx_pcc_xpencil, dm%dpcc, dm%rc, 1, IPENCIL(1))

!!    compute qydy_pcc_xpencil
      qyix_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk2
      appc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3))     => fl%wk3
      call Get_x_midp_C2P_3D(fl%qy, qyix_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
      call transpose_x_to_y (qyix_ppc_xpencil, appc_ypencil, dm%dppc)
      fbcy_p4c(1:npccy(1),1:4,1:npccy(3)) => fl%wkbc2
      if(is_fbcy_velo_required) then
        call extract_dirichlet_fbcy(fbcy_p4c, appc_ypencil, dm%dppc, dm, is_reversed = .true.)
      else
        nx = npccy(1); nz = npccy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_p4c(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if
      qydy_pcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
      apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3))     => fl%wk4
      call Get_y_1der_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, fbcy_p4c)
      call transpose_y_to_x(apcc_ypencil, qydy_pcc_xpencil, dm%dpcc)

!!    compute qzdz_pcc_xpencil
      qzix_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk3
      apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))     => fl%wk4
      qzix_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk5
      call Get_x_midp_C2P_3D(fl%qz, qzix_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
      call transpose_x_to_y(qzix_pcp_xpencil, apcp_ypencil, dm%dpcp)
      call transpose_y_to_z(apcp_ypencil, qzix_pcp_zpencil, dm%dpcp)
      fbcz_pc4(1:npcpz(1),1:npcpz(2),1:4) => fl%wkbc2
      if(is_fbcz_velo_required) then
        call extract_dirichlet_fbcz(fbcz_pc4, qzix_pcp_zpencil, dm%dpcp)
      else
        nx = npcpz(1); ny = npcpz(2)
        !$acc parallel loop collapse(3) default(present)
        do k=1,4; do j=1,ny; do i=1,nx
          fbcz_pc4(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if
      qzdz_pcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk3
      apcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3))     => fl%wk4
      call Get_z_1der_P2C_3D(qzix_pcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, fbcz_pc4)
      apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3))     => fl%wk5
      call transpose_z_to_y(apcc_zpencil, apcc_ypencil, dm%dpcc )
      call transpose_y_to_x(apcc_ypencil, qzdz_pcc_xpencil, dm%dpcc)

!!    compute fbcx_div_4cc, qxdx + 1/r * qydy + 1/r qzdz
      nx = npccx(1); ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        qzdz_pcc_xpencil(i,j,k) = qzdz_pcc_xpencil(i,j,k) + qydy_pcc_xpencil(i,j,k) + qxdx_pcc_xpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
      if(dm%icoordinate == ICYLINDRICAL) &
      call multiple_cylindrical_rn(qzdz_pcc_xpencil, dm%dpcc, dm%rci, 1, IPENCIL(1))
      fbcx_div_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc2
      call extract_dirichlet_fbcx(fbcx_div_4cc, qzdz_pcc_xpencil, dm%dpcc)

!!    compute fbcx_4cc
      call extract_dirichlet_fbcx(fbcx_4cc, qxdx_pcc_xpencil, dm%dpcc)

!!    compute fbcx_mu_4cc
      fbcx_mu_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc3
      if(.not.dm%is_thermo) then
        ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,4
         fbcx_mu_4cc(i,j,k) = ONE
        end do; end do; end do
        !$acc end parallel loop
      else
        fbcx_4cc1(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc4
        ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,4
          fbcx_4cc1(i,j,k) = dm%fbcx_ftp(i,j,k)%m
        end do; end do; end do
        !$acc end parallel loop
        apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
        apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk2
        call Get_x_midp_C2P_3D(fl%mVisc, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4cc1)
        call transpose_x_to_y(apcc_xpencil, apcc_ypencil, dm%dpcc)
        call extract_dirichlet_fbcx(fbcx_mu_4cc, apcc_xpencil, dm%dpcc)
      end if
      ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = TWO * fbcx_mu_4cc(i,j,k) * (fbcx_4cc(i,j,k) - ONE_THIRD * fbcx_div_4cc(i,j,k))
      end do; end do; end do
      !$acc end parallel loop
    else
      ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

!!  compute div_ccc_xpencil
    div_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3))  => fl%wk1
    call compute_div_ccc_xpencil(div_ccc_xpencil, fl%wk2, fl%wk3, fl%wk4)

!!  compute mu_ccc_xpencil
    mu_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3))   => fl%wk2
    call compute_mu_ccc_xpencil(mu_ccc_xpencil)

    qxdx_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk3
    call compute_qxdx_ccc_xpencil(qxdx_ccc_xpencil)

    accc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3))     => fl%wk1
    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accc_xpencil(i,j,k) = TWO * mu_ccc_xpencil(i,j,k) * (qxdx_ccc_xpencil(i,j,k) &
                          - ONE_THIRD * div_ccc_xpencil(i,j,k))
    end do; end do; end do
    !$acc end parallel loop

    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
    call Get_x_1der_C2P_3D(accc_xpencil, apcc_xpencil, dm, dm%iAccuracy, mbcx_tau1, fbcx_4cc) 

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mx_rhs(i,j,k) = fl%mx_rhs(i,j,k) + apcc_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!!------------------------------------------------------------------------------
!!!  X-mom diffusion term 2/3 at (i', j, k)
!!!  diff-y-m1 = 1/r  * d[muixy * (qydx + r * qxdy)]/dy
!!!------------------------------------------------------------------------------
    qxdy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk1
    call compute_qxdy_ppc_ypencil(qxdy_ppc_ypencil, fl%wk2)

    qydx_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk2
    qydx_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk3
    call Get_x_1der_C2P_3D(fl%qy, qydx_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
    call transpose_x_to_y(qydx_ppc_xpencil, qydx_ppc_ypencil, dm%dppc)

    muixy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk3
    if(.not. dm%is_thermo) then
      nx = nppcy(1); ny = nppcy(2); nz = nppcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        muixy_ppc_ypencil(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      muixy_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk4
      call compute_muixy_ppc_xpencil(muixy_ppc_xpencil, fl%wk3, fl%wk5, fl%wkbc1, fl%wkbc2, fl%wkbc3)
      call transpose_x_to_y(muixy_ppc_xpencil, muixy_ppc_ypencil, dm%dppc)
    end if

    appc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk1
    if(dm%outlet_sponge_layer(1) > MINP) then
      nx = nppcy(1); ny = nppcy(2); nz = nppcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        appc_ypencil(i,j,k) = qydx_ppc_ypencil(i,j,k) * muixy_ppc_ypencil(i,j,k) &
            + qxdy_ppc_ypencil(i,j,k) * (muixy_ppc_ypencil(i,j,k) + safe_divide(fl%rre_sponge_p(i),fl%rre))
      end do; end do; end do
      !$acc end parallel loop
    else
      nx = nppcy(1); ny = nppcy(2); nz = nppcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        appc_ypencil(i,j,k) = (qxdy_ppc_ypencil(i,j,k) + qydx_ppc_ypencil(i,j,k)) * muixy_ppc_ypencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    fbcy_p4c(1:npccy(1),1:4,1:npccy(3)) => fl%wkbc1
    if(is_fbcy_velo_required) then
      call extract_dirichlet_fbcy(fbcy_p4c, appc_ypencil, dm%dppc, dm, is_reversed = .true.)
    else
      nx = nppcy(1); nz = nppcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_p4c(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk2
!!  FIXME: The following operation is where the difference happens up to the last 2-3 digits
!!        between the CPU and GPU versions. appc_ypencil only have tiny difference (last digit)
!!        fbcy_p4c is completely the same. The difference is amplified via mesh scaling factor in
!!        subroutine Prepare_TDMA_1deri_P2C_RHS_array: dd
    call Get_y_1der_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, mbcy_tau1, fbcy_p4c)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))

    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mx_rhs(i,j,k) = fl%mx_rhs(i,j,k) + apcc_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  X-mom diffusion term 3/3 at (i', j, k)
!!  diff-z-m1 = 1/r * d[muixz * (qzdx + 1/r * qxdz)]/dz
!!------------------------------------------------------------------------------
    qxdz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk1
    call compute_qxdz_pcp_zpencil(qxdz_pcp_zpencil, fl%wk2, fl%wk3)

    qzdx_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk2
    qzdx_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk3
    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))     => fl%wk4
    call Get_x_1der_C2P_3D(fl%qz, qzdx_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    call transpose_x_to_y(qzdx_pcp_xpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_z(apcp_ypencil, qzdx_pcp_zpencil, dm%dpcp)

    muixz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk3
    call compute_muixz_pcp_zpencil(muixz_pcp_zpencil, fl%wk4, fl%wk5, fl%wkbc1)

    apcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk1
    if(dm%outlet_sponge_layer(1) > MINP) then
      nx = npcpz(1); ny = npcpz(2); nz = npcpz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        apcp_zpencil(i,j,k) = qzdx_pcp_zpencil(i,j,k) * muixz_pcp_zpencil(i,j,k) &
            + qxdz_pcp_zpencil(i,j,k) * (muixz_pcp_zpencil(i,j,k) + safe_divide(fl%rre_sponge_p(i),fl%rre))
      end do; end do; end do
      !$acc end parallel loop
    else
      nx = npcpz(1); ny = npcpz(2); nz = npcpz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        apcp_zpencil(i,j,k) = (qxdz_pcp_zpencil(i,j,k) + qzdx_pcp_zpencil(i,j,k)) * muixz_pcp_zpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    fbcz_pc4(1:npccz(1),1:npccz(2),1:4) => fl%wkbc1
    if(is_fbcz_velo_required) then
      call extract_dirichlet_fbcz(fbcz_pc4, apcp_zpencil, dm%dpcp)
    else
      nx = npcpz(1); ny = npcpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_pc4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    apcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => fl%wk2
    call Get_z_1der_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, mbcz_tau1, fbcz_pc4)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 1, IPENCIL(3))

    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk3
    call transpose_z_to_y(apcc_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x(apcc_ypencil, apcc_xpencil, dm%dpcc)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mx_rhs(i,j,k) = fl%mx_rhs(i,j,k) + apcc_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!!------------------------------------------------------------------------------
!!  pressure gradients and all other body forces
!!  X-mom pressure gradient:
!!  p-m1 = - dpdx
!!------------------------------------------------------------------------------
    mx_rhs_pfc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1

!!  TODO: not needed? to be removed
    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      mx_rhs_pfc_xpencil(i,j,k) = ZERO
    end do; end do; end do
    !$acc end parallel loop

    call Get_x_1der_C2P_3D(fl%pres, mx_rhs_pfc_xpencil, dm, dm%iAccuracy, dm%ibcx_pr, dm%fbcx_pr)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      mx_rhs_pfc_xpencil(i,j,k) = -mx_rhs_pfc_xpencil(i,j,k) * dm%sigma1p
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  X-mom gravity in x direction, x-pencil
!!------------------------------------------------------------------------------
    if(dm%is_thermo .and. (fl%igravity == idir .or. fl%igravity == -idir)) then
      fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc1
      ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = dm%fbcx_ftp(i,j,k)%d
      end do; end do; end do
      !$acc end parallel loop

      apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
      call Get_x_midp_C2P_3D(fl%dDens, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4cc)

      nx = npccx(1); ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mx_rhs_pfc_xpencil(i,j,k) =  mx_rhs_pfc_xpencil(i,j,k) + fl%fgravity(idir) * apcc_xpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

!!------------------------------------------------------------------------------
!!  X-mom Lorentz Force in x direction, x-pencil
!!------------------------------------------------------------------------------
    if(dm%is_mhd) then
      nx = npccx(1); ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mx_rhs_pfc_xpencil(i,j,k) = mx_rhs_pfc_xpencil(i,j,k) + fl%lrfx(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

!!  To build up final x-momentum rhs in fractional steps
    call Calculate_momentum_fractional_step(fl%mx_rhs0, fl%mx_rhs, mx_rhs_pfc_xpencil, dm%dpcc, dm, isub)

!!==============================================================================
!!  y-momentum equation
!!==============================================================================
    idir = 2

    !$acc kernels default(present)
    fl%my_rhs          = ZERO
    !$acc end kernels

!!------------------------------------------------------------------------------
!!  Y-mom convection term 1/4 at (i, j', k)
!!  conv-x-m2 = - d(gxiy * qyix)/dx
!!------------------------------------------------------------------------------
    qyix_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk1
    call Get_x_midp_C2P_3D(fl%qy, qyix_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)

    appc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3))     => fl%wk2
    if (.not. dm%is_thermo) then
      qxiy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk2
      call compute_qxiy_ppc_ypencil(qxiy_ppc_ypencil, fl%wk3)
      qxiy_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk3
      call transpose_y_to_x(qxiy_ppc_ypencil, qxiy_ppc_xpencil, dm%dppc)
      nx = nppcx(1); ny = nppcx(2); nz = nppcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        appc_xpencil(i,j,k) = -qxiy_ppc_xpencil(i,j,k) * qyix_ppc_xpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      gxiy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk2
      call compute_gxiy_ppc_ypencil(gxiy_ppc_ypencil, fl%wk3)
      gxiy_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk3
      call transpose_y_to_x (gxiy_ppc_ypencil, gxiy_ppc_xpencil, dm%dppc)
      nx = nppcx(1); ny = nppcx(2); nz = nppcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        appc_xpencil(i,j,k) = -gxiy_ppc_xpencil(i,j,k) * qyix_ppc_xpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    fbcx_4pc(1:4,1:nppcx(2),1:nppcx(3)) => fl%wkbc1
    if (is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4pc, appc_xpencil, dm%dppc)
    else
      ny = nppcx(2); nz = nppcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4pc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk1
    call Get_x_1der_P2C_3D(appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, mbcx_cov2, fbcx_4pc)

    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) + acpc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Y-mom convection term 2/4 at (i, j', k)
!!  conv-y-m2 = - d(gyiy * qyriy)/dy
!!------------------------------------------------------------------------------
    fbcy_c4c1(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
    fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3))  => fl%wkbc2
    nx = ncpcy(1); nz = ncpcy(3)
    if(is_fbcy_velo_required) then
      if(dm%icoordinate == ICYLINDRICAL) then
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c1(i,j,k) = dm%fbcy_qyr(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      else
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c1(i,j,k) = dm%fbcy_qy(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
      if ( .not. dm%is_thermo) then
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c(i,j,k) = -fbcy_c4c1(i,j,k) * dm%fbcy_qy(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      else
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c(i,j,k) = -fbcy_c4c1(i,j,k) * dm%fbcy_gy(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    else
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    if(dm%icoordinate == ICYLINDRICAL) then
      qyr_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))       => fl%wk1
      call compute_qyr_ypencil(qyr_ypencil, fl%wk2)
      qyriy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
      call Get_y_midp_P2C_3D(qyr_ypencil, qyriy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
      accc_ypencil1 => qyriy_ccc_ypencil
    else
      qyiy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))  => fl%wk2
      call compute_qyiy_ccc_ypencil(qyiy_ccc_ypencil, fl%wk1)
      accc_ypencil1     => qyiy_ccc_ypencil
    end if

    if (.not. dm%is_thermo) then
      qyiy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))  => fl%wk1
      call compute_qyiy_ccc_ypencil(qyiy_ccc_ypencil, fl%wk3)
      accc_ypencil => qyiy_ccc_ypencil
    else
      gyiy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))  => fl%wk1
      call compute_gyiy_ccc_ypencil(gyiy_ccc_ypencil, fl%wk3)
      accc_ypencil => gyiy_ccc_ypencil
    end if

    nx = ncccy(1); ny = ncccy(2); nz = ncccy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accc_ypencil(i,j,k) = -accc_ypencil1(i,j,k) * accc_ypencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
    call Get_y_1der_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcy_cov2, fbcy_c4c)

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk1
    call transpose_y_to_x (acpc_ypencil, acpc_xpencil, dm%dcpc)

    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) + acpc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Y-mom convection term 3/4 at (i, j', k)
!!  conv-z-m2 = - d(gziy * qyriz)/dz
!!------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      qyr_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))       => fl%wk1
      call compute_qyr_ypencil(qyr_ypencil, fl%wk2)
      acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3))      => fl%wk2
      call transpose_y_to_z(qyr_ypencil, acpc_zpencil, dm%dcpc)
      qyriz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk3
      call Get_z_midp_C2P_3D(acpc_zpencil, qyriz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qyr)
      acpp_zpencil1     => qyriz_cpp_zpencil
    else
      qyiz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))  => fl%wk3
      call compute_qyiz_cpp_zpencil(qyiz_cpp_zpencil, fl%wk1, fl%wk2)
      acpp_zpencil1     => qyiz_cpp_zpencil
    end if

    if (.not. dm%is_thermo) then
      qziy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk2
      call compute_qziy_cpp_zpencil(qziy_cpp_zpencil, fl%wk1, fl%wk4)
      acpp_zpencil     => qziy_cpp_zpencil
    else
      gziy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk2
      call compute_gziy_cpp_zpencil(gziy_cpp_zpencil, fl%wk1, fl%wk4)
      acpp_zpencil     => gziy_cpp_zpencil
    end if

    nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      acpp_zpencil(i,j,k) = -acpp_zpencil1(i,j,k) * acpp_zpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    fbcz_cp4(1:ncppz(1),1:ncppz(2),1:4)  => fl%wkbc1
    if(is_fbcz_velo_required) then
      call extract_dirichlet_fbcz(fbcz_cp4, acpp_zpencil, dm%dcpp)
    else
      nx = ncppz(1); ny = ncppz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_cp4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => fl%wk1
    call Get_z_1der_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, mbcz_cov2, fbcz_cp4)

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk2
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
    call transpose_z_to_y(acpc_zpencil, acpc_ypencil, dm%dcpc)
    call transpose_y_to_x(acpc_ypencil, acpc_xpencil, dm%dcpc)

    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) + acpc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Y-mom convection term 4/4 at (i, j', k)
!!  conv-r-m2 = (gziz * qziz)^y
!!------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
      if(is_fbcy_velo_required) then
        qziy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk1
        call compute_qziy_cpp_zpencil(qziy_cpp_zpencil, fl%wk2, fl%wk3)
        acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3))     => fl%wk2
        call Get_z_midp_P2C_3D(qziy_cpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
        acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))     => fl%wk3
        call transpose_z_to_y(acpc_zpencil, acpc_ypencil, dm%dcpc)
        fbcy_c4c1(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc2
        call extract_dirichlet_fbcy(fbcy_c4c1, acpc_ypencil, dm%dcpc, dm, is_reversed = .true.)
        if (.not. dm%is_thermo) then
          nx = ncpcy(1); nz = ncpcy(3)
          !$acc parallel loop collapse(3) default(present)
          do k=1,nz; do j=1,4; do i=1,nx
            fbcy_c4c(i,j,k) = fbcy_c4c1(i,j,k)
          end do; end do; end do
          !$acc end parallel loop
        else
          gziy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk1
          call compute_gziy_cpp_zpencil(gziy_cpp_zpencil, fl%wk2, fl%wk3)
          acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3))     => fl%wk2
          call Get_z_midp_P2C_3D(gziy_cpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
          acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))     => fl%wk3
          call transpose_z_to_y(acpc_zpencil, acpc_ypencil, dm%dcpc)
          call extract_dirichlet_fbcy(fbcy_c4c, acpc_ypencil, dm%dcpc, dm, is_reversed = .true.)
        end if
        nx = ncpcy(1); nz = ncpcy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c(i,j,k) = fbcy_c4c(i,j,k) * fbcy_c4c1(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncpcy(1); nz = ncpcy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if

      qziz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
      call compute_qziz_ccc_zpencil(qziz_ccc_zpencil, fl%wk2, fl%wk3)
      qziz_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
      call transpose_z_to_y(qziz_ccc_zpencil, qziz_ccc_ypencil, dm%dccc)
      if (.not. dm%is_thermo) then
        accc_ypencil1 => qziz_ccc_ypencil
      else
        gziz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
        call compute_gziz_ccc_zpencil(gziz_ccc_zpencil, fl%wk3, fl%wk4)
        gziz_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk3
        call transpose_z_to_y(gziz_ccc_zpencil, gziz_ccc_ypencil, dm%dccc)
        accc_ypencil1 => gziz_ccc_ypencil
      end if

      accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk1
      nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        accc_ypencil(i,j,k) = accc_ypencil1(i,j,k) * qziz_ccc_ypencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop

      acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
      call Get_y_midp_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcr_cov2, fbcy_c4c)
      call axis_estimating_radial_xpx(acpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(2), is_reversed = .true.)

      acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk1
      call transpose_y_to_x(acpc_ypencil, acpc_xpencil, dm%dcpc)

      nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) + acpc_xpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop

    end if

!!------------------------------------------------------------------------------
!!  Y-mom diffusion term 1/4 at (i, j', k)
!!  diff-x-m2 = d[muixy * (qydx + r * qxdy)]/dx
!!------------------------------------------------------------------------------
    muixy_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk1
    call compute_muixy_ppc_xpencil(muixy_ppc_xpencil, fl%wk2, fl%wk3, fl%wkbc1, fl%wkbc2, fl%wkbc3)

    qxdy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk2
    call compute_qxdy_ppc_ypencil(qxdy_ppc_ypencil, fl%wk3)
    qxdy_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk3
    call transpose_y_to_x(qxdy_ppc_ypencil, qxdy_ppc_xpencil, dm%dppc)

    qydx_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk2
    call Get_x_1der_C2P_3D(fl%qy, qydx_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)

    appc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3))     => fl%wk1
    nx = nppcx(1); ny = nppcx(2); nz = nppcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      appc_xpencil(i,j,k) = (qxdy_ppc_xpencil(i,j,k) + qydx_ppc_xpencil(i,j,k)) * muixy_ppc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    fbcx_4pc(1:4,1:nppcx(2),1:nppcx(3)) => fl%wkbc1
    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4pc, appc_xpencil, dm%dppc)
    else
      ny = nppcx(2); nz = nppcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4pc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk2
    call Get_x_1der_P2C_3D(appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, mbcx_tau2, fbcx_4pc)

    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) + acpc_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Y-mom diffusion term 2/4 at (i, j', k)
!!  diff-y-m2 = d[r * 2 * mu * (qyrdy - 1/3 * div)]/dy
!!------------------------------------------------------------------------------
    fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
    if(is_fbcy_velo_required) then
      if(dm%icoordinate == ICYLINDRICAL) then
        qyr_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))  => fl%wk1
        call compute_qyr_ypencil(qyr_ypencil, fl%wk2)
        accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
        call Get_y_midp_P2C_3D(qyr_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
        acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk1
        call Get_y_1der_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
        call extract_dirichlet_fbcy(fbcy_c4c, acpc_ypencil, dm%dcpc, dm, is_reversed = .true.)
      else
        qydy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk1
        call compute_qydy_cpc_ypencil(qydy_cpc_ypencil, fl%wk2, fl%wk3)
        call extract_dirichlet_fbcy(fbcy_c4c, qydy_cpc_ypencil, dm%dcpc, dm, is_reversed = .true.)
      end if

      fbcy_div_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc2
!      call compute_fbcy_div_c4c(fbcy_div_c4c, fl%wk1, fl%wk2, fl%wk3, fl%wkbc3, fl%wkbc4)
      call compute_fbcy_div_c4c_1(fbcy_div_c4c, fl%wk1, fl%wk2, fl%wk3, fl%wkbc3, fl%wkbc4, fl%wkbc5)

      fbcy_mu_c4c(1:ncpcy(1),1:4,1:ncpcy(3))  => fl%wkbc3
      call compute_fbcy_mu_c4c(fbcy_mu_c4c, fl%wk1, fl%wk2, fl%wkbc4)

      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = (fbcy_c4c(i,j,k) - ONE_THIRD * fbcy_div_c4c(i,j,k)) * TWO * fbcy_mu_c4c(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
      if(dm%icoordinate == ICYLINDRICAL) then
        call multiple_cylindrical_rn_x4x(fbcy_c4c, dm%dcpc, dm%rp, 1, IPENCIL(2))
      end if
    else
      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    div_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk2
    call compute_div_ccc_xpencil(div_ccc_xpencil, fl%wk1, fl%wk3, fl%wk4)
    div_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk1
    call transpose_x_to_y (div_ccc_xpencil, div_ccc_ypencil, dm%dccc)

    if(dm%icoordinate == ICYLINDRICAL) then
      qyr_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))       => fl%wk2
      call compute_qyr_ypencil(qyr_ypencil, fl%wk3)
      qyrdy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk3
      call Get_y_1der_P2C_3D(qyr_ypencil, qyrdy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
      accc_ypencil1     => qyrdy_ccc_ypencil
    else
      qydy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))  => fl%wk3
      call compute_qydy_ccc_ypencil(qydy_ccc_ypencil, fl%wk2)
      accc_ypencil1     => qydy_ccc_ypencil
    end if

    mu_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
    call compute_mu_ccc_ypencil(mu_ccc_ypencil)

    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk1
    if(dm%outlet_sponge_layer(1) > MINP) then
      nx = ncccy(1); ny = ncccy(2); nz = ncccy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        accc_ypencil(i,j,k) = -ONE_THIRD * div_ccc_ypencil(i,j,k) * TWO * mu_ccc_ypencil(i,j,k) &
            + accc_ypencil1(i,j,k) * TWO * (mu_ccc_ypencil(i,j,k) + safe_divide(fl%rre_sponge_c(i),fl%rre))
      end do; end do; end do
      !$acc end parallel loop
    else
      nx = ncccy(1); ny = ncccy(2); nz = ncccy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        accc_ypencil(i,j,k) = (accc_ypencil1(i,j,k) - ONE_THIRD * div_ccc_ypencil(i,j,k)) * TWO * mu_ccc_ypencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accc_ypencil, dm%dccc, dm%rc, 1, IPENCIL(2))

    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
    call Get_y_1der_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcy_tau2, fbcy_c4c)

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk1
    call transpose_y_to_x(acpc_ypencil, acpc_xpencil, dm%dcpc)

    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) + acpc_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!  Note: calculation of qydy_cpc_ypencil, qydy_ccc_ypencil, div_ccc_xpencil have a machine
!!  level difference from those of the CPU version - they all have P2C and C2P operations in 
!!  y-direction involved. (disappear for cylindrical)

!!------------------------------------------------------------------------------
!!  Y-mom diffusion term 3/4 at (i, j', k)
!!  diff-z-m2 = d[muiyz * (qzdy + qyr2dz - qzriy]/dz
!!------------------------------------------------------------------------------
    muiyz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))  => fl%wk1
    call compute_muiyz_cpp_zpencil(muiyz_cpp_zpencil, fl%wk2, fl%wk3, fl%wkbc1)

    acpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))       => fl%wk1
    if(dm%icoordinate == ICYLINDRICAL) then
      qyr2dz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk2
      call compute_qyr2dz_cpp_zpencil(qyr2dz_cpp_zpencil, fl%wk3, fl%wk4)
      qzriy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))  => fl%wk3
      call compute_qzriy_cpp_ypencil(qzriy_cpp_ypencil, fl%wk4, fl%wk5)
      qzriy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))  => fl%wk4
      call transpose_y_to_z(qzriy_cpp_ypencil, qzriy_cpp_zpencil, dm%dcpp)
      qzdy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))   => fl%wk3
      call compute_qzdy_cpp_ypencil(qzdy_cpp_ypencil, fl%wk5)
      qzdy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))   => fl%wk5
      call transpose_y_to_z(qzdy_cpp_ypencil, qzdy_cpp_zpencil, dm%dcpp)
      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_zpencil(i,j,k) = (qzdy_cpp_zpencil(i,j,k) + qyr2dz_cpp_zpencil(i,j,k) &
              - qzriy_cpp_zpencil(i,j,k)) * (muiyz_cpp_zpencil(i,j,k) + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_zpencil(i,j,k) = (qzdy_cpp_zpencil(i,j,k) + qyr2dz_cpp_zpencil(i,j,k) &
                                 - qzriy_cpp_zpencil(i,j,k)) * muiyz_cpp_zpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    else
      qydz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk2
      call compute_qydz_cpp_zpencil(qydz_cpp_zpencil, fl%wk3, fl%wk4)
      qzdy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))   => fl%wk3
      call compute_qzdy_cpp_ypencil(qzdy_cpp_ypencil, fl%wk4)
      qzdy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))   => fl%wk4
      call transpose_y_to_z(qzdy_cpp_ypencil, qzdy_cpp_zpencil, dm%dcpp)
      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_zpencil(i,j,k) = (qzdy_cpp_zpencil(i,j,k) + qydz_cpp_zpencil(i,j,k)) &
                              * (muiyz_cpp_zpencil(i,j,k) + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_zpencil(i,j,k) = (qzdy_cpp_zpencil(i,j,k) + qydz_cpp_zpencil(i,j,k)) * muiyz_cpp_zpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    end if

    fbcz_cp4(1:ncppz(1),1:ncppz(2),1:4) => fl%wkbc1
    if(is_fbcz_velo_required) then
      call extract_dirichlet_fbcz(fbcz_cp4, acpp_zpencil, dm%dcpp)
    else
      nx = ncppz(1); ny = ncppz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_cp4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => fl%wk2
    call Get_z_1der_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, mbcz_tau2, fbcz_cp4)

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk1
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
    call transpose_z_to_y(acpc_zpencil, acpc_ypencil, dm%dcpc)
    call transpose_y_to_x(acpc_ypencil, acpc_xpencil, dm%dcpc)

    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) + acpc_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!   Note: last digit different appears at muiyz_cpp_zpencil for cylindrical with thermal

!!------------------------------------------------------------------------------
!!  Y-mom diffusion term 4/4 at (i, j', k)
!!  diff-r-m2 = - [2 mu * (1/r * qzdz - 1/3 * div + 1/r * qyriy)]^y
!!------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
      if(is_fbcy_velo_required) then
        qzrdz_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk1
        call compute_qzrdz_cpc_ypencil(qzrdz_cpc_ypencil, fl%wk2, fl%wk3, fl%wk4)
        qyr2_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))      => fl%wk2
        call compute_qyr2_ypencil(qyr2_ypencil, fl%wk3)
        fbcy_c4c1(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc2
        call extract_dirichlet_fbcy(fbcy_c4c, qzrdz_cpc_ypencil, dm%dcpc, dm, is_reversed = .false.)
        call extract_dirichlet_fbcy(fbcy_c4c1, qyr2_ypencil, dm%dcpc, dm, is_reversed = .true.)
        nx = ncpcy(2); nz = ncpcy(3)
        fbcy_div_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc3
        call compute_fbcy_div_c4c(fbcy_div_c4c, fl%wk1, fl%wk2, fl%wk3, fl%wkbc4, fl%wkbc5)
        fbcy_mu_c4c(1:ncpcy(1),1:4,1:ncpcy(3))  => fl%wkbc4
        call compute_fbcy_mu_c4c(fbcy_mu_c4c, fl%wk1, fl%wk2, fl%wkbc5)
        nx = ncpcy(2); nz = ncpcy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c(i,j,k) = (fbcy_c4c(i,j,k) + fbcy_c4c1(i,j,k) - ONE_THIRD * fbcy_div_c4c(i,j,k)) &
                             * TWO * fbcy_mu_c4c(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncpcy(2); nz = ncpcy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4c(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if

      accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))     => fl%wk1
      div_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3))  => fl%wk2
      call compute_div_ccc_xpencil(div_ccc_xpencil, fl%wk1, fl%wk3, fl%wk4)
      div_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))  => fl%wk1
      call transpose_x_to_y(div_ccc_xpencil, div_ccc_ypencil, dm%dccc)

      qzdz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk3
      call compute_qzdz_ccc_zpencil(qzdz_ccc_zpencil, fl%wk2, fl%wk4)
      qzdz_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
      call transpose_z_to_y(qzdz_ccc_zpencil, qzdz_ccc_ypencil, dm%dccc)
      call multiple_cylindrical_rn(qzdz_ccc_ypencil, dm%dccc, dm%rci, 1, IPENCIL(2))

      qyr_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))       => fl%wk4
      call compute_qyr_ypencil(qyr_ypencil, fl%wk3)
      qyriy_ccc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
      call Get_y_midp_P2C_3D(qyr_ypencil, qyriy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
      call multiple_cylindrical_rn(qyriy_ccc_ypencil, dm%dccc, dm%rci, 1, IPENCIL(2))

      mu_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))    => fl%wk4
      call compute_mu_ccc_ypencil(mu_ccc_ypencil)

      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncccy(1); ny = ncccy(2); nz = ncccy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accc_ypencil(i,j,k) = -ONE_THIRD * div_ccc_ypencil(i,j,k) * TWO * mu_ccc_ypencil(i,j,k)   &
              + (qzdz_ccc_ypencil(i,j,k) + qyriy_ccc_ypencil(i,j,k)) * TWO * (mu_ccc_ypencil(i,j,k) &
              + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncccy(1); ny = ncccy(2); nz = ncccy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accc_ypencil(i,j,k) = (qzdz_ccc_ypencil(i,j,k) + qyriy_ccc_ypencil(i,j,k) - ONE_THIRD * div_ccc_ypencil(i,j,k)) &
                              * TWO * mu_ccc_ypencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if

      acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
      call Get_y_midp_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcr_tau2, fbcy_c4c)

      acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk1
      call transpose_y_to_x(acpc_ypencil, acpc_xpencil, dm%dcpc)

      nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%my_rhs(i,j,k) = fl%my_rhs(i,j,k) - acpc_xpencil(i,j,k) * fl%rre
      end do; end do; end do
      !$acc end parallel loop

    end if

!!------------------------------------------------------------------------------
!!  pressure gradients and all other body forces
!!  Y-mom pressure gradient in y direction, Y-pencil, d(sigma_1 p)
!!  p-m2 = - r * dpdy
!!------------------------------------------------------------------------------
    my_rhs_pfc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk1

    pres_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))       => fl%wk2
    call transpose_x_to_y(fl%pres, pres_ypencil, dm%dccc)
    call Get_y_1der_C2P_3D(pres_ypencil, my_rhs_pfc_ypencil, dm, dm%iAccuracy, dm%ibcy_pr, dm%fbcy_pr)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(my_rhs_pfc_ypencil, dm%dcpc, dm%rp, 1, IPENCIL(2))

    nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      my_rhs_pfc_ypencil(i,j,k) = -my_rhs_pfc_ypencil(i,j,k) * dm%sigma1p
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!  Y-mom gravity in y direction ( not r direction ), y-pencil
!!------------------------------------------------------------------------------
    if(dm%is_thermo) then
      acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))   => fl%wk2
      !$acc kernels default(present)
      acpc_ypencil(:,:,:) = ZERO
      !$acc end kernels

      if(dm%icoordinate == ICYLINDRICAL) then
        accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk3
        if(fl%igravity == idir .or. fl%igravity == -idir) then
          call gravity_decomposition_to_rz(fl%dDens, fl%igravity, fl%fgravity(idir), &
                          acpc_ypencil, accp_zpencil, fl%wk4, fl%wk5, fl%wkbc1, dm)
        else if(fl%igravity == idir+1 .or. fl%igravity == -idir-1) then
          call gravity_decomposition_to_rz(fl%dDens, fl%igravity, fl%fgravity(idir+1), &
                            acpc_ypencil, accp_zpencil, fl%wk4, fl%wk5, fl%wkbc1, dm)
        end if
      else
        if(fl%igravity == idir .or. fl%igravity == -idir) then
          fbcy_c4c (1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
          nx = ncpcy(1); nz = ncpcy(3)
          !$acc parallel loop collapse(3) default(present)
          do k=1,nz; do j=1,4; do i=1,nx
            fbcy_c4c(i,j,k) = dm%fbcy_ftp(i,j,k)%d
          end do; end do; end do
          !$acc end parallel loop
          dDens_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk3
          call transpose_x_to_y (fl%dDens, dDens_ypencil, dm%dccc)
          call Get_y_midp_C2P_3D(dDens_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)

          nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
          !$acc parallel loop collapse(3) default(present)
          do k=1,nz; do j=1,ny; do i=1,nx
            acpc_ypencil(i,j,k) = fl%fgravity(idir) * acpc_ypencil(i,j,k)
          end do; end do; end do
          !$acc end parallel loop
        end if
      end if

      nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        my_rhs_pfc_ypencil(i,j,k) =  my_rhs_pfc_ypencil(i,j,k) + acpc_ypencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    my_rhs_pfc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk2
    call transpose_y_to_x(my_rhs_pfc_ypencil, my_rhs_pfc_xpencil, dm%dcpc)

!!------------------------------------------------------------------------------
!!  Y-mom Lorentz Force in y direction, x-pencil
!!------------------------------------------------------------------------------
    if(dm%is_mhd) then
      nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        my_rhs_pfc_xpencil(i,j,k) = my_rhs_pfc_xpencil(i,j,k) + fl%lrfy(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

!!  To build up final y-momentum rhs in fractional steps
    call Calculate_momentum_fractional_step(fl%my_rhs0, fl%my_rhs, my_rhs_pfc_xpencil, dm%dcpc, dm, isub)

!!==============================================================================
!!  z-momentum equation
!!==============================================================================
    idir = 3

    !$acc kernels default(present)
    fl%mz_rhs          = ZERO
    !$acc end kernels

!!------------------------------------------------------------------------------
!!  Z-mom convection term 1/3 at (i, j, k')
!!  conv-x-m3 = - d(gxiz * qzix)/dx
    if (.not. dm%is_thermo) then
      qxiz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk1
      call compute_qxiz_pcp_zpencil(qxiz_pcp_zpencil, fl%wk2, fl%wk3)
      qxiz_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk2
      apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))     => fl%wk3
      call transpose_z_to_y(qxiz_pcp_zpencil, apcp_ypencil, dm%dpcp)
      call transpose_y_to_x(apcp_ypencil, qxiz_pcp_xpencil, dm%dpcp)
      apcp_xpencil     => qxiz_pcp_xpencil
    else
      gxiz_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk2
      call compute_gxiz_pcp_xpencil(gxiz_pcp_xpencil, fl%wk1, fl%wk3)
      apcp_xpencil     => gxiz_pcp_xpencil
    end if

    qzix_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3))   => fl%wk1
    call Get_x_midp_C2P_3D(fl%qz, qzix_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    nx = npcpx(1); ny = npcpx(2); nz = npcpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      apcp_xpencil(i,j,k) = - apcp_xpencil(i,j,k) * qzix_pcp_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    fbcx_4cp(1:4,1:npcpx(2),1:npcpx(3)) => fl%wkbc1
    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4cp, apcp_xpencil, dm%dpcp)
    else
      ny = npcpx(2); nz = npcpx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cp(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk1
    call Get_x_1der_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, mbcx_cov3, fbcx_4cp)

    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mz_rhs(i,j,k) = fl%mz_rhs(i,j,k) + accp_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Z-mom convection term 2/3 at (i, j, k')
!!  conv-y-m3 = - 1/r^2 * d( r * qziy * gyiz)/dy
!!------------------------------------------------------------------------------
    qziy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))   => fl%wk1
    call compute_qziy_cpp_ypencil(qziy_cpp_ypencil, fl%wk2)
    if (.not. dm%is_thermo) then
      qyiz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk3
      call compute_qyiz_cpp_zpencil(qyiz_cpp_zpencil, fl%wk2, fl%wk4)
      qyiz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk2
      call transpose_z_to_y(qyiz_cpp_zpencil, qyiz_cpp_ypencil, dm%dcpp)
      acpp_ypencil     => qyiz_cpp_ypencil
    else
      gyiz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk2
      call compute_gyiz_cpp_ypencil(gyiz_cpp_ypencil, fl%wk3, fl%wk4)
      acpp_ypencil     => gyiz_cpp_ypencil
    end if

    nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      acpp_ypencil(i,j,k) = -acpp_ypencil(i,j,k) * qziy_cpp_ypencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(acpp_ypencil, dm%dcpp, dm%rp, 1, IPENCIL(2)) 

    fbcy_c4p(1:ncppy(1),1:4,1:ncppy(3)) => fl%wkbc1
    if(is_fbcy_velo_required) then
      call extract_dirichlet_fbcy(fbcy_c4p, acpp_ypencil, dm%dcpp, dm, is_reversed = .true.)
    else
      nx = ncppy(1); nz = ncppy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4p(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk1
    call Get_y_1der_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcy_cov3, fbcy_c4p)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accp_ypencil, dm%dccp, dm%rci, 2, IPENCIL(2)) 

    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk2
    call transpose_y_to_x(accp_ypencil, accp_xpencil, dm%dccp)

    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mz_rhs(i,j,k) = fl%mz_rhs(i,j,k) + accp_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!   FIXME: the last digit difference between the CPU and GPU versions
!!          (for cylindrical the difference disappears)

!!------------------------------------------------------------------------------
!!  Z-mom convection term 3/3 at (i, j, k')
!!  conv-z-m3 = 1/r  * d(gziz * qziz)/dz
!!------------------------------------------------------------------------------
    fbcz_cc4(1:nccpz(1),1:nccpz(2),1:4) => fl%wkbc1
    if(is_fbcz_velo_required) then
      if (.not. dm%is_thermo) then
        nx = nccpz(1); ny = nccpz(2)
        !$acc parallel loop collapse(3) default(present)
        do k=1,4; do j=1,ny; do i=1,nx
          fbcz_cc4(i,j,k) = dm%fbcz_qz(i,j,k) * dm%fbcz_qz(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = nccpz(1); ny = nccpz(2)
        !$acc parallel loop collapse(3) default(present)
        do k=1,4; do j=1,ny; do i=1,nx
          fbcz_cc4(i,j,k) = dm%fbcz_gz(i,j,k) * dm%fbcz_qz(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    else
      nx = nccpz(1); ny = nccpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_cc4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))      => fl%wk1
    qziz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))  => fl%wk1
    call compute_qziz_ccc_zpencil(qziz_ccc_zpencil, fl%wk2, fl%wk3)
    if (.not. dm%is_thermo) then
      accc_zpencil1 => qziz_ccc_zpencil
    else
      gziz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk2
      call compute_gziz_ccc_zpencil(gziz_ccc_zpencil, fl%wk3, fl%wk4)
      accc_zpencil1 => gziz_ccc_zpencil
    end if
    nx = ncccz(1); ny = ncccz(2); nz = ncccz(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accc_zpencil(i,j,k) = -accc_zpencil1(i,j,k) * qziz_ccc_zpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk2
    call Get_z_1der_C2P_3D(accc_zpencil, accp_zpencil, dm, dm%iAccuracy, mbcz_cov3, fbcz_cc4)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accp_zpencil, dm%dccp, dm%rci, 1, IPENCIL(3))

    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk1
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk3
    call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x(accp_ypencil, accp_xpencil, dm%dccp)

    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mz_rhs(i,j,k) = fl%mz_rhs(i,j,k) + accp_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Z-mom diffusion term 1/4  at (i, j, k')
!!  diff-x-m3 = d[muixz * (qzdx + 1/r * qxdz)]/dx
!!------------------------------------------------------------------------------
    muixz_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk1
    muixz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk2
    call compute_muixz_pcp_zpencil(muixz_pcp_zpencil, fl%wk1, fl%wk3, fl%wkbc1)
    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))      => fl%wk3
    call transpose_z_to_y(muixz_pcp_zpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_x(apcp_ypencil, muixz_pcp_xpencil, dm%dpcp)

    qxdz_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3))  => fl%wk2
    qxdz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3))  => fl%wk3
    call compute_qxdz_pcp_zpencil(qxdz_pcp_zpencil, fl%wk2, fl%wk4)

    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3))      => fl%wk4
    call transpose_z_to_y(qxdz_pcp_zpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_x(apcp_ypencil, qxdz_pcp_xpencil, dm%dpcp)

    qzdx_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk3
    call Get_x_1der_C2P_3D(fl%qz, qzdx_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)

    apcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3))     => fl%wk1
    nx = npcpx(1); ny = npcpx(2); nz = npcpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      apcp_xpencil(i,j,k) = (qzdx_pcp_xpencil(i,j,k) + qxdz_pcp_xpencil(i,j,k)) * muixz_pcp_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    fbcx_4cp(1:4,1:npcpx(2),1:npcpx(3)) => fl%wkbc1
    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4cp, apcp_xpencil, dm%dpcp)
    else
      ny = npcpx(2); nz = npcpx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cp(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk2
    call Get_x_1der_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, mbcx_tau3, fbcx_4cp)

    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mz_rhs(i,j,k) = fl%mz_rhs(i,j,k) + accp_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Z-mom diffusion term 2/4 at (i, j, k')
!!  diff-y-m3 = 1/r * d[muiyz * (r * qzdy + qyrdz - qziy]/dy
!!------------------------------------------------------------------------------
    muiyz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk1
    muiyz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk2
    call compute_muiyz_cpp_zpencil(muiyz_cpp_zpencil, fl%wk1, fl%wk3, fl%wkbc1)
    call transpose_z_to_y(muiyz_cpp_zpencil, muiyz_cpp_ypencil, dm%dcpp)

    qzdy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))  => fl%wk2
    call compute_qzdy_cpp_ypencil(qzdy_cpp_ypencil, fl%wk3)

    acpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))      => fl%wk1
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(qzdy_cpp_ypencil, dm%dcpp, dm%rp, 1, IPENCIL(2))
      qyrdz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk3
      call compute_qyrdz_cpp_ypencil(qyrdz_cpp_ypencil, fl%wk4, fl%wk5)
      qziy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))  => fl%wk4
      call compute_qziy_cpp_ypencil(qziy_cpp_ypencil, fl%wk5)
      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_ypencil(i,j,k) = (qzdy_cpp_ypencil(i,j,k) + qyrdz_cpp_ypencil(i,j,k) - qziy_cpp_ypencil(i,j,k)) &
                              * (muiyz_cpp_ypencil(i,j,k) + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_ypencil(i,j,k) = (qzdy_cpp_ypencil(i,j,k) + qyrdz_cpp_ypencil(i,j,k) - qziy_cpp_ypencil(i,j,k)) &
                              * muiyz_cpp_ypencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    else
      qydz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk3
      qydz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk4
      call compute_qydz_cpp_zpencil(qydz_cpp_zpencil, fl%wk3, fl%wk5)
      call transpose_z_to_y(qydz_cpp_zpencil, qydz_cpp_ypencil, dm%dcpp)
      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_ypencil(i,j,k) = (qzdy_cpp_ypencil(i,j,k) + qydz_cpp_ypencil(i,j,k)) * (muiyz_cpp_ypencil(i,j,k) &
                              + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_ypencil(i,j,k) = (qzdy_cpp_ypencil(i,j,k) + qydz_cpp_ypencil(i,j,k)) * muiyz_cpp_ypencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    end if

    fbcy_c4p(1:ncppy(1),1:4,1:ncppy(3)) => fl%wkbc1
    if(is_fbcy_velo_required) then
      call extract_dirichlet_fbcy(fbcy_c4p, acpp_ypencil, dm%dcpp, dm, is_reversed = .true.)
    else
      nx = ncppy(1); nz = ncppy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4p(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk2
    call Get_y_1der_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcy_tau3, fbcy_c4p)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accp_ypencil, dm%dccp, dm%rci, 1, IPENCIL(2))

    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk1
    call transpose_y_to_x(accp_ypencil, accp_xpencil, dm%dccp)

    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mz_rhs(i,j,k) = fl%mz_rhs(i,j,k) + accp_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!   FIXME: check the cause of last digit difference between the CPU and GPU versions
!!   difference from CPU at acpp_ypencil is only limited to the last digit
!!   difference from CPU at accp_ypencil has propagates to the second last digit 
!!   Note: disappear for cylindrical
!!         reappear at muiyz_cpp_ypencil for cylindrical with thermal

!!------------------------------------------------------------------------------
!!  Z-mom diffusion term 3/4 at (i, j, k')
!!  diff-z-m3 = 1/r2 * d[2 * mu * (qzdz - r/3 * div + qyriy)]/dz
!!------------------------------------------------------------------------------
    if(is_fbcz_velo_required) then
      fbcz_cc4(1:nccpz(1),1:nccpz(2),1:4)     => fl%wkbc1
      fbcz_div_cc4(1:nccpz(1),1:nccpz(2),1:4) => fl%wkbc2
      call compute_fbcz_div_cc4(fbcz_div_cc4, fl%wk1, fl%wk2, fl%wk3, fl%wk4, fl%wk5, fl%wkbc1)

      qziz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
      call compute_qziz_ccc_zpencil(qziz_ccc_zpencil, fl%wk2, fl%wk3)
      qzdz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk2
      call Get_z_1der_C2P_3D(qziz_ccc_zpencil, qzdz_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)
      call extract_dirichlet_fbcz(fbcz_cc4, qzdz_ccp_zpencil, dm%dccp)

      fbcz_mu_cc4(1:nccpz(1),1:nccpz(2),1:4) => fl%wkbc3
      call compute_fbcz_mu_cc4(fbcz_mu_cc4, fl%wk1, fl%wk2, fl%wkbc4)

      if(dm%icoordinate == ICYLINDRICAL) then
        call multiple_cylindrical_rn_xx4(fbcz_div_cc4, dm%dccp, dm%rc, 1, IPENCIL(3))
        if(dm%outlet_sponge_layer(1) > MINP) then
          nx = nccpz(1); ny = nccpz(2)
          !$acc parallel loop collapse(3) default(present)
          do k=1,4; do j=1,ny; do i=1,nx
            fbcz_cc4(i,j,k)  = -ONE_THIRD * fbcz_div_cc4(i,j,k) * TWO * fbcz_mu_cc4(i,j,k) &
                             + (fbcz_cc4(i,j,k) + dm%fbcz_qyr(i,j,k)) * TWO * (fbcz_mu_cc4(i,j,k) &
                             + safe_divide(fl%rre_sponge_c(i),fl%rre))
          end do; end do; end do
          !$acc end parallel loop
        else
          nx = nccpz(1); ny = nccpz(2)
          !$acc parallel loop collapse(3) default(present)
          do k=1,4; do j=1,ny; do i=1,nx
            fbcz_cc4(i,j,k)  = (fbcz_cc4(i,j,k) - ONE_THIRD * fbcz_div_cc4(i,j,k) + dm%fbcz_qyr(i,j,k)) &
                             * TWO * fbcz_mu_cc4(i,j,k)
          end do; end do; end do
          !$acc end parallel loop
        end if
      else
        if(dm%outlet_sponge_layer(1) > MINP) then
          nx = nccpz(1); ny = nccpz(2)
          !$acc parallel loop collapse(3) default(present)
          do k=1,4; do j=1,ny; do i=1,nx
            fbcz_cc4(i,j,k) = -ONE_THIRD * fbcz_div_cc4(i,j,k) * TWO * fbcz_mu_cc4(i,j,k) &
                            + (fbcz_cc4(i,j,k)) * TWO * (fbcz_mu_cc4(i,j,k) &
                            + safe_divide(fl%rre_sponge_c(i),fl%rre))
          end do; end do; end do
          !$acc end parallel loop
        else
          nx = nccpz(1); ny = nccpz(2)
          !$acc parallel loop collapse(3) default(present)
          do k=1,4; do j=1,ny; do i=1,nx
            fbcz_cc4(i,j,k) = (fbcz_cc4(i,j,k) - ONE_THIRD * fbcz_div_cc4(i,j,k)) * TWO * fbcz_mu_cc4(i,j,k)
          end do; end do; end do
          !$acc end parallel loop
        end if
      end if
    else
      nx = nccpz(1); ny = nccpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_cc4(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    div_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk1
    div_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk2
    call compute_div_ccc_xpencil(div_ccc_xpencil, fl%wk1, fl%wk3, fl%wk4)

    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))    => fl%wk3
    call transpose_x_to_y(div_ccc_xpencil, accc_ypencil, dm%dccc)
    call transpose_y_to_z(accc_ypencil, div_ccc_zpencil, dm%dccc)

    qzdz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk2
    call compute_qzdz_ccc_zpencil(qzdz_ccc_zpencil, fl%wk3, fl%wk4)

    accc_zpencil (1:ncccz(1),1:ncccz(2),1:ncccz(3))    => fl%wk1
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(div_ccc_zpencil, dm%dccc, dm%rc, 1, IPENCIL(3))
      qyriy_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk3
      call compute_qyriy_ccc_zpencil(qyriy_ccc_zpencil, fl%wk4, fl%wk5)
      mu_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk4
      call compute_mu_ccc_zpencil(mu_ccc_zpencil, fl%wk5)
      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncccz(1); ny = ncccz(2); nz = ncccz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accc_zpencil(i,j,k) = -ONE_THIRD * div_ccc_zpencil(i,j,k) * TWO * mu_ccc_zpencil(i,j,k) &
                              + (qzdz_ccc_zpencil(i,j,k) + qyriy_ccc_zpencil(i,j,k)) * TWO &
                              * (mu_ccc_zpencil(i,j,k) + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncccz(1); ny = ncccz(2); nz = ncccz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accc_zpencil(i,j,k) = (qzdz_ccc_zpencil(i,j,k) - ONE_THIRD * div_ccc_zpencil(i,j,k) &
                                 + qyriy_ccc_zpencil(i,j,k)) * TWO * mu_ccc_zpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    else
      mu_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk4
      call compute_mu_ccc_zpencil(mu_ccc_zpencil, fl%wk5)
      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncccz(1); ny = ncccz(2); nz = ncccz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accc_zpencil(i,j,k) = -ONE_THIRD * div_ccc_zpencil(i,j,k) * TWO * mu_ccc_zpencil(i,j,k) &
                              + qzdz_ccc_zpencil(i,j,k) * TWO * (mu_ccc_zpencil(i,j,k) &
                              + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncccz(1); ny = ncccz(2); nz = ncccz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accc_zpencil(i,j,k) = (qzdz_ccc_zpencil(i,j,k) - ONE_THIRD * div_ccc_zpencil(i,j,k)) &
                              * TWO * mu_ccc_zpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    end if

    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk2
    call Get_z_1der_C2P_3D(accc_zpencil, accp_zpencil, dm, dm%iAccuracy, mbcz_tau3, fbcz_cc4)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accp_zpencil, dm%dccp, dm%rci, 2, IPENCIL(3))

    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk1
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk3
    call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccc)
    call transpose_y_to_x(accp_ypencil, accp_xpencil, dm%dccc)

    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      fl%mz_rhs(i,j,k) = fl%mz_rhs(i,j,k) + accp_xpencil(i,j,k) * fl%rre
    end do; end do; end do
    !$acc end parallel loop

!!   FIXME: check the cause of last digit difference between the CPU and GPU versions
!!   difference from CPU at acpp_ypencil is only limited to the last digit
!!   difference from CPU at accp_ypencil has propagates to the second last digit 
!!   Note: disappear for cylindrical
!!         reappear at mu_ccc_zpencil for cylindrical with thermal

!!------------------------------------------------------------------------------
!!  Z-mom diffusion term 4/4 at (i, j, k')
!!  diff-r-m3 = 1/r * [muiyz * (qzdy + 1/r2 * qydz - qzriy)]^y  (method 1)
!!            = 1/r * [muiyz * (qzdy + 1/r2 * qydz)]^y - 1/r^2 * muiz * qz (method 2)
!!------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      muiyz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk2
      call compute_muiyz_cpp_zpencil(muiyz_cpp_zpencil, fl%wk1, fl%wk3, fl%wkbc1)
      muiyz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk1
      call transpose_z_to_y(muiyz_cpp_zpencil, muiyz_cpp_ypencil, dm%dcpp)
      qzriy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))  => fl%wk2
      call compute_qzriy_cpp_ypencil(qzriy_cpp_ypencil, fl%wk3, fl%wk4)
      qyr2dz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk3
      qyr2dz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk4
      call compute_qyr2dz_cpp_zpencil(qyr2dz_cpp_zpencil, fl%wk3, fl%wk5)
      call transpose_z_to_y(qyr2dz_cpp_zpencil, qyr2dz_cpp_ypencil, dm%dcpp)
      qzdy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))   => fl%wk4
      call compute_qzdy_cpp_ypencil(qzdy_cpp_ypencil, fl%wk5)
      acpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))       => fl%wk1
      if(dm%outlet_sponge_layer(1) > MINP) then
        nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_ypencil(i,j,k) = (qzdy_cpp_ypencil(i,j,k) + qyr2dz_cpp_ypencil(i,j,k) &
                                 - qzriy_cpp_ypencil(i,j,k)) * (muiyz_cpp_ypencil(i,j,k) &
                                 + safe_divide(fl%rre_sponge_c(i),fl%rre))
        end do; end do; end do
        !$acc end parallel loop
      else
        nx = ncppy(1); ny = ncppy(2); nz = ncppy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpp_ypencil(i,j,k) = (qzdy_cpp_ypencil(i,j,k) + qyr2dz_cpp_ypencil(i,j,k) &
                                 - qzriy_cpp_ypencil(i,j,k)) * muiyz_cpp_ypencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if

      fbcy_c4p(1:ncppy(1),1:4,1:ncppy(3)) => fl%wkbc1
      if(is_fbcy_velo_required) then
        call extract_dirichlet_fbcy(fbcy_c4p, acpp_ypencil, dm%dcpp, dm, is_reversed = .true.)
      else
        nx = ncppy(1); nz = ncppy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,4; do i=1,nx
          fbcy_c4p(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if

      accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk2
      call Get_y_midp_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcr_tau3, fbcy_c4p)
      if(dm%icoordinate == ICYLINDRICAL) &
      call multiple_cylindrical_rn(accp_ypencil, dm%dccp, dm%rci, 1, IPENCIL(2))

      accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk2
      call transpose_y_to_x(accp_ypencil, accp_xpencil, dm%dccp)

      nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%mz_rhs(i,j,k) = fl%mz_rhs(i,j,k) + accp_xpencil(i,j,k) * fl%rre
      end do; end do; end do
      !$acc end parallel loop

!!    Note: last digit difference appears at muiyz_cpp_ypencil for cylindrical with thermal

    end if

!!------------------------------------------------------------------------------
!!  pressure gradients and all other body forces
!!  Z-mom pressure gradient in z direction, Z-pencil, d(sigma_1 p)
!!  p-m3 = -1/r * dpdz
!!------------------------------------------------------------------------------
    mz_rhs_pfc_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk1

    pres_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))  => fl%wk2
    pres_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))  => fl%wk3
    call transpose_x_to_y(fl%pres, pres_ypencil, dm%dccc)
    call transpose_y_to_z(pres_ypencil, pres_zpencil, dm%dccc)

    call Get_z_1der_C2P_3D(pres_zpencil, mz_rhs_pfc_zpencil, dm, dm%iAccuracy, dm%ibcz_pr, dm%fbcz_pr)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(mz_rhs_pfc_zpencil, dm%dccp, dm%rci, 1, IPENCIL(3))

    nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      mz_rhs_pfc_zpencil(i,j,k) = -mz_rhs_pfc_zpencil(i,j,k) * dm%sigma1p
    end do; end do; end do
    !$acc end parallel loop

!!------------------------------------------------------------------------------
!!  Z-mom gravity in z direction, Z-pencil
!!------------------------------------------------------------------------------
    if(dm%is_thermo) then
      accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3))   => fl%wk2
      !$acc kernels default(present)
      accp_zpencil(:,:,:) = ZERO
      !$acc end kernels

      if(dm%icoordinate == ICYLINDRICAL) then
        acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
        if(fl%igravity == idir .or. fl%igravity == -idir) then
          call gravity_decomposition_to_rz(fl%dDens, fl%igravity, fl%fgravity(idir), &
                          acpc_ypencil, accp_zpencil, fl%wk4, fl%wk5, fl%wkbc1, dm)
        else if(fl%igravity == idir-1 .or. fl%igravity == -idir+1) then
          call gravity_decomposition_to_rz(fl%dDens, fl%igravity, fl%fgravity(idir-1), &
                          acpc_ypencil, accp_zpencil, fl%wk4, fl%wk5, fl%wkbc1, dm)
        end if
      else
        if(fl%igravity == idir .or. fl%igravity == -idir) then
          fbcz_cc4(1:nccpz(1),1:nccpz(2),1:4)    => fl%wkbc1
          nx = nccpy(1); ny = nccpz(2)
          !$acc parallel loop collapse(3) default(present)
          do k=1,4; do j=1,ny; do i=1,nx
            fbcz_cc4(i,j,k) = dm%fbcz_ftp(i,j,k)%d
          end do; end do; end do
          !$acc end parallel loop
          dDens_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk2
          dDens_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk3
          call transpose_x_to_y(fl%dDens, dDens_ypencil, dm%dccc)
          call transpose_y_to_z(dDens_ypencil, dDens_zpencil, dm%dccc)
          call Get_z_midp_C2P_3D(dDens_zpencil, accp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp, fbcz_cc4)
          nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
          !$acc parallel loop collapse(3) default(present)
          do k=1,nz; do j=1,ny; do i=1,nx
            accp_zpencil(i,j,k) = fl%fgravity(idir) * accp_zpencil(i,j,k) 
          end do; end do; end do
          !$acc end parallel loop
        end if
      end if

      nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mz_rhs_pfc_zpencil(i,j,k) =  mz_rhs_pfc_zpencil(i,j,k) + accp_zpencil(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    mz_rhs_pfc_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk2
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))       => fl%wk3
    call transpose_z_to_y(mz_rhs_pfc_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x(accp_ypencil, mz_rhs_pfc_xpencil, dm%dccp)

!!------------------------------------------------------------------------------
!!  Z-mom Lorentz Force in z direction, x-pencil
!!------------------------------------------------------------------------------
    if(dm%is_mhd) then
      nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mz_rhs_pfc_xpencil(i,j,k) = mz_rhs_pfc_xpencil(i,j,k) + fl%lrfz(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

!!  To build up final y-momentum rhs in fractional steps
    call Calculate_momentum_fractional_step(fl%mz_rhs0, fl%mz_rhs, mz_rhs_pfc_xpencil, dm%dccp, dm, isub)

!!===============================================================================
!!  x-pencil : flow drive terms (source terms) in periodic Streamwise flow
!!===============================================================================
    !$acc data create(rhsx_bulk, rhsz_bulk)
    if (fl%idriven == IDRVF_X_MASSFLUX) then
      call Get_volumetric_average_3d(dm, dm%dpcc, fl%mx_rhs, rhsx_bulk, SPACE_AVERAGE, "mx_rhs")
      !$acc update device(rhsx_bulk)
      !$acc kernels default(present)
      fl%mx_rhs(:,:,:) = fl%mx_rhs(:,:,:) - rhsx_bulk
      !$acc end kernels
    else if (fl%idriven == IDRVF_X_DPDX) then
      rhsx_bulk = - HALF * fl%drvfc * dm%tAlpha(isub) * dm%dt
      !$acc update device(rhsx_bulk)
      !$acc kernels default(present)
      fl%mx_rhs(:,:,:) = fl%mx_rhs(:,:,:) - rhsx_bulk
      !$acc end kernels
    else if (fl%idriven == IDRVF_X_TAUW) then
      ! to add on wall-cells only
      if(nrank == 0) call Print_error_msg('Function to be added soon.')
    else if (fl%idriven == IDRVF_Z_MASSFLUX) then
      call Get_volumetric_average_3d(dm, dm%dccp, fl%mz_rhs, rhsz_bulk, SPACE_AVERAGE, "mz_rhs")
      !$acc update device(rhsz_bulk)
      !$acc kernels default(present)
      fl%mz_rhs(:,:,:) = fl%mz_rhs(:,:,:) - rhsz_bulk
      !$acc end kernels
    else if (fl%idriven == IDRVF_Z_DPDZ) then
      rhsz_bulk = - HALF * fl%drvfc * dm%tAlpha(isub) * dm%dt
      !$acc update device(rhsz_bulk)
      !$acc kernels default(present)
      fl%mz_rhs(:,:,:) = fl%mz_rhs(:,:,:) - rhsz_bulk
      !$acc end kernels
    else if (fl%idriven == IDRVF_Z_TAUW) then
      ! to add on wall-cells only
      if(nrank == 0) call Print_error_msg('Function to be added soon.')
    else
    end if
    !$acc end data

  return

  contains

!!compute qxiy_ppc_ypencil
  subroutine compute_qxiy_ppc_ypencil(qxiy_ppc_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qxiy_ppc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: qx_pcc_ypencil

    qx_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => wk1
    call transpose_x_to_y (fl%qx, qx_pcc_ypencil, dm%dpcc)
    call Get_y_midp_C2P_3D(qx_pcc_ypencil, qxiy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_qx)

  end subroutine

!!compute gxiy_ppc_ypencil
  subroutine compute_gxiy_ppc_ypencil(gxiy_ppc_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: gxiy_ppc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: gx_pcc_ypencil
    
    gx_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => wk1
    call transpose_x_to_y (fl%gx, gx_pcc_ypencil, dm%dpcc)
    call Get_y_midp_C2P_3D(gx_pcc_ypencil, gxiy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_gx)

  end subroutine

!!compute qxdy_ppc_ypencil
  subroutine compute_qxdy_ppc_ypencil(qxdy_ppc_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qxdy_ppc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: qx_pcc_ypencil

    qx_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => wk1
    call transpose_x_to_y(fl%qx, qx_pcc_ypencil, dm%dpcc)
    call Get_y_1der_C2P_3D(qx_pcc_ypencil, qxdy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_qx)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(qxdy_ppc_ypencil, dm%dppc, dm%rp, 1, IPENCIL(2))

  end subroutine

!!compute qxdz_pcp_zpencil
  subroutine compute_qxdz_pcp_zpencil(qxdz_pcp_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qxdz_pcp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qx_pcc_ypencil, qx_pcc_zpencil

    qx_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => wk1
    qx_pcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => wk2
    call transpose_x_to_y(fl%qx, qx_pcc_ypencil, dm%dpcc)
    call transpose_y_to_z(qx_pcc_ypencil, qx_pcc_zpencil, dm%dpcc)
    call Get_z_1der_C2P_3D(qx_pcc_zpencil, qxdz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_qx)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(qxdz_pcp_zpencil, dm%dpcp, dm%rci, 1, IPENCIL(3))

  end subroutine

!!compute qxiz_pcp_zpencil
  subroutine compute_qxiz_pcp_zpencil(qxiz_pcp_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qxiz_pcp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qx_pcc_ypencil, qx_pcc_zpencil

    qx_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => wk1
    qx_pcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => wk2
    call transpose_x_to_y (fl%qx, qx_pcc_ypencil, dm%dpcc)
    call transpose_y_to_z (qx_pcc_ypencil, qx_pcc_zpencil, dm%dpcc)
    call Get_z_midp_C2P_3D(qx_pcc_zpencil, qxiz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_qx)

  end subroutine

!!compute gxiz_pcp_xpencil
  subroutine compute_gxiz_pcp_xpencil(gxiz_pcp_xpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: gxiz_pcp_xpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: apcc_ypencil, gx_pcc_zpencil, gxiz_pcp_zpencil

    gx_pcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3))   => wk1
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3))     => wk2
    call transpose_x_to_y (fl%gx, apcc_ypencil, dm%dpcc)
    call transpose_y_to_z (apcc_ypencil, gx_pcc_zpencil, dm%dpcc)
    gxiz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => wk2
    call Get_z_midp_C2P_3D(gx_pcc_zpencil, gxiz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_gx)
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3))     => wk1
    call transpose_z_to_y (gxiz_pcp_zpencil, apcc_ypencil, dm%dpcp)
    call transpose_y_to_x (apcc_ypencil, gxiz_pcp_xpencil, dm%dpcp)

  end subroutine

!!compute_qyr_ypencil
  subroutine compute_qyr_ypencil(qyr_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qyr_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: qyr_xpencil

    qyr_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => wk1
    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      qyr_xpencil(i,j,k) = fl%qy(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    call multiple_cylindrical_rn(qyr_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1))
    call transpose_x_to_y(qyr_xpencil, qyr_ypencil, dm%dcpc)

    call axis_estimating_radial_xpx(qyr_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(2), is_reversed = .true.)

  end subroutine

!!compute qyiy_ccc_ypencil
  subroutine compute_qyiy_ccc_ypencil(qyiy_ccc_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qyiy_ccc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: qy_cpc_ypencil

    qy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    call transpose_x_to_y(fl%qy, qy_cpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(qy_cpc_ypencil, qyiy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)

  end subroutine

!!compute qyriy_ccc_zpencil
  subroutine compute_qyriy_ccc_zpencil(qyriy_ccc_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qyriy_ccc_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qyr_ypencil, qyriy_ccc_ypencil

    qyr_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))       => wk1
    call compute_qyr_ypencil(qyr_ypencil, wk2)
    qyriy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => wk2
    call Get_y_midp_P2C_3D(qyr_ypencil, qyriy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
    call transpose_y_to_z(qyriy_ccc_ypencil, qyriy_ccc_zpencil, dm%dccc)

  end subroutine

!!compute qyr2dz_cpp_zpencil
  subroutine compute_qyr2dz_cpp_zpencil(qyr2dz_cpp_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qyr2dz_cpp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: acpc_xpencil, qyr2_ypencil

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => wk1
    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      acpc_xpencil(i,j,k) = fl%qy(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 2, IPENCIL(1))
    qyr2_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
    call transpose_x_to_y(acpc_xpencil, qyr2_ypencil, dm%dcpc)
    call axis_estimating_radial_xpx(qyr2_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(2), is_reversed = .true.)
    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
    call transpose_y_to_z(qyr2_ypencil, acpc_zpencil, dm%dcpc)
    call Get_z_1der_C2P_3D(acpc_zpencil, qyr2dz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qyr)

  end subroutine

!!compute qyr2_ypencil
  subroutine compute_qyr2_ypencil(qyr2_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qyr2_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: acpc_xpencil

    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => wk1
    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      acpc_xpencil(i,j,k) = fl%qy(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 2, IPENCIL(1))
    call transpose_x_to_y(acpc_xpencil, qyr2_ypencil, dm%dcpc)
    call axis_estimating_radial_xpx(qyr2_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(2), is_reversed = .true.)

  end subroutine

!!compute qyrdz_cpp_ypencil
  subroutine compute_qyrdz_cpp_ypencil(qyrdz_cpp_ypencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qyrdz_cpp_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qyr_ypencil, acpc_zpencil, acpp_zpencil

    qyr_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))  => wk1
    call compute_qyr_ypencil(qyr_ypencil, wk2)
    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk2
    call transpose_y_to_z(qyr_ypencil, acpc_zpencil, dm%dcpc)
    acpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => wk1
    call Get_z_1der_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qyr)
    call transpose_z_to_y(acpp_zpencil, qyrdz_cpp_ypencil, dm%dcpp)

  end subroutine

!!compute gyiy_ccc_ypencil
  subroutine compute_gyiy_ccc_ypencil(gyiy_ccc_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: gyiy_ccc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: gy_cpc_ypencil

    gy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    call transpose_x_to_y(fl%gy, gy_cpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(gy_cpc_ypencil, gyiy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_gy)

  end subroutine

!!compute qyiz_cpp_zpencil
  subroutine compute_qyiz_cpp_zpencil(qyiz_cpp_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qyiz_cpp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qy_cpc_ypencil, qy_cpc_zpencil

    qy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    qy_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk2
    call transpose_x_to_y(fl%qy, qy_cpc_ypencil, dm%dcpc)
    call transpose_y_to_z(qy_cpc_ypencil, qy_cpc_zpencil, dm%dcpc)
    call Get_z_midp_C2P_3D(qy_cpc_zpencil, qyiz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy)

  end subroutine

!!compute qydy_cpc_ypencil
  subroutine compute_qydy_cpc_ypencil(qydy_cpc_ypencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qydy_cpc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qy_cpc_ypencil, qy_ccc_ypencil

    qy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    call transpose_x_to_y (fl%qy, qy_cpc_ypencil, dm%dcpc)
    qy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => wk2
    call Get_y_midp_P2C_3D(qy_cpc_ypencil, qy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
    call Get_y_1der_C2P_3D(qy_ccc_ypencil, qydy_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)

  end subroutine

!!compute qydz_cpp_zpencil
  subroutine compute_qydz_cpp_zpencil(qydz_cpp_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qydz_cpp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qy_cpc_ypencil, qy_cpc_zpencil

    qy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    qy_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk2
    call transpose_x_to_y(fl%qy, qy_cpc_ypencil, dm%dcpc)
    call transpose_y_to_z(qy_cpc_ypencil, qy_cpc_zpencil, dm%dcpc)
    call Get_z_1der_C2P_3D(qy_cpc_zpencil, qydz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy)

  end subroutine

!!compute gyiz_cpp_ypencil
  subroutine compute_gyiz_cpp_ypencil(gyiz_cpp_ypencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: gyiz_cpp_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: acpc_ypencil, gy_cpc_zpencil, gyiz_cpp_zpencil

    gy_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))   => wk2
    call transpose_x_to_y(fl%gy, acpc_ypencil, dm%dcpc)
    call transpose_y_to_z(acpc_ypencil, gy_cpc_zpencil, dm%dcpc)
    gyiz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => wk2
    call Get_z_midp_C2P_3D(gy_cpc_zpencil, gyiz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_gy)
    call transpose_z_to_y(gyiz_cpp_zpencil, gyiz_cpp_ypencil, dm%dcpp)

  end subroutine

!!compute qziy_cpp_ypencil
  subroutine compute_qziy_cpp_ypencil(qziy_cpp_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qziy_cpp_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: qz_ccp_ypencil

    qz_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk1
    call transpose_x_to_y(fl%qz, qz_ccp_ypencil, dm%dccp)
    call Get_y_midp_C2P_3D(qz_ccp_ypencil, qziy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz)

  end subroutine

!!compute qzriy_cpp_ypencil
  subroutine compute_qzriy_cpp_ypencil(qzriy_cpp_ypencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qzriy_cpp_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: accp_xpencil, accp_ypencil

    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => wk1
    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accp_xpencil(i,j,k) = fl%qz(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rci, 1, IPENCIL(1))
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk2
    call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%dccp)
    call Get_y_midp_C2P_3D(accp_ypencil, qzriy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qzr)

  end subroutine

!!compute qziy_cpp_zpencil
  subroutine compute_qziy_cpp_zpencil(qziy_cpp_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qziy_cpp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qziy_cpp_ypencil

    qziy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => wk1
    call compute_qziy_cpp_ypencil(qziy_cpp_ypencil, wk2)
    call transpose_y_to_z(qziy_cpp_ypencil, qziy_cpp_zpencil, dm%dcpp)

  end subroutine

!!compute gziy_cpp_zpencil
  subroutine compute_gziy_cpp_zpencil(gziy_cpp_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: gziy_cpp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: gz_ccp_ypencil, gziy_cpp_ypencil

    gz_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))   => wk1
    call transpose_x_to_y(fl%gz, gz_ccp_ypencil, dm%dccp)
    gziy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => wk2
    call Get_y_midp_C2P_3D(gz_ccp_ypencil, gziy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_gz)
    call transpose_y_to_z(gziy_cpp_ypencil, gziy_cpp_zpencil, dm%dcpp)

  end subroutine

!!compute qziz_ccc_zpencil
  subroutine compute_qziz_ccc_zpencil(qziz_ccc_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qziz_ccc_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qz_ccp_ypencil, qz_ccp_zpencil

    qz_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk1
    qz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => wk2
    call transpose_x_to_y(fl%qz, qz_ccp_ypencil, dm%dccp)
    call transpose_y_to_z(qz_ccp_ypencil, qz_ccp_zpencil, dm%dccp)
    call Get_z_midp_P2C_3D(qz_ccp_zpencil, qziz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)

  end subroutine

!!compute gziz_ccc_zpencil
  subroutine compute_gziz_ccc_zpencil(gziz_ccc_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: gziz_ccc_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: gz_ccp_ypencil, gz_ccp_zpencil

    gz_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk1
    gz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => wk2
    call transpose_x_to_y(fl%gz, gz_ccp_ypencil, dm%dccp)
    call transpose_y_to_z(gz_ccp_ypencil, gz_ccp_zpencil, dm%dccp)
    call Get_z_midp_P2C_3D(gz_ccp_zpencil, gziz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gz)

  end subroutine

!!compute qzdy_cpp_ypencil
  subroutine compute_qzdy_cpp_ypencil(qzdy_cpp_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qzdy_cpp_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: qz_ccp_ypencil

    qz_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk1
    call transpose_x_to_y(fl%qz, qz_ccp_ypencil, dm%dccp)
    call Get_y_1der_C2P_3D(qz_ccp_ypencil, qzdy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz)

  end subroutine

!!compute qzrdz_cpc_ypencil
  subroutine compute_qzrdz_cpc_ypencil(qzrdz_cpc_ypencil, wk1, wk2, wk3)

    real(WP), pointer, dimension(:,:,:) :: qzrdz_cpc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2, wk3

    real(WP), pointer, dimension(:,:,:) :: qzriy_cpp_ypencil, qzriy_cpp_zpencil, qzrdz_cpc_zpencil

    qzriy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => wk1
    call compute_qzriy_cpp_ypencil(qzriy_cpp_ypencil, wk2, wk3)
    qzriy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => wk2
    call transpose_y_to_z(qzriy_cpp_ypencil, qzriy_cpp_zpencil, dm%dcpp)
    qzrdz_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
    call Get_z_1der_P2C_3D(qzriy_cpp_zpencil, qzrdz_cpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
    call transpose_z_to_y(qzrdz_cpc_zpencil, qzrdz_cpc_ypencil, dm%dcpc)

  end subroutine

!!compute div_ccc_xpencil
  subroutine compute_div_ccc_xpencil(div_ccc_xpencil, wk1, wk2, wk3)

    real(WP), pointer, dimension(:,:,:) :: div_ccc_xpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2, wk3

    real(WP), pointer, dimension(:,:,:) :: qxdx_ccc_xpencil, qydy_ccc_xpencil, &
                                           qydy_ccc_ypencil, qzdz_ccc_xpencil, &
                                           qzdz_ccc_ypencil, qzdz_ccc_zpencil

    qzdz_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => wk1
    qzdz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => wk2
    qzdz_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => wk3
    call compute_qzdz_ccc_zpencil(qzdz_ccc_zpencil, wk2, wk3)
    call transpose_z_to_y (qzdz_ccc_zpencil, qzdz_ccc_ypencil, dm%dccc)
    call transpose_y_to_x (qzdz_ccc_ypencil, qzdz_ccc_xpencil, dm%dccc)

    qydy_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => wk2
    qydy_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => wk3
    call compute_qydy_ccc_ypencil(qydy_ccc_ypencil, wk2)
    call transpose_y_to_x (qydy_ccc_ypencil, qydy_ccc_xpencil, dm%dccc)

    qxdx_ccc_xpencil(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => wk3
    call compute_qxdx_ccc_xpencil(qxdx_ccc_xpencil)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(qxdx_ccc_xpencil, dm%dccc, dm%rc, 1, IPENCIL(1))

    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      div_ccc_xpencil(i,j,k) = qxdx_ccc_xpencil(i,j,k) + qydy_ccc_xpencil(i,j,k) + qzdz_ccc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(div_ccc_xpencil, dm%dccc, dm%rci, 1, IPENCIL(1))

  end subroutine

!!compute qxdx_ccc_xpencil
  subroutine compute_qxdx_ccc_xpencil(qxdx_ccc_xpencil)

    real(WP), pointer, dimension(:,:,:) :: qxdx_ccc_xpencil

    call Get_x_1der_P2C_3D(fl%qx, qxdx_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)

  end subroutine

!!compute qydy_ccc_ypencil
  subroutine compute_qydy_ccc_ypencil(qydy_ccc_ypencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: qydy_ccc_ypencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: qy_cpc_ypencil

    qy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    call transpose_x_to_y (fl%qy, qy_cpc_ypencil, dm%dcpc)
    call Get_y_1der_P2C_3D(qy_cpc_ypencil, qydy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)

  end subroutine

!!compute qzdz_ccc_zpencil
  subroutine compute_qzdz_ccc_zpencil(qzdz_ccc_zpencil, wk1, wk2)

    real(WP), pointer, dimension(:,:,:) :: qzdz_ccc_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2

    real(WP), pointer, dimension(:,:,:) :: qz_ccp_ypencil, qz_ccp_zpencil

    qz_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk1
    qz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => wk2
    call transpose_x_to_y(fl%qz, qz_ccp_ypencil, dm%dccp)
    call transpose_y_to_z(qz_ccp_ypencil, qz_ccp_zpencil, dm%dccp)
    call Get_z_1der_P2C_3D(qz_ccp_zpencil, qzdz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)

  end subroutine

!!compute mu_ccc_xpencil
  subroutine compute_mu_ccc_xpencil(mu_ccc_xpencil)

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_xpencil

    if(.not. dm%is_thermo) then
      nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mu_ccc_xpencil(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mu_ccc_xpencil(i,j,k) = fl%mVisc(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

  end subroutine

!!compute mu_ccc_ypencil
  subroutine compute_mu_ccc_ypencil(mu_ccc_ypencil)

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_ypencil

    if(.not. dm%is_thermo) then
      nx = ncccy(1); ny = ncccy(2); nz = ncccy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mu_ccc_ypencil(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      call transpose_x_to_y(fl%mVisc, mu_ccc_ypencil, dm%dccc)
    end if

  end subroutine

!!compute mu_ccc_ypencil
  subroutine compute_mu_ccc_zpencil(mu_ccc_zpencil, wk1)

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_ypencil

    if(.not. dm%is_thermo) then
      nx = ncccz(1); ny = ncccz(2); nz = ncccz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        mu_ccc_zpencil(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      mu_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => wk1
      call transpose_x_to_y(fl%mVisc, mu_ccc_ypencil, dm%dccc)
      call transpose_y_to_z(mu_ccc_ypencil, mu_ccc_zpencil, dm%dccc)
    end if

  end subroutine

!!compute muixy_ppc_xpencil
  subroutine compute_muixy_ppc_xpencil(muixy_ppc_xpencil, wk1, wk2, wkbc1, wkbc2, wkbc3)

    real(WP), pointer, dimension(:,:,:) :: muixy_ppc_xpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2
    real(WP), pointer, contiguous, dimension(:) :: wkbc1, wkbc2, wkbc3

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_xpencil,   mu_ccc_ypencil,  &
                                           muiy_cpc_xpencil, muiy_cpc_ypencil
    real(WP), pointer, dimension(:,:,:) :: fbcx_4cc, fbcy_c4c, fbcx_4pc


    if(.not. dm%is_thermo) then
      nx = nppcx(1); ny = nppcx(2); nz = nppcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        muixy_ppc_xpencil(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => wkbc1
      ny = ncccx(2); nz = ncccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = dm%fbcx_ftp(i,j,k)%m
      end do; end do; end do
      !$acc end parallel loop

      fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc2
      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = dm%fbcy_ftp(i,j,k)%m
      end do; end do; end do
      !$acc end parallel loop

      mu_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))   => wk1
      muiy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
      call transpose_x_to_y(fl%mVisc, mu_ccc_ypencil, dm%dccc)
      call Get_y_midp_C2P_3D(mu_ccc_ypencil, muiy_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)

      muiy_cpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => wk1
      call transpose_y_to_x(muiy_cpc_ypencil, muiy_cpc_xpencil, dm%dcpc)

      fbcx_4pc(1:4,1:nppcx(2),1:nppcx(3)) => wkbc3
      if(is_fbcx_velo_required) then
        call get_fbcx_ftp_4pc(fbcx_4cc, fbcx_4pc, dm)
      else
        ny = nppcx(2); nz = nppcx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,4
          fbcx_4pc(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if
      call Get_x_midp_C2P_3D(muiy_cpc_xpencil, muixy_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4pc)
    end if

  end subroutine

!!compute muixz_pcp_zpencil
  subroutine compute_muixz_pcp_zpencil(muixz_pcp_zpencil, wk1, wk2, wkbc1)

    real(WP), pointer, dimension(:,:,:) :: muixz_pcp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2
    real(WP), pointer, contiguous, dimension(:) :: wkbc1

    real(WP), pointer, dimension(:,:,:) :: muix_pcc_xpencil, muix_pcc_ypencil, muix_pcc_zpencil
    real(WP), pointer, dimension(:,:,:) :: fbcx_4cc

    if(.not. dm%is_thermo) then
      nx = npcpz(1); ny = npcpz(2); nz = npcpz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        muixz_pcp_zpencil(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => wkbc1
      ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cc(i,j,k) = dm%fbcx_ftp(i,j,k)%m
      end do; end do; end do
      !$acc end parallel loop

      muix_pcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => wk1
      call Get_x_midp_C2P_3D(fl%mVisc, muix_pcc_xpencil, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4cc)
      muix_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => wk2
      call transpose_x_to_y(muix_pcc_xpencil, muix_pcc_ypencil, dm%dpcc)
      muix_pcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => wk1
      call transpose_y_to_z(muix_pcc_ypencil, muix_pcc_zpencil, dm%dpcc)
      call Get_z_midp_C2P_3D(muix_pcc_zpencil, muixz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp)
    end if

  end subroutine

!!compute muiyz_cpp_zpencil
  subroutine compute_muiyz_cpp_zpencil(muiyz_cpp_zpencil, wk1, wk2, wkbc1)

    real(WP), pointer, dimension(:,:,:) :: muiyz_cpp_zpencil
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2
    real(WP), pointer, contiguous, dimension(:) :: wkbc1

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_ypencil, acpc_ypencil, acpc_zpencil
    real(WP), pointer, dimension(:,:,:) :: fbcy_c4c

    if(.not. dm%is_thermo) then
      nx = ncppz(1); ny = ncppz(2); nz = ncppz(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        muiyz_cpp_zpencil(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc1
      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = dm%fbcy_ftp(i,j,k)%m
      end do; end do; end do
      !$acc end parallel loop

      mu_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => wk1
      call transpose_x_to_y(fl%mVisc, mu_ccc_ypencil, dm%dccc)
      acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))   => wk2
      call Get_y_midp_C2P_3D(mu_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)
      acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3))   => wk1
      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)
      call Get_z_midp_C2P_3D(acpc_zpencil, muiyz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp)
    end if

  end subroutine

!!compute fbcy_mu_c4c
  subroutine compute_fbcy_mu_c4c(fbcy_mu_c4c, wk1, wk2, wkbc1)

    real(WP), pointer, dimension(:,:,:) :: fbcy_mu_c4c
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2
    real(WP), pointer, contiguous, dimension(:) :: wkbc1

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_ypencil, muiy_cpc_ypencil
    real(WP), pointer, dimension(:,:,:) :: fbcy_c4c

    if(.not. dm%is_thermo) then
      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_mu_c4c(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc1
      nx = ncpcy(1); nz = ncpcy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4c(i,j,k) = dm%fbcy_ftp(i,j,k)%m
      end do; end do; end do
      !$acc end parallel loop

      mu_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))   => wk1
      call transpose_x_to_y(fl%mVisc, mu_ccc_ypencil, dm%dccc)
      muiy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
      call Get_y_midp_C2P_3D(mu_ccc_ypencil, muiy_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)
      call extract_dirichlet_fbcy(fbcy_mu_c4c, muiy_cpc_ypencil, dm%dcpc, dm, is_reversed = .false.)
    end if

  end subroutine

!!compute fbcy_div_c4c
  subroutine compute_fbcy_div_c4c(fbcy_div_c4c, wk1, wk2, wk3, wkbc1, wkbc2)

    real(WP), pointer, dimension(:,:,:) :: fbcy_div_c4c
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2, wk3
    real(WP), pointer, contiguous, dimension(:) :: wkbc1, wkbc2

    real(WP), pointer, dimension(:,:,:) :: qxiy_ppc_ypencil,  qxiy_ppc_xpencil,  qxdx_cpc_xpencil,  &
                                           qxdx_cpc_ypencil,  qydy_cpc_ypencil,  qzriy_cpp_ypencil, &
                                           qzriy_cpp_zpencil, qzrdz_cpc_zpencil, qzrdz_cpc_ypencil, &
                                           qziy_cpp_ypencil,  qziy_cpp_zpencil,  qzdz_cpc_ypencil,  &
                                           qzdz_cpc_zpencil,  acpc_ypencil
    real(WP), pointer, dimension(:,:,:) :: fbcx_4pc, fbcz_cp4, fbcy_div_c4c1

!!  compute qxdx_cpc_ypencil
    qxiy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => wk1
    call compute_qxiy_ppc_ypencil(qxiy_ppc_ypencil, wk2)
    qxiy_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => wk2
    call transpose_y_to_x(qxiy_ppc_ypencil, qxiy_ppc_xpencil, dm%dppc)
    fbcx_4pc(1:4,1:nppcx(2),1:nppcx(3)) => wkbc1
    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4pc, qxiy_ppc_xpencil, dm%dppc)
    else
      ny = nppcx(2); nz = nppcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4pc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if
    qxdx_cpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => wk1
    call Get_x_1der_P2C_3D(qxiy_ppc_xpencil, qxdx_cpc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, fbcx_4pc)
    qxdx_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
    call transpose_x_to_y(qxdx_cpc_xpencil, qxdx_cpc_ypencil, dm%dcpc)
    fbcy_div_c4c1(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc2
    call extract_dirichlet_fbcy(fbcy_div_c4c1, qxdx_cpc_ypencil, dm%dcpc, dm, is_reversed = .false.)

    nx = ncpcy(1); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,4; do i=1,nx
      fbcy_div_c4c(i,j,k) = fbcy_div_c4c1(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    qydy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    call compute_qydy_cpc_ypencil(qydy_cpc_ypencil, wk2, wk3)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(qydy_cpc_ypencil, dm%dcpc, dm%rpi, 1, IPENCIL(2))
    call axis_estimating_radial_xpx(qydy_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(2), is_reversed = .true.)
    call extract_dirichlet_fbcy(fbcy_div_c4c1, qydy_cpc_ypencil, dm%dcpc, dm, is_reversed = .true.)

    nx = ncpcy(1); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,4; do i=1,nx
      fbcy_div_c4c(i,j,k) = fbcy_div_c4c(i,j,k) + fbcy_div_c4c1(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    if(dm%icoordinate == ICYLINDRICAL) then
      qzriy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => wk1
      call compute_qzriy_cpp_ypencil(qzriy_cpp_ypencil, wk2, wk3)
      qzriy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => wk2
      call transpose_y_to_z(qzriy_cpp_ypencil, qzriy_cpp_zpencil, dm%dcpp)
      qzrdz_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
      call Get_z_1der_P2C_3D(qzriy_cpp_zpencil, qzrdz_cpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
      qzrdz_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
      call transpose_z_to_y(qzrdz_cpc_zpencil, qzrdz_cpc_ypencil, dm%dcpc)
      acpc_ypencil      => qzrdz_cpc_ypencil
    else
      qziy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))  => wk1
      call compute_qziy_cpp_ypencil(qziy_cpp_ypencil, wk2)
      qziy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))  => wk2
      call transpose_y_to_z(qziy_cpp_ypencil, qziy_cpp_zpencil, dm%dcpp)
      fbcz_cp4(1:ncppz(1),1:ncppz(2),1:4) => wkbc2
      if(is_fbcz_velo_required) then
        call extract_dirichlet_fbcz(fbcz_cp4, qziy_cpp_zpencil, dm%dcpp)
      else
        nx = nppcz(1); ny = nppcz(2)
        !$acc parallel loop collapse(3) default(present)
        do k=1,4; do j=1,ny; do i=1,nx
          fbcz_cp4(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if
      qzdz_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
      call Get_z_1der_P2C_3D(qziy_cpp_zpencil, qzdz_cpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, fbcz_cp4)
      qzdz_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
      call transpose_z_to_y(qzdz_cpc_zpencil, qzdz_cpc_ypencil, dm%dcpc)
      acpc_ypencil     => qzdz_cpc_ypencil
    end if

    call extract_dirichlet_fbcy(fbcy_div_c4c1, acpc_ypencil, dm%dcpc, dm, is_reversed = .false.)

    nx = ncpcy(1); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,4; do i=1,nx
      fbcy_div_c4c(i,j,k) = fbcy_div_c4c(i,j,k) + fbcy_div_c4c1(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

  end subroutine

!!compute fbcy_div_c4c (fewer kernels than compute_fbcy_div_c4c)
  subroutine compute_fbcy_div_c4c_1(fbcy_div_c4c, wk1, wk2, wk3, wkbc1, wkbc2, wkbc3)

    real(WP), pointer, dimension(:,:,:) :: fbcy_div_c4c
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2, wk3
    real(WP), pointer, contiguous, dimension(:) :: wkbc1, wkbc2, wkbc3

    real(WP), pointer, dimension(:,:,:) :: qxiy_ppc_ypencil,  qxiy_ppc_xpencil,  qxdx_cpc_xpencil,  &
                                           qxdx_cpc_ypencil,  qydy_cpc_ypencil,  qzriy_cpp_ypencil, &
                                           qzriy_cpp_zpencil, qzrdz_cpc_zpencil, qzrdz_cpc_ypencil, &
                                           qziy_cpp_ypencil,  qziy_cpp_zpencil,  qzdz_cpc_ypencil,  &
                                           qzdz_cpc_zpencil,  acpc_ypencil
    real(WP), pointer, dimension(:,:,:) :: fbcx_4pc, fbcz_cp4, fbcy_div_c4c1, fbcy_div_c4c2, fbcy_div_c4c3

!!  compute qxdx_cpc_ypencil
    qxiy_ppc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => wk1
    call compute_qxiy_ppc_ypencil(qxiy_ppc_ypencil, wk2)
    qxiy_ppc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => wk2
    call transpose_y_to_x(qxiy_ppc_ypencil, qxiy_ppc_xpencil, dm%dppc)
    fbcx_4pc(1:4,1:nppcx(2),1:nppcx(3)) => wkbc1
    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4pc, qxiy_ppc_xpencil, dm%dppc)
    else
      ny = nppcx(2); nz = nppcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4pc(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if
    qxdx_cpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => wk1
    call Get_x_1der_P2C_3D(qxiy_ppc_xpencil, qxdx_cpc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, fbcx_4pc)
    qxdx_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
    call transpose_x_to_y(qxdx_cpc_xpencil, qxdx_cpc_ypencil, dm%dcpc)
    fbcy_div_c4c1(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc2
    call extract_dirichlet_fbcy(fbcy_div_c4c1, qxdx_cpc_ypencil, dm%dcpc, dm, is_reversed = .false.)

    qydy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk1
    call compute_qydy_cpc_ypencil(qydy_cpc_ypencil, wk2, wk3)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(qydy_cpc_ypencil, dm%dcpc, dm%rpi, 1, IPENCIL(2))
    call axis_estimating_radial_xpx(qydy_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(2), is_reversed = .true.)
    fbcy_div_c4c2(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc3
    call extract_dirichlet_fbcy(fbcy_div_c4c2, qydy_cpc_ypencil, dm%dcpc, dm, is_reversed = .true.)

    if(dm%icoordinate == ICYLINDRICAL) then
      qzriy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => wk1
      call compute_qzriy_cpp_ypencil(qzriy_cpp_ypencil, wk2, wk3)
      qzriy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => wk2
      call transpose_y_to_z(qzriy_cpp_ypencil, qzriy_cpp_zpencil, dm%dcpp)
      qzrdz_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
      call Get_z_1der_P2C_3D(qzriy_cpp_zpencil, qzrdz_cpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
      qzrdz_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
      call transpose_z_to_y(qzrdz_cpc_zpencil, qzrdz_cpc_ypencil, dm%dcpc)
      acpc_ypencil      => qzrdz_cpc_ypencil
    else
      qziy_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3))  => wk1
      call compute_qziy_cpp_ypencil(qziy_cpp_ypencil, wk2)
      qziy_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3))  => wk2
      call transpose_y_to_z(qziy_cpp_ypencil, qziy_cpp_zpencil, dm%dcpp)
      fbcz_cp4(1:ncppz(1),1:ncppz(2),1:4) => wkbc1
      if(is_fbcz_velo_required) then
        call extract_dirichlet_fbcz(fbcz_cp4, qziy_cpp_zpencil, dm%dcpp)
      else
        nx = ncppz(1); ny = ncppz(2)
        !$acc parallel loop collapse(3) default(present)
        do k=1,4; do j=1,ny; do i=1,nx
          fbcz_cp4(i,j,k) = MAXP
        end do; end do; end do
        !$acc end parallel loop
      end if
      qzdz_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => wk1
      call Get_z_1der_P2C_3D(qziy_cpp_zpencil, qzdz_cpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, fbcz_cp4)
      qzdz_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => wk2
      call transpose_z_to_y(qzdz_cpc_zpencil, qzdz_cpc_ypencil, dm%dcpc)
      acpc_ypencil     => qzdz_cpc_ypencil
    end if

    fbcy_div_c4c3(1:ncpcy(1),1:4,1:ncpcy(3)) => wkbc1
    call extract_dirichlet_fbcy(fbcy_div_c4c3, acpc_ypencil, dm%dcpc, dm, is_reversed = .false.)

    nx = ncpcy(1); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,4; do i=1,nx
      fbcy_div_c4c(i,j,k) = fbcy_div_c4c1(i,j,k) + fbcy_div_c4c2(i,j,k) + fbcy_div_c4c3(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

  end subroutine

!!compute fbcz_mu_cc4
  subroutine compute_fbcz_mu_cc4(fbcz_mu_cc4, wk1, wk2, wkbc1)

    real(WP), pointer, dimension(:,:,:) :: fbcz_mu_cc4
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2
    real(WP), pointer, contiguous, dimension(:) :: wkbc1

    real(WP), pointer, dimension(:,:,:) :: mu_ccc_ypencil, mu_ccc_zpencil, muiz_ccp_zpencil

    if(.not. dm%is_thermo) then
      nx = nccpz(1); ny = nccpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,4; do j=1,ny; do i=1,nx
        fbcz_mu_cc4(i,j,k) = ONE
      end do; end do; end do
      !$acc end parallel loop
    else
      fbcz_cc4(1:nccpz(1),1:nccpz(2),1:4) => wkbc1
      nx = nccpz(1); ny = nccpz(2)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcz_cc4(i,j,k) = dm%fbcz_ftp(i,j,k)%m
      end do; end do; end do
      !$acc end parallel loop
      mu_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))   => wk1
      mu_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))   => wk2
      call transpose_x_to_y(fl%mVisc, mu_ccc_ypencil, dm%dccc)
      call transpose_x_to_y(mu_ccc_ypencil, mu_ccc_zpencil, dm%dccc)
      muiz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => wk2
      call Get_z_midp_C2P_3D(mu_ccc_zpencil, muiz_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp, fbcz_cc4)
      call extract_dirichlet_fbcz(fbcz_mu_cc4, muiz_ccp_zpencil, dm%dccp)
    end if

  end subroutine

!!compute fbcz_div_cc4
  subroutine compute_fbcz_div_cc4(fbcz_div_cc4, wk1, wk2, wk3, wk4, wk5, wkbc1)

    real(WP), pointer, dimension(:,:,:) :: fbcz_div_cc4
    real(WP), pointer, contiguous, dimension(:) :: wk1, wk2, wk3, wk4, wk5
    real(WP), pointer, contiguous, dimension(:) :: wkbc1

    real(WP), pointer, dimension(:,:,:) :: qxiz_pcp_zpencil, qxiz_pcp_xpencil, qxdx_ccp_xpencil, &
                                           qxdx_ccp_zpencil, qy_zpencil,       qyiz_cpp_zpencil, &
                                           qyiz_cpp_ypencil, qydy_ccp_ypencil, qydy_ccp_zpencil, &
                                           qziz_ccc_zpencil, qzdz_ccp_zpencil, accp_zpencil,     &
                                           apcp_ypencil,     accp_ypencil,     acpc_ypencil

!!  compute qxdx_ccp_zpencil
    qxiz_pcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => wk1
    qxiz_pcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => wk2
    call compute_qxiz_pcp_zpencil(qxiz_pcp_zpencil, wk1, wk3)
    apcp_ypencil(1:npcpz(1),1:npcpz(2),1:npcpz(3))     => wk3
    call transpose_z_to_y(qxiz_pcp_zpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_x(apcp_ypencil, qxiz_pcp_xpencil, dm%dpcp)

    fbcx_4cp(1:4,1:npcpx(2),1:npcpx(3)) => wkbc1
    if(is_fbcx_velo_required) then
      call extract_dirichlet_fbcx(fbcx_4cp, qxiz_pcp_xpencil, dm%dpcp)
    else
      ny = npcpx(2); nz = npcpx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,4
        fbcx_4cp(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    qxdx_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => wk1
    qxdx_ccp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => wk2
    call Get_x_1der_P2C_3D(qxiz_pcp_xpencil, qxdx_ccp_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, fbcx_4cp)
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))     => wk3
    call transpose_x_to_y(qxdx_ccp_xpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_z(accp_ypencil, qxdx_ccp_zpencil, dm%dccp)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(qxdx_ccp_zpencil, dm%dccp, dm%rc, 1, IPENCIL(3))

!!  compute qydy_ccp_zpencil
    qy_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3))       => wk2
    call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
    call transpose_y_to_z(acpc_ypencil, qy_zpencil, dm%dcpc)

    qyiz_cpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => wk3
    call Get_z_midp_C2P_3D(qy_zpencil, qyiz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy)
    qyiz_cpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => wk2
    call transpose_z_to_y(qyiz_cpp_zpencil, qyiz_cpp_ypencil, dm%dcpp)

    fbcy_c4p(1:ncppy(1),1:4,1:ncppy(3)) => wkbc1
    if(is_fbcy_velo_required) then
      call extract_dirichlet_fbcy(fbcy_c4p, qyiz_cpp_ypencil, dm%dcpp, dm, is_reversed = .true.)
    else
      nx = ncppy(1); nz = ncppy(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,4; do i=1,nx
        fbcy_c4p(i,j,k) = MAXP
      end do; end do; end do
      !$acc end parallel loop
    end if

    qydy_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => wk3
    call Get_y_1der_P2C_3D(qyiz_cpp_ypencil, qydy_ccp_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, fbcy_c4p)
    qydy_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => wk2
    call transpose_y_to_z(qydy_ccp_ypencil, qydy_ccp_zpencil, dm%dccp)

    qzdz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => wk3
    qziz_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => wk4
    call compute_qziz_ccc_zpencil(qziz_ccc_zpencil, wk3, wk5)
    call Get_z_1der_C2P_3D(qziz_ccc_zpencil, qzdz_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)

    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3))     => wk1
    nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accp_zpencil(i,j,k) = qxdx_ccp_zpencil(i,j,k) + qydy_ccp_zpencil(i,j,k) + qzdz_ccp_zpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(accp_zpencil, dm%dccp, dm%rci, 1, IPENCIL(3))
    call extract_dirichlet_fbcz(fbcz_div_cc4, accp_zpencil, dm%dccp)

  end subroutine

  end subroutine Compute_momentum_rhs
!==========================================================================================================
!==========================================================================================================

!==========================================================================================================
!==========================================================================================================
  subroutine Correct_massflux(fl, phi_ccc, dm, isub)
    use udf_type_mod
    use input_general_mod
    use operations
    use parameters_constant_mod
    use cylindrical_rn_mod
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in) :: dm
    integer,        intent(in) :: isub
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(in ) :: phi_ccc

    integer, dimension(3) ::  ncccx, ncccy, ncccz, npccx, ncpcx, ncpcy, nccpx, nccpy, nccpz
    integer :: nx, ny, nz, i, j, k
    real(WP), pointer, dimension(:,:,:):: dphidx_pcc, dphidy_cpc, dphidz_ccp
    real(WP), pointer, dimension(:,:,:):: phi_ccc_ypencil, dphidy_cpc_ypencil, dphidz_ccp_ypencil
    real(WP), pointer, dimension(:,:,:):: pphi_ccc_zpencil, dphidz_ccp_zpencil

    ncccx = dm%dccc%xsz
    ncccy = dm%dccc%ysz
    ncccz = dm%dccc%zsz
    npccx = dm%dpcc%xsz
    ncpcx = dm%dcpc%xsz
    ncpcy = dm%dcpc%ysz
    nccpx = dm%dccp%xsz
    nccpy = dm%dccp%ysz
    nccpz = dm%dccp%zsz
!----------------------------------------------------------------------------------------------------------
!   x-pencil, qx = qx - dt * alpha * d(phi_ccc)/dx
!----------------------------------------------------------------------------------------------------------
    dphidx_pcc(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    call Get_x_1der_C2P_3D(phi_ccc, dphidx_pcc, dm, dm%iAccuracy, dm%ibcx_pr, dm%fbcx_pr)

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    if(dm%is_thermo) then
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%gx(i,j,k) = fl%gx(i,j,k) - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidx_pcc(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%qx(i,j,k) = fl%qx(i,j,k) - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidx_pcc(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if
!----------------------------------------------------------------------------------------------------------
!   y-pencil, qy = qy - dt * alpha * r * d(phi_ccc)/dy
!----------------------------------------------------------------------------------------------------------
    phi_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))    => fl%wk1
    call transpose_x_to_y (phi_ccc, phi_ccc_ypencil, dm%dccc)
    dphidy_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
    call Get_y_1der_C2P_3D(phi_ccc_ypencil, dphidy_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_pr, dm%fbcy_pr)
    !call axis_estimating_radial_xpx(dphidy_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(1))
    dphidy_cpc(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk3
    call transpose_y_to_x (dphidy_cpc_ypencil, dphidy_cpc, dm%dcpc)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(dphidy_cpc, dm%dcpc, dm%rp, 1, IPENCIL(1))

    nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
    if(dm%is_thermo) then
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%gy(i,j,k) = fl%gy(i,j,k) - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidy_cpc(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%qy(i,j,k) = fl%qy(i,j,k) - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidy_cpc(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if
!----------------------------------------------------------------------------------------------------------
!   z-pencil, qz = qz - dt * alpha * 1/r * d(phi_ccc)/dz (if qz = uz)
!----------------------------------------------------------------------------------------------------------
    pphi_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))   => fl%wk2
    call transpose_y_to_z (phi_ccc_ypencil, pphi_ccc_zpencil, dm%dccc)
    dphidz_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk1
    call Get_z_1der_C2P_3D(pphi_ccc_zpencil, dphidz_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_pr, dm%fbcz_pr)
    dphidz_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk2
    call transpose_z_to_y (dphidz_ccp_zpencil, dphidz_ccp_ypencil, dm%dccp)
    dphidz_ccp(1:nccpx(1),1:nccpx(2),1:nccpx(3))         => fl%wk1
    call transpose_y_to_x (dphidz_ccp_ypencil, dphidz_ccp,         dm%dccp)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(dphidz_ccp, dm%dccp, dm%rci, 1, IPENCIL(1))

    nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
    if(dm%is_thermo) then
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%gz(i,j,k) = fl%gz(i,j,k) - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidz_ccp(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        fl%qz(i,j,k) = fl%qz(i,j,k) - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidz_ccp(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    return
  end subroutine Correct_massflux

!==========================================================================================================
  subroutine solve_pressure_poisson(fl, dm, isub)
    use udf_type_mod
    use parameters_constant_mod
    use poisson_interface_mod
    use continuity_eq_mod
    use cylindrical_rn_mod
    use find_max_min_ave_mod
    use typeconvert_mod
    use solver_tools_mod
    implicit none
    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl
    integer,        intent(in)    :: isub
    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div
    real(WP) :: coeff, pres_bulk, mass_imbalance(8)

#ifdef DEBUG_STEPS
    if (nrank == 0) then
      call Print_debug_inline_msg("Calculating the RHS of Poisson equation ...")
    end if
#endif
    !------------------------------------------------------------------
    ! 1. Build RHS components of Poisson equation
    !------------------------------------------------------------------
    ! RHS = d(rho)/dt + div(rho u)
    ! For cylindrical coordinate:
    ! Poisson eq is: r^2 * d2/dx2 + r * d(r * d/dy)/dy + d2/dz2  = &
    !                r^2 * div + r^2 * d(rho)/dt
    !------------------------------------------------------------------
    ! time scaling coefficient
    coeff = ONE / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)
    !$acc data create(div)
    call Get_divergence_flow(div, fl, dm)
    ! final RHS
    if (dm%is_thermo) then
      !$acc kernels default(present)
      fl%pcor(:,:,:) = fl%drhodt(:,:,:) + div(:,:,:)
      !$acc end kernels
    else
      !$acc kernels default(present)
      fl%pcor(:,:,:) = div(:,:,:)
      !$acc end kernels
    end if
    !$acc end data

    !write(*,*) 'mass_imbalance before fft1', mass_imbalance
    if(is_global_mass_correction) then
      call check_global_mass_balance(mass_imbalance, fl%drhodt, dm)
      !$acc data copyin(mass_imbalance)
      !$acc kernels default(present)
      fl%pcor(:,:,:) = fl%pcor(:,:,:) - mass_imbalance(8)/dm%vol
      !$acc end kernels
      !$acc end data
      !write(*,*) 'mass_imbalance before fft', mass_imbalance
    end if
    !call check_global_mass_balance(mass_imbalance, fl%drhodt- mass_imbalance(8)/dm%vol, dm)
    !write(*,*) 'mass_imbalance before fftX', mass_imbalance

    if (dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(fl%pcor, dm%dccc, dm%rc, 2, IPENCIL(1))
    end if
    !$acc kernels default(present)
    fl%pcor(:,:,:) = fl%pcor(:,:,:) * coeff
    !$acc end kernels
    !------------------------------------------------------------------
    ! 2. Solve Poisson equation
    !------------------------------------------------------------------
    call solve_fft_poisson(fl%pcor, fl, dm)
    ! remove drift
    call Get_volumetric_average_3d(dm, dm%dccc, fl%pcor, pres_bulk, SPACE_AVERAGE, "phi")
    !if(nrank==0) write(*,*) 'shifted phi:', pres_bulk
    !$acc kernels default(present)
    fl%pcor(:,:,:) = fl%pcor(:,:,:) - pres_bulk
    !$acc end kernels
    !------------------------------------------------------------------
    ! 3. Pressure correction
    !------------------------------------------------------------------
    !$acc kernels default(present)
    fl%pres(:,:,:) = fl%pres(:,:,:) + fl%pcor(:,:,:)
    !$acc end kernels
    ! Remove pressure drift (zero-mean constraint)
    call Get_volumetric_average_3d(dm, dm%dccc, fl%pres, pres_bulk, SPACE_AVERAGE, "pressure")
    !if(nrank==0) write(*,*) 'shifted pres:', pres_bulk
    !$acc kernels default(present)
    fl%pres(:,:,:) = fl%pres(:,:,:) - pres_bulk
    !$acc end kernels

#ifdef DEBUG_STEPS
    ! call wrt_3d_pt_debug(fl%pres, dm%dccc, fl%iteration, isub, 'pr_updated')
    write(*,*) 'solved_phi:', fl%pcor(2, 1:4, 2)
#endif

    return
  end subroutine solve_pressure_poisson
! !==========================================================================================================
!   subroutine solve_poisson_x2z(fl, dm, isub)
!     use udf_type_mod
!     use parameters_constant_mod
!     use decomp_2d_poisson
!     use transpose_extended_mod
!     use continuity_eq_mod

!     implicit none
!     type(t_domain), intent( in    ) :: dm
!     type(t_flow),   intent( inout ) :: fl                  
!     integer,        intent( in    ) :: isub

!     real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: rhs
!     real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: rhs_ypencil
!     real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: rhs_zpencil

  

!     real(WP), dimension( dm%dccc%zst(1) : dm%dccc%zen(1), &
!                          dm%dccc%zst(2) : dm%dccc%zen(2), &
!                          dm%dccc%zst(3) : dm%dccc%zen(3) ) :: rhs_zpencil_ggg
!     !integer :: i, j, k, jj, ii

! !==========================================================================================================
! ! RHS of Poisson Eq.
! !==========================================================================================================
!     fl%pcor_zpencil_ggg = ZERO
! !----------------------------------------------------------------------------------------------------------
! ! $d\rho / dt$ at cell centre
! !----------------------------------------------------------------------------------------------------------
!     if (dm%is_thermo) then
!       rhs = ZERO
!       rhs_ypencil = ZERO
!       rhs_zpencil = ZERO
!       rhs_zpencil_ggg = ZERO
!       call Calculate_drhodt(dm, tm%dDens, tm%dDens0, tm%dDensm2, rhs)
!       call transpose_x_to_y(rhs,         rhs_ypencil)
!       call transpose_y_to_z(rhs_ypencil, rhs_zpencil)
!       call zpencil_index_llg2ggg(rhs_zpencil, rhs_zpencil_ggg, dm%dccc)
!       fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg + rhs_zpencil_ggg
!     end if
! !----------------------------------------------------------------------------------------------------------
! ! $d(\rho u_i)) / dx_i $ at cell centre
! !----------------------------------------------------------------------------------------------------------
!     rhs_zpencil_ggg  = ZERO
!     if (dm%is_thermo) then
!       call Get_divergence_x2z(fl%gx, fl%gy, fl%gz, rhs_zpencil_ggg, dm)
!     else
!       call Get_divergence_x2z(fl%qx, fl%qy, fl%qz, rhs_zpencil_ggg, dm)
!     end if
!     fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg + rhs_zpencil_ggg
!     fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)
! !==========================================================================================================
! !   solve Poisson
! !==========================================================================================================
!     call poisson(fl%pcor_zpencil_ggg)
! !==========================================================================================================
! !   convert back RHS from zpencil ggg to xpencil gll
! !==========================================================================================================
!     call zpencil_index_ggg2llg(fl%pcor_zpencil_ggg, rhs_zpencil, dm%dccc)
!     call transpose_z_to_y (rhs_zpencil, rhs_ypencil, dm%dccc)
!     call transpose_y_to_x (rhs_ypencil, fl%pcor,     dm%dccc)

!     return
!   end subroutine
!==========================================================================================================
!> \brief To update the provisional u or rho u.
!>
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  fl            flow field
!> \param[inout]  dm            domain
!> \param[in]     isub         RK sub-iteration
!==========================================================================================================
  subroutine Solve_momentum_eq(fl, dm, isub)
    use bc_convective_outlet_mod
    use boundary_conditions_mod
    use continuity_eq_mod
    use convert_primary_conservative_mod
    use eq_energy_mod
    use find_max_min_ave_mod
    use io_restart_mod
    use mpi_mod
    use parameters_constant_mod
    use solver_tools_mod
    use typeconvert_mod
    use udf_type_mod
#ifdef DEBUG_STEPS
    use io_tools_mod
    use typeconvert_mod
    use wtformat_mod
#endif
    implicit none
    ! arguments
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(inout) :: dm
    integer,        intent(in)    :: isub
    ! local variables
    logical :: flg_bc_conv
    integer :: i
    real(WP) :: pres_bulk, mass_imbalance(8)

    ! set up thermo info based on different time stepping
    ! drho/dt is updated here as it is used in balancing mass in convective outlet
    if(dm%is_thermo) then
      call Calculate_drhodt(fl, dm, isub)
    end if
    ! to set up convective outlet b.c.
    call update_convective_outlet_flow(fl, dm, isub)
    ! Main Momentum RHS
    ! to calculate the rhs of the momenturn equation in stepping method
    call Compute_momentum_rhs(fl, dm, isub)
    if ( .not. dm%is_thermo) then 
      ! to update intermediate (\hat{q}) or (\hat{g})
      !$acc kernels default(present)
      fl%qx(:, :, :) = fl%qx(:, :, :) + fl%mx_rhs(:, :, :)
      fl%qy(:, :, :) = fl%qy(:, :, :) + fl%my_rhs(:, :, :)
      fl%qz(:, :, :) = fl%qz(:, :, :) + fl%mz_rhs(:, :, :)
      !$acc end kernels
      call enforce_velo_from_fbc(dm, fl%qx, fl%qy, fl%qz, dm%fbcx_qx, dm%fbcy_qy, dm%fbcz_qz)
    else if ( dm%is_thermo) then 
      ! to update intermediate (\hat{q}) or (\hat{g})
      !$acc kernels default(present)
      fl%gx(:, :, :) = fl%gx(:, :, :) + fl%mx_rhs(:, :, :)
      fl%gy(:, :, :) = fl%gy(:, :, :) + fl%my_rhs(:, :, :)
      fl%gz(:, :, :) = fl%gz(:, :, :) + fl%mz_rhs(:, :, :)
      !$acc end kernels
      call enforce_velo_from_fbc(dm, fl%gx, fl%gy, fl%gz, dm%fbcx_gx, dm%fbcy_gy, dm%fbcz_gz)
    else
      call Print_error_msg("Error in velocity updating")
    end if
    if(dm%icase == ICASE_PIPE) call update_fbcy_cc_flow_halo(fl, dm)
    !----------------------------------------------------------------
    ! Poisson eq and velocity correction
    !----------------------------------------------------------------
    if(is_RK_proj(isub)) then
      ! to solve Poisson equation
      !if(nrank == 0) call Print_debug_inline_msg("  Solving Poisson Equation ...")
      !call solve_poisson_x2z(fl, dm, isub)
      call solve_pressure_poisson(fl, dm, isub) ! test show above two methods gave the same results.

      ! to update velocity/massflux correction
      !if(nrank == 0) call Print_debug_inline_msg("  Updating velocity/mass flux ...")
      call Correct_massflux(fl, fl%pcor, dm, isub)
      if ( .not. dm%is_thermo) then
        call enforce_velo_from_fbc(dm, fl%qx, fl%qy, fl%qz, dm%fbcx_qx, dm%fbcy_qy, dm%fbcz_qz)
      else
        call enforce_velo_from_fbc(dm, fl%gx, fl%gy, fl%gz, dm%fbcx_gx, dm%fbcy_gy, dm%fbcz_gz)
      end if
      if(dm%icase == ICASE_PIPE) call update_fbcy_cc_flow_halo(fl, dm)

#ifdef DEBUG_STEPS
      call check_global_mass_balance(mass_imbalance, fl%drhodt, dm)
      write(*,*) 'mass_imbalance after correction', mass_imbalance
#endif
#ifdef DEBUG_STEPS
      if(dm%is_thermo) then
      call wrt_3d_pt_debug(fl%gx, dm%dpcc,   fl%iteration, isub, 'gx_updated') ! debug_ww
      call wrt_3d_pt_debug(fl%gy, dm%dcpc,   fl%iteration, isub, 'gy_updated') ! debug_ww
      call wrt_3d_pt_debug(fl%gz, dm%dccp,   fl%iteration, isub, 'gz_updated') ! debug_ww
      end if
      call Check_element_mass_conservation(fl, dm, opt_isub=isub) 
#endif
    end if
    !----------------------------------------------------------------
    ! to update velocity from gx gy gz 
    !----------------------------------------------------------------
    if(dm%is_thermo) then
      call convert_primary_conservative(fl, dm, fl%dDens, IG2Q, IALL, fl%qx, fl%qy, fl%qz, fl%gx, fl%gy, fl%gz)
    end if

#ifdef DEBUG_STEPS
  call wrt_3d_pt_debug(fl%qx, dm%dpcc,   fl%iteration, isub, 'qx_updated') ! debug_ww
  call wrt_3d_pt_debug(fl%qy, dm%dcpc,   fl%iteration, isub, 'qy_updated') ! debug_ww
  call wrt_3d_pt_debug(fl%qz, dm%dccp,   fl%iteration, isub, 'qz_updated') ! debug_ww
  if(nrank == 0) then
    call Print_debug_inline_msg("Conservative parameters have been updated.")
    ! write(*,*) 'updated qx', fl%qx(1:4, 1, 1)
    ! write(*,*) 'updated qx', fl%qx(1, 1:4, 1)
  end if
#endif

    return
  end subroutine Solve_momentum_eq

end module eq_momentum_mod
