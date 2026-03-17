module continuity_eq_mod
  use operations
  use decomp_2d

  public :: Get_divergence_vector
  public :: Get_divergence_flow
  public :: Get_divergence_vel_x2z
  public :: Check_element_mass_conservation
contains
!==========================================================================================================
!==========================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    div          div(q) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Get_divergence_flow(div, fl, dm)
    use parameters_constant_mod
    use udf_type_mod
    use solver_tools_mod
    use cylindrical_rn_mod
    implicit none

    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent (out) :: div

    integer, dimension(3) ::  ncccx, ncccy, ncccz, npccx, ncpcx, ncpcy, nccpx, nccpy, nccpz
    integer :: nx, ny, nz, i, j, k

    real(WP), pointer, dimension(:,:,:):: qx, qy ,qz
    real(WP), pointer, dimension(:,:,:):: div0x, div0y, div0z, div0_ypencil, div0_zpencil
    real(WP), pointer, dimension(:,:,:):: qy_ypencil, qz_ypencil, qz_zpencil

!   wk1, wk2 and wk3 used as storage, wk4 and wk5 used as workspace buffers

    ncccx = dm%dccc%xsz
    ncccy = dm%dccc%ysz
    ncccz = dm%dccc%zsz
    npccx = dm%dpcc%xsz
    ncpcx = dm%dcpc%xsz
    ncpcy = dm%dcpc%ysz
    nccpx = dm%dccp%xsz
    nccpy = dm%dccp%ysz
    nccpz = dm%dccp%zsz

    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      div(i,j,k) = ZERO
    end do; end do; end do
    !$acc end parallel loop

    div0x(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk1
    div0y(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk2
    div0z(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk3

    qx(1:npccx(1),1:npccx(2),1:npccx(3))    => fl%wk4
    if(dm%is_thermo) then
      nx = npccx(1); ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        qx(i,j,k) = fl%gx(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      nx = npccx(1); ny = npccx(2); nz = npccx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        qx(i,j,k) = fl%qx(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

!!   Note: The following has very similar perfromance, but the above version
!!   is still preferred approach. However, do not use complex structures within
!!   automatic acc kernels/end kernels, as the performance may be very poor.
!    nx = npccx(1); ny = npccx(2); nz = npccx(3)
!    !$acc parallel loop collapse(3) default(present)
!    do k=1,nz; do j=1,ny; do i=1,nx
!      if(dm%is_thermo) then
!        qx(i,j,k) = fl%gx(i,j,k)
!      else
!        qx(i,j,k) = fl%qx(i,j,k)
!      end if
!    end do; end do; end do
!    !$acc end parallel loop

    call Get_x_1der_P2C_3D(qx, div0x, dm, dm%iAccuracy, dm%ibcx_qx(:), dm%fbcx_qx)

    qy(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk4
    if(dm%is_thermo) then
      nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        qy(i,j,k) = fl%gy(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        qy(i,j,k) = fl%qy(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    qy_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))   => fl%wk5
    call transpose_x_to_y(qy, qy_ypencil, dm%dcpc)
    div0_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk4
    call Get_y_1der_P2C_3D(qy_ypencil, div0_ypencil, dm, dm%iAccuracy, dm%ibcy_qy(:), dm%fbcy_qy)
    call transpose_y_to_x(div0_ypencil, div0y, dm%dccc)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(div0y, dm%dccc, dm%rci, 1, IPENCIL(1))

    qz(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk4
    if(dm%is_thermo) then
      nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        qz(i,j,k) = fl%gz(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    else
      nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
      !$acc parallel loop collapse(3) default(present)
      do k=1,nz; do j=1,ny; do i=1,nx
        qz(i,j,k) = fl%qz(i,j,k)
      end do; end do; end do
      !$acc end parallel loop
    end if

    qz_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))   => fl%wk5
    call transpose_x_to_y(qz, qz_ypencil, dm%dccp)
    qz_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3))   => fl%wk4
    call transpose_y_to_z(qz_ypencil, qz_zpencil, dm%dccp)
    div0_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk5
    call Get_z_1der_P2C_3D(qz_zpencil, div0_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)
    div0_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk4
    call transpose_z_to_y(div0_zpencil, div0_ypencil, dm%dccc)
    call transpose_y_to_x(div0_ypencil, div0z,         dm%dccc)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(div0z, dm%dccc, dm%rci, 1, IPENCIL(1))

    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
        div(i,j,k) = div(i,j,k) + div0x(i,j,k) + div0y(i,j,k) + div0z(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
!!  Note: In practice, it is not necessary to initialise div to zero. However,
!!  if the above line is changed to div(i,j,k) = div0x(i,j,k) + div0y(i,j,k) + div0z(i,j,k),
!!  and remove the initialisation, this would cause a calculation difference
!!  at the last digit.

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Get_divergence_vector(ux, uy, uz, div, fl, dm)
    use parameters_constant_mod
    use udf_type_mod
    use cylindrical_rn_mod
    implicit none

    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (in ) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (in ) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (in ) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent (out) :: div

    integer, dimension(3) ::  ncccx, ncccy, ncccz, npccx, ncpcx, ncpcy, nccpx, nccpy, nccpz
    integer :: nx, ny, nz, i, j, k

    real(WP), pointer, dimension(:,:,:):: div0x, div0y, div0z, div0_ypencil, div0_zpencil
    real(WP), pointer, dimension(:,:,:):: uy_ypencil, uz_ypencil, uz_zpencil

!   wk1, wk2 and wk3 used as storage, wk4 and wk5 used as workspace buffers

    ncccx = dm%dccc%xsz
    ncccy = dm%dccc%ysz
    ncccz = dm%dccc%zsz
    npccx = dm%dpcc%xsz
    ncpcx = dm%dcpc%xsz
    ncpcy = dm%dcpc%ysz
    nccpx = dm%dccp%xsz
    nccpy = dm%dccp%ysz
    nccpz = dm%dccp%zsz

    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      div(i,j,k) = ZERO
    end do; end do; end do
    !$acc end parallel loop

    div0x(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk1
    div0y(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk2
    div0z(1:ncccx(1),1:ncccx(2),1:ncccx(3)) => fl%wk3
!----------------------------------------------------------------------------------------------------------
! operation in x pencil, du/dx
!----------------------------------------------------------------------------------------------------------
    call Get_x_1der_P2C_3D(ux, div0x, dm, dm%iAccuracy, dm%ibcx_qx(:), dm%fbcx_qx)
    !write(*,*) 'div, x', div0(1, 1, 1), div0(2, 2, 2), div0(8, 8, 8)!, div0(16, 8, 8), div0(32, 8, 8)
!----------------------------------------------------------------------------------------------------------
! operation in y pencil, dqy/dy * (1/r)
!----------------------------------------------------------------------------------------------------------
    uy_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))   => fl%wk5
    call transpose_x_to_y(uy, uy_ypencil, dm%dcpc)
    div0_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk4
    call Get_y_1der_P2C_3D(uy_ypencil, div0_ypencil, dm, dm%iAccuracy, dm%ibcy_qy(:), dm%fbcy_qy)
    call transpose_y_to_x(div0_ypencil, div0y, dm%dccc)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(div0y, dm%dccc, dm%rci, 1, IPENCIL(1))
    !write(*,*) 'div, y', div0(1, 1, 1), div0(2, 2, 2), div0(8, 8, 8)!, div0(16, 8, 8), div0(32, 8, 8)
!----------------------------------------------------------------------------------------------------------
! operation in z pencil, dw/dz * (1/r)
!----------------------------------------------------------------------------------------------------------
    uz_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))   => fl%wk5
    call transpose_x_to_y(uz,         uz_ypencil, dm%dccp)
    uz_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3))   => fl%wk4
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
    div0_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk5
    call Get_z_1der_P2C_3D(uz_zpencil, div0_zpencil, dm, dm%iAccuracy, dm%ibcz_qz(:), dm%fbcz_qz)
    div0_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk4
    call transpose_z_to_y(div0_zpencil, div0_ypencil, dm%dccc)
    call transpose_y_to_x(div0_ypencil, div0z, dm%dccc)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(div0z, dm%dccc, dm%rci, 1, IPENCIL(1))

    nx = ncccx(1); ny = ncccx(2); nz = ncccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
        div(i,j,k) = div(i,j,k) + div0x(i,j,k) + div0y(i,j,k) + div0z(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    !write(*,*) 'div, z', div0(1, 1, 1), div0(2, 2, 2), div0(8, 8, 8)!, div0(16, 8, 8), div0(32, 8, 8)
    !write(*,*) 'divall', div0(1, 1, 1), div(8, 8, 8)

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Get_divergence_vel_x2z(ux, uy, uz, div_zpencil_ggg, dm)
    use parameters_constant_mod
    use udf_type_mod
    use transpose_extended_mod
    use cylindrical_rn_mod
    use poisson_interface_mod
    use cylindrical_rn_mod
    use decomp_extended_mod
    implicit none

    type(t_domain), intent (in) :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (in) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (in) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (in) :: uz
    real(WP), dimension(dm%dccc%zst(1) : dm%dccc%zen(1), &
                        dm%dccc%zst(2) : dm%dccc%zen(2), &
                        dm%dccc%zst(3) : dm%dccc%zen(3)), intent (out) :: div_zpencil_ggg

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div0
    real(WP), dimension(dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3)) :: div0_ypencil
    real(WP), dimension(dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3)) :: div0_zpencil
    real(WP), dimension(dm%dccc%yst(1) : dm%dccc%yen(1), &
                        dm%dccc%yst(2) : dm%dccc%yen(2), &
                        dm%dccc%ysz(3))                  :: div0_ypencil_ggl
    real(WP), dimension(dm%dccc%yst(1) : dm%dccc%yen(1), &
                        dm%dccc%yst(2) : dm%dccc%yen(2), &
                        dm%dccc%ysz(3))                  :: div_ypencil_ggl
    real(WP), dimension(dm%dccc%zst(1) : dm%dccc%zen(1), &
                        dm%dccc%zst(2) : dm%dccc%zen(2), &
                        dm%dccc%zst(3) : dm%dccc%zen(3)) :: div0_zpencil_ggg

    real(WP), dimension(dm%dcpc%ysz(1),                  dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: uy_ypencil
    !real(WP), dimension(dm%dccp%yst(1) : dm%dccp%yen(1), dm%dccp%ysz(2), dm%dccp%ysz(3)) :: uz_ypencil_ggl

    real(WP), dimension(dm%dccp%ysz(1),                  dm%dccp%ysz(2),                  dm%dccp%ysz(3)) :: uz_ypencil
    real(WP), dimension(dm%dccp%zsz(1),                  dm%dccp%zsz(2),                  dm%dccp%zsz(3)) :: uz_zpencil

!----------------------------------------------------------------------------------------------------------
! operation in x pencil, du/dx
!----------------------------------------------------------------------------------------------------------
    div0 = ZERO
    div0_ypencil_ggl = ZERO
    div_ypencil_ggl = ZERO
    call Get_x_1der_P2C_3D(ux, div0, dm, dm%iAccuracy, dm%ibcx_qx(:), dm%fbcx_qx)
    call transpose_x_to_y(div0, div0_ypencil_ggl, dm%dccc)
    div_ypencil_ggl = div0_ypencil_ggl
!----------------------------------------------------------------------------------------------------------
! operation in y pencil, dv/dy * (1/r)
!----------------------------------------------------------------------------------------------------------
    uy_ypencil = ZERO
    div0_ypencil = ZERO
    div0_ypencil_ggl = ZERO
    call transpose_x_to_y(uy, uy_ypencil, dm%dcpc)
    call Get_y_1der_P2C_3D(uy_ypencil, div0_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(div0_ypencil, dm%dccc, dm%rci, 1, IPENCIL(2))
    call ypencil_index_lgl2ggl(div0_ypencil, div0_ypencil_ggl, dm%dccc)
    div_ypencil_ggl = div_ypencil_ggl + div0_ypencil_ggl
    call transpose_y_to_z(div_ypencil_ggl, div_zpencil_ggg, dm%dccc)
!----------------------------------------------------------------------------------------------------------
! operation in z pencil, dw/dz * (1/r)
!----------------------------------------------------------------------------------------------------------
    uz_ypencil = ZERO
    uz_zpencil = ZERO
    div0_zpencil = ZERO
    div0_zpencil_ggg = ZERO
    call transpose_x_to_y(uz,         uz_ypencil, dm%dccp)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
    call Get_z_1der_P2C_3D(uz_zpencil, div0_zpencil, dm, dm%iAccuracy, dm%ibcz_qz(:), dm%fbcz_qz)
    if(dm%icoordinate == ICYLINDRICAL) &
    call multiple_cylindrical_rn(div0_zpencil, dm%dccc, dm%rci, 1, IPENCIL(3))
    call zpencil_index_llg2ggg(div0_zpencil, div0_zpencil_ggg, dm%dccc)
    div_zpencil_ggg = div_zpencil_ggg + div0_zpencil_ggg

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Check_element_mass_conservation(fl, dm, opt_isub, opt_str)
    use precision_mod
    use udf_type_mod
    use input_general_mod    
    use parameters_constant_mod
    use math_mod                
    use mpi_mod
    use solver_tools_mod
    use wtformat_mod
    use find_max_min_ave_mod
    use typeconvert_mod
    implicit none

    type(t_domain), intent( in) :: dm
    type(t_flow),   intent( inout) :: fl
    integer, intent(in), optional :: opt_isub
    character(*), intent(in), optional :: opt_str

    character(32) :: str
    integer :: n, nlayer, isub
    real(WP) :: mm(2), mm0
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div
    !----------------------------------------------------------------
    ! safe-proof
    !----------------------------------------------------------------
    if(present(opt_str)) then
      str = trim(opt_str)//'_iter_'//int2str(fl%iteration)
    else
      str = 'at iter = '//int2str(fl%iteration)
    end if

    if(present(opt_isub)) then
      isub = opt_isub
    else
      isub = 0
    end if
    !----------------------------------------------------------------
    ! Calculate mass conservation residual
    !----------------------------------------------------------------
    !$acc data create(div)
    ! $d\rho / dt$ at cell centre
    ! $d(\rho u_i)) / dx_i $ at cell centre
    call Get_divergence_flow(div, fl, dm)

    if (dm%is_thermo) then
      !$acc kernels default(present)
      div(:,:,:) = div(:,:,:) + fl%drhodt(:,:,:)
      !$acc end kernels
    end if

#ifdef DEBUG_STEPS
    if(MOD(fl%iteration, dm%visu_nfre) == 0) &
    call write_visu_any3darray(div, 'divU', 'debug'//trim(str), dm%dccc, dm, fl%iteration)
#endif
    !----------------------------------------------------------------
    ! Find Max. mass conservation residual
    !----------------------------------------------------------------
    n = dm%dccc%xsz(1)
    mm0 = fl%mcon(1)
    fl%mcon = ZERO

    if(dm%is_periodic(1)) then
      nlayer = 0
    else
      nlayer = 4
      call Find_max_min_3d(div(1:nlayer, :, :), opt_calc='MAXI', &
            opt_work=mm, opt_name="Mass Consv. (inlet  4) =")
      fl%mcon(2) = mm(2)
      call Find_max_min_3d(div(n-nlayer+1:n, :, :), opt_calc='MAXI', &
            opt_work=mm, opt_name="Mass Consv. (outlet 4) =")
      fl%mcon(3) = mm(2)
    end if
    call Find_max_min_3d(div(nlayer+1:n-nlayer, :, :), opt_calc='MAXI', &
        opt_work=mm, opt_name="Mass Consv. (bulk    ) =")
    fl%mcon(1) = mm(2)
    fl%mcon(4) = safe_divide(fl%mcon(1)-mm0, dabs(mm0))
    if(nrank==0) write(*, '(A,1F9.2,A)') ' ', fl%mcon(4)*100.0_WP, '%'
    ! terminate code once too large
    if(nrank == 0) then
      if(fl%mcon(1) > 1.0_WP .and. fl%iteration > 10000 ) &
      call Print_error_msg("Mass conservation is not strictly satisfied at the machine precision level.")
    end if
    if (nrank == 0) then
      write (*, wrtfmt1el) 'global mass flux imbalance =', fl%tt_mass_change
    end if
    !$acc end data

    !----------------------------------------------------------------
    ! turn on numerical tricks based on mass conservation residual
    !----------------------------------------------------------------
    if(fl%iteration >= 1) then
      if(dm%is_conv_outlet(1) .or. dm%is_conv_outlet(3)) then
        if(.not. is_damping_drhodt) then
          if(fl%mcon(1) > 1.0e-1_WP) then
            is_damping_drhodt = .true.
            if(nrank==0) call Print_warning_msg('drho/dt damping function is on.')
          end if
        end if
      else
        if(.not. is_global_mass_correction) then
          if(fl%tt_mass_change > 1.0e-7_WP) then
            is_global_mass_correction = .true. ! scaled convective b.c. has already met this.
            if(nrank==0) call Print_warning_msg('is_global_mass_correction is True for RHS of Pression Poisson Eq.')
          end if
        end if
      end if
    end if
    if(fl%mcon(4) > ONE) then 
      if(nrank==0) call Print_warning_msg('Mass conservation residual increased by 100%!')
    end if
#ifdef DEBUG_STEPS
    if(nrank == 0) then
      write (*, *) "  Check Mass Conservation:", fl%iteration, isub, fl%mcon(1:4) 
    end if
#endif
    return
  end subroutine Check_element_mass_conservation

end module continuity_eq_mod
