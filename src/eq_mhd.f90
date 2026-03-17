module mhd_mod
! Note: This MHD solver is potential solver only.
!       Assumed: the induced magnetic field is negligible.
  use parameters_constant_mod
  implicit none 

  private :: cross_production_mhd
  public  :: initialise_mhd
  public  :: cleanup_device_mem_mhd
  public  :: compute_Lorentz_force
  public  :: check_current_conservation
contains
!==========================================================================================================
  subroutine initialise_mhd(fl, mh, dm)
    use udf_type_mod
    use math_mod
    use mpi_mod
    use print_msg_mod
    implicit none 
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_mhd),    intent(inout) :: mh

    if(nrank==0) call Print_debug_start_msg('Initialising MHD ...')

    if(mh%is_NStuart) mh%NHartmn = sqrt_wp( ONE/fl%rre * mh%NStuart)
    if(mh%is_NHartmn) mh%NStuart = mh%NHartmn * mh%NHartmn * fl%rre
    mh%iterfrom = fl%iterfrom

    mh%ibcx_bx(:) = dm%ibcx_qx(:)
    mh%ibcx_by(:) = dm%ibcx_qy(:)
    mh%ibcx_bz(:) = dm%ibcx_qz(:)
    mh%ibcy_bx(:) = dm%ibcy_qx(:)
    mh%ibcy_by(:) = dm%ibcy_qy(:)
    mh%ibcy_bz(:) = dm%ibcy_qz(:)
    mh%ibcz_bx(:) = dm%ibcz_qx(:)
    mh%ibcz_by(:) = dm%ibcz_qy(:)
    mh%ibcz_bz(:) = dm%ibcz_qz(:)

    mh%ibcx_jx(:) = dm%ibcx_qx(:)
    mh%ibcx_jy(:) = dm%ibcx_qy(:)
    mh%ibcx_jz(:) = dm%ibcx_qz(:)
    mh%ibcy_jx(:) = dm%ibcy_qx(:)
    mh%ibcy_jy(:) = dm%ibcy_qy(:)
    mh%ibcy_jz(:) = dm%ibcy_qz(:)
    mh%ibcz_jx(:) = dm%ibcz_qx(:)
    mh%ibcz_jy(:) = dm%ibcz_qy(:)
    mh%ibcz_jz(:) = dm%ibcz_qz(:)

    mh%ibcx_ep(:) = dm%ibcx_pr(:)
    mh%ibcy_ep(:) = dm%ibcy_pr(:)
    mh%ibcz_ep(:) = dm%ibcz_pr(:)
    !$acc update device(mh)
!----------------------------------------------------------------------------------------------------------
!   allocate variables
!----------------------------------------------------------------------------------------------------------
    call alloc_x(mh%ep, dm%dccc); mh%ep = ZERO

    call alloc_x(mh%jx, dm%dpcc); mh%jx = ZERO
    call alloc_x(mh%jy, dm%dcpc); mh%jy = ZERO
    call alloc_x(mh%jz, dm%dccp); mh%jz = ZERO

    call alloc_x(mh%bx, dm%dpcc); mh%bx = ZERO
    call alloc_x(mh%by, dm%dcpc); mh%by = ZERO
    call alloc_x(mh%bz, dm%dccp); mh%bz = ZERO

    call alloc_x(fl%lrfx, dm%dpcc); fl%lrfx = ZERO
    call alloc_x(fl%lrfy, dm%dcpc); fl%lrfy = ZERO
    call alloc_x(fl%lrfz, dm%dccp); fl%lrfz = ZERO

    !$acc enter data create(mh%ep)
    !$acc enter data create(mh%jx, mh%jy, mh%jz)
    !$acc enter data create(mh%bx, mh%by, mh%bz)
    !$acc enter data create(fl%lrfx, fl%lrfy, fl%lrfz)

    allocate( mh%fbcx_jx(             4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil
    allocate( mh%fbcy_jx(dm%dpcc%ysz(1),              4, dm%dpcc%ysz(3)) )! default y pencil
    allocate( mh%fbcz_jx(dm%dpcc%zsz(1), dm%dpcc%zsz(2),              4) )! default z pencil

    allocate( mh%fbcx_jy(             4, dm%dcpc%xsz(2), dm%dcpc%xsz(3)) )! default x pencil
    allocate( mh%fbcy_jy(dm%dcpc%ysz(1),              4, dm%dcpc%ysz(3)) )! default y pencil
    allocate( mh%fbcz_jy(dm%dcpc%zsz(1), dm%dcpc%zsz(2),              4) )! default z pencil

    allocate( mh%fbcx_jz(             4, dm%dccp%xsz(2), dm%dccp%xsz(3)) )! default x pencil
    allocate( mh%fbcy_jz(dm%dccp%ysz(1),              4, dm%dccp%ysz(3)) )! default y pencil
    allocate( mh%fbcz_jz(dm%dccp%zsz(1), dm%dccp%zsz(2),              4) )! default z pencil

    allocate( mh%fbcx_bx(             4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil
    allocate( mh%fbcy_bx(dm%dpcc%ysz(1),              4, dm%dpcc%ysz(3)) )! default y pencil
    allocate( mh%fbcz_bx(dm%dpcc%zsz(1), dm%dpcc%zsz(2),              4) )! default z pencil

    allocate( mh%fbcx_by(             4, dm%dcpc%xsz(2), dm%dcpc%xsz(3)) )! default x pencil
    allocate( mh%fbcy_by(dm%dcpc%ysz(1),              4, dm%dcpc%ysz(3)) )! default y pencil
    allocate( mh%fbcz_by(dm%dcpc%zsz(1), dm%dcpc%zsz(2),              4) )! default z pencil

    allocate( mh%fbcx_bz(             4, dm%dccp%xsz(2), dm%dccp%xsz(3)) )! default x pencil
    allocate( mh%fbcy_bz(dm%dccp%ysz(1),              4, dm%dccp%ysz(3)) )! default y pencil
    allocate( mh%fbcz_bz(dm%dccp%zsz(1), dm%dccp%zsz(2),              4) )! default z pencil

    allocate( mh%fbcx_ep(             4, dm%dccc%xsz(2), dm%dccc%xsz(3)) )! default x pencil
    allocate( mh%fbcy_ep(dm%dccc%ysz(1),              4, dm%dccc%ysz(3)) )! default y pencil
    allocate( mh%fbcz_ep(dm%dccc%zsz(1), dm%dccc%zsz(2),              4) )! default z pencil

    !$acc enter data create(mh%fbcx_jx, mh%fbcy_jx, mh%fbcz_jx)
    !$acc enter data create(mh%fbcx_jy, mh%fbcy_jy, mh%fbcz_jy)
    !$acc enter data create(mh%fbcx_jz, mh%fbcy_jz, mh%fbcz_jz)
    !$acc enter data create(mh%fbcx_bx, mh%fbcy_bx, mh%fbcz_bx)
    !$acc enter data create(mh%fbcx_by, mh%fbcy_by, mh%fbcz_by)
    !$acc enter data create(mh%fbcx_bz, mh%fbcy_bz, mh%fbcz_bz)
    !$acc enter data create(mh%fbcx_ep, mh%fbcy_ep, mh%fbcz_ep)
!----------------------------------------------------------------------------------------------------------
!   Since B=(a, b, c) is a uniform field (constant in space), thus all derivatives of B are zero.
!   This means: the the current density introduced by this magnetic field must be zero everywhere, 
!   including at the boundaries.
!----------------------------------------------------------------------------------------------------------
    !$acc kernels default(present)
    mh%bx(:, :, :) = mh%B_static(1)
    mh%by(:, :, :) = mh%B_static(2)
    mh%bz(:, :, :) = mh%B_static(3)

    mh%jx(:, :, :) = ZERO
    mh%jy(:, :, :) = ZERO
    mh%jz(:, :, :) = ZERO

    mh%ep(:, :, :) = ZERO
    !$acc end kernels
!----------------------------------------------------------------------------------------------------------
! Boundary for static magnetic field
!----------------------------------------------------------------------------------------------------------
    !$acc kernels default(present)
    mh%fbcx_bx(:, :, :) = mh%B_static(1)
    mh%fbcy_bx(:, :, :) = mh%B_static(1)
    mh%fbcz_bx(:, :, :) = mh%B_static(1)
    mh%fbcx_by(:, :, :) = mh%B_static(2)
    mh%fbcy_by(:, :, :) = mh%B_static(2)
    mh%fbcz_by(:, :, :) = mh%B_static(2)
    mh%fbcx_bz(:, :, :) = mh%B_static(3)
    mh%fbcy_bz(:, :, :) = mh%B_static(3)
    mh%fbcz_bz(:, :, :) = mh%B_static(3)
    !$acc end kernels
!----------------------------------------------------------------------------------------------------------
! Boundary for current density
!----------------------------------------------------------------------------------------------------------
    !$acc kernels default(present)
    mh%fbcx_jx(:, :, :) = ZERO
    mh%fbcy_jx(:, :, :) = ZERO
    mh%fbcz_jx(:, :, :) = ZERO
    mh%fbcx_jy(:, :, :) = ZERO
    mh%fbcy_jy(:, :, :) = ZERO
    mh%fbcz_jy(:, :, :) = ZERO
    mh%fbcx_jz(:, :, :) = ZERO
    mh%fbcy_jz(:, :, :) = ZERO
    mh%fbcz_jz(:, :, :) = ZERO
    !$acc end kernels
!----------------------------------------------------------------------------------------------------------
! Boundary for electrical potential
!----------------------------------------------------------------------------------------------------------
    !$acc kernels default(present)
    mh%fbcx_ep(:, :, :) = ZERO
    mh%fbcy_ep(:, :, :) = ZERO
    mh%fbcz_ep(:, :, :) = ZERO
    !$acc end kernels

    !call write_visu_mhd(mh, fl, dm, 'initial_mhd')

    if(nrank==0) call Print_debug_end_msg()
    return
  end subroutine
!==========================================================================================================
  subroutine cleanup_device_mem_mhd(fl, mh, dm)
    use udf_type_mod

    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_mhd),    intent(inout) :: mh
    type(t_domain), intent(in)    :: dm

    !$acc exit data delete(mh%ep)
    !$acc exit data delete(mh%jx, mh%jy, mh%jz)
    !$acc exit data delete(mh%bx, mh%by, mh%bz)
    !$acc exit data delete(fl%lrfx, fl%lrfy, fl%lrfz)

    !$acc exit data delete(mh%fbcx_jx, mh%fbcy_jx, mh%fbcz_jx)
    !$acc exit data delete(mh%fbcx_jy, mh%fbcy_jy, mh%fbcz_jy)
    !$acc exit data delete(mh%fbcx_jz, mh%fbcy_jz, mh%fbcz_jz)
    !$acc exit data delete(mh%fbcx_bx, mh%fbcy_bx, mh%fbcz_bx)
    !$acc exit data delete(mh%fbcx_by, mh%fbcy_by, mh%fbcz_by)
    !$acc exit data delete(mh%fbcx_bz, mh%fbcy_bz, mh%fbcz_bz)
    !$acc exit data delete(mh%fbcx_ep, mh%fbcy_ep, mh%fbcz_ep)

    return
  end subroutine
!==========================================================================================================
  subroutine cross_production_mhd(fl, mh, ab_cross_x, ab_cross_y, ab_cross_z, str, dm) ! to add cylindrical
    use udf_type_mod
    use operations
    use print_msg_mod
    use decomp_2d
    implicit none 
    type(t_flow), intent(in) :: fl
    type(t_mhd),  intent(in) :: mh
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(out) :: ab_cross_x
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(out) :: ab_cross_y
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(out) :: ab_cross_z
    character(8), intent(in) :: str

    integer :: n, iacc
    integer :: i, j, k
    integer :: nx, ny, nz
    integer, dimension(3) ::  ncccx, ncccy, ncccz, &
                              npccx, npccy, npccz, &
                              ncpcx, ncpcy, ncpcz, &
                              nccpx, nccpy, nccpz, &
                              nppcx, nppcy, nppcz, &
                              npcpx, npcpy, npcpz, &
                              ncppx, ncppy, ncppz

    real(WP), pointer, dimension(:,:,:) ::  accc_xpencil, &
                                            accc_ypencil, &
                                            accc_zpencil, &
                                            accp_ypencil, &
                                            acpc_ypencil, &
                                            ax_cpc_ypencil, &
                                            bx_cpc_ypencil, &
                                            az_cpc_ypencil, &
                                            bz_cpc_ypencil, &
                                            accp_zpencil, &
                                            ax_ccp_zpencil, &
                                            bx_ccp_zpencil, &
                                            ay_ccp_zpencil, &
                                            by_ccp_zpencil, &
                                            apcc_xpencil, &
                                            ay_pcc_xpencil, &
                                            by_pcc_xpencil, &
                                            az_pcc_xpencil, &
                                            bz_pcc_xpencil, &
                                            accp_xpencil, &
                                            acpc_xpencil, &
                                            apcp_xpencil, &
                                            appc_xpencil, &
                                            appc_ypencil, &
                                            apcc_ypencil, &
                                            apcp_ypencil, &
                                            acpp_ypencil, &
                                            acpc_zpencil, &
                                            acpp_zpencil, &
                                            apcc_zpencil, &
                                            apcp_zpencil

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

    iacc = dm%iAccuracy

    if(trim(str) /= 'ub_cross' .and. trim(str) /= 'jb_cross') then
      call Print_error_msg('The required cross production is not supported.')
    end if
!----------------------------------------------------------------------------------------------------------
! compute ab_cross_x on staggered grid
!----------------------------------------------------------------------------------------------------------
! ay_cpc_xpencil to ay_pcc_xpencil
    appc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk2
    if(trim(str) == 'ub_cross') then
      call Get_x_midp_C2P_3D(fl%qy, appc_xpencil, dm, iacc, dm%ibcx_qy(:), dm%fbcx_qy(:, :, :))
    else
      call Get_x_midp_C2P_3D(mh%jy, appc_xpencil, dm, iacc, mh%ibcx_jy(:), mh%fbcx_jy(:, :, :))
    end if
    appc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk1
    call transpose_x_to_y (appc_xpencil, appc_ypencil, dm%dppc)
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk2
    if(trim(str) == 'ub_cross') then
      call Get_y_midp_P2C_3D(appc_ypencil, apcc_ypencil, dm, iacc, dm%ibcy_qy(:))
    else
      call Get_y_midp_P2C_3D(appc_ypencil, apcc_ypencil, dm, iacc, mh%ibcy_jy(:))
    end if
    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)
    ay_pcc_xpencil => apcc_xpencil

! bz_ccp_xpencil to bz_pcc_xpencil
    apcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk3
    call Get_x_midp_C2P_3D(mh%bz, apcp_xpencil, dm, iacc, mh%ibcx_bz(:), mh%fbcx_bz(:, :, :)) 
    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3)) => fl%wk2
    call transpose_x_to_y (apcp_xpencil, apcp_ypencil, dm%dpcp)
    apcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk3
    call transpose_y_to_z (apcp_ypencil, apcp_zpencil, dm%dpcp)
    apcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => fl%wk2
    call Get_z_midp_P2C_3D(apcp_zpencil, apcc_zpencil, dm, iacc, mh%ibcx_bz(:))
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk3
    call transpose_z_to_y (apcc_zpencil, apcc_ypencil, dm%dpcc)
    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk2
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)
    bz_pcc_xpencil => apcc_xpencil

! az_ccp_xpencil to az_pcc_xpencil
    apcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk4
    if(trim(str) == 'ub_cross') then
      call Get_x_midp_C2P_3D(fl%qz, apcp_xpencil, dm, iacc, dm%ibcx_qz(:), dm%fbcx_qz(:, :, :))
    else
      call Get_x_midp_C2P_3D(mh%jz, apcp_xpencil, dm, iacc, mh%ibcx_jz(:), mh%fbcx_jz(:, :, :))
    end if
    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3)) => fl%wk3
    call transpose_x_to_y (apcp_xpencil, apcp_ypencil, dm%dpcp)
    apcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk4
    call transpose_y_to_z (apcp_ypencil, apcp_zpencil, dm%dpcp)
    apcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => fl%wk3
    if(trim(str) == 'ub_cross') then
      call Get_z_midp_P2C_3D(apcp_zpencil, apcc_zpencil, dm, iacc, dm%ibcz_qz(:))
    else
      call Get_z_midp_P2C_3D(apcp_zpencil, apcc_zpencil, dm, iacc, mh%ibcz_jz(:))
    end if
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk4
    call transpose_z_to_y (apcc_zpencil, apcc_ypencil, dm%dpcc)
    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk3
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)
    az_pcc_xpencil => apcc_xpencil

! by_cpc_xpencil to by_pcc_xpencil
    appc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk5
    call Get_x_midp_C2P_3D(mh%by, appc_xpencil, dm, iacc, mh%ibcx_by(:), mh%fbcx_by(:, :, :))
    appc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk4
    call transpose_x_to_y (appc_xpencil, appc_ypencil, dm%dppc)
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk5
    call Get_y_midp_P2C_3D(appc_ypencil, apcc_ypencil, dm, iacc, mh%ibcy_by(:))
    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk4
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)
    by_pcc_xpencil => apcc_xpencil

    nx = npccx(1); ny = npccx(2); nz = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      ab_cross_x(i,j,k) = ay_pcc_xpencil(i,j,k) * bz_pcc_xpencil(i,j,k) &
                        - az_pcc_xpencil(i,j,k) * by_pcc_xpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

!----------------------------------------------------------------------------------------------------------
! compute ab_cross_y on staggered grid
!----------------------------------------------------------------------------------------------------------
! az_ccp_xpencil to az_cpc_ypencil
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk1
    acpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk2
    if(trim(str) == 'ub_cross') then
      call transpose_x_to_y (fl%qz, accp_ypencil, dm%dccp)
      call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, iacc, dm%ibcy_qz(:), dm%fbcy_qz(:, :, :))
    else
      call transpose_x_to_y (mh%jz, accp_ypencil, dm%dccp)
      call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, iacc, mh%ibcy_jz(:), mh%fbcy_jz(:, :, :))
    end if
    acpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk1
    call transpose_y_to_z (acpp_ypencil, acpp_zpencil, dm%dcpp)
    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => fl%wk2
    if(trim(str) == 'ub_cross') then
      call Get_z_midp_P2C_3D(acpp_zpencil, acpc_zpencil, dm, iacc, dm%ibcz_qz(:))
    else
      call Get_z_midp_P2C_3D(acpp_zpencil, acpc_zpencil, dm, iacc, mh%ibcz_jz(:))
    end if
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk1
    call transpose_z_to_y (acpc_zpencil, acpc_ypencil, dm%dcpc)
    az_cpc_ypencil => acpc_ypencil

! bx_pcc_xpencil to bx_cpc_ypencil
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk2
    call transpose_x_to_y (mh%bx, apcc_ypencil, dm%dpcc)
    appc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk3
    call Get_y_midp_C2P_3D(apcc_ypencil, appc_ypencil, dm, iacc, mh%ibcy_bx(:), mh%fbcy_bx(:, :, :))
    appc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk2
    call transpose_y_to_x (appc_ypencil, appc_xpencil, dm%dppc)
    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk3
    call Get_x_midp_P2C_3D(appc_xpencil, acpc_xpencil, dm, iacc, mh%ibcx_bx(:))
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
    call transpose_x_to_y (acpc_xpencil, acpc_ypencil, dm%dcpc)
    bx_cpc_ypencil => acpc_ypencil

! ax_pcc_xpencil to ax_cpc_ypencil
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk3
    appc_ypencil(1:nppcy(1),1:nppcy(2),1:nppcy(3)) => fl%wk4
    if(trim(str) == 'ub_cross') then
      call transpose_x_to_y (fl%qx, apcc_ypencil, dm%dpcc)
      call Get_y_midp_C2P_3D(apcc_ypencil, appc_ypencil, dm, iacc, dm%ibcy_qx(:), dm%fbcy_qx(:, :, :))
    else
      call transpose_x_to_y (mh%jx, apcc_ypencil, dm%dpcc)
      call Get_y_midp_C2P_3D(apcc_ypencil, appc_ypencil, dm, iacc, mh%ibcy_jx(:), mh%fbcy_jx(:, :, :))
    end if
    appc_xpencil(1:nppcx(1),1:nppcx(2),1:nppcx(3)) => fl%wk3
    call transpose_y_to_x (appc_ypencil, appc_xpencil, dm%dppc)
    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk4
    if(trim(str) == 'ub_cross') then
      call Get_x_midp_P2C_3D(appc_xpencil, acpc_xpencil, dm, iacc, dm%ibcx_qx(:))
    else
      call Get_x_midp_P2C_3D(appc_xpencil, acpc_xpencil, dm, iacc, mh%ibcx_jx(:))
    end if
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
    call transpose_x_to_y (acpc_xpencil, acpc_ypencil, dm%dcpc)
    ax_cpc_ypencil => acpc_ypencil

! bz_ccp_xpencil to bz_cpc_ypencil
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk4
    call transpose_x_to_y (mh%bz, accp_ypencil, dm%dccp)
    acpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk5
    call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, iacc, mh%ibcy_bz(:), mh%fbcy_bz(:, :, :)) 
    acpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk4
    call transpose_y_to_z (acpp_ypencil, acpp_zpencil, dm%dcpp)
    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => fl%wk5
    call Get_z_midp_P2C_3D(acpp_zpencil, acpc_zpencil, dm, iacc, mh%ibcz_bz(:))
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk4
    call transpose_z_to_y (acpc_zpencil, acpc_ypencil, dm%dcpc)
    bz_cpc_ypencil => acpc_ypencil

    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk5
    nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      acpc_ypencil(i,j,k) = az_cpc_ypencil(i,j,k) * bx_cpc_ypencil(i,j,k) &
                          - ax_cpc_ypencil(i,j,k) * bz_cpc_ypencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    call transpose_y_to_x(acpc_ypencil, ab_cross_y,   dm%dcpc)

!----------------------------------------------------------------------------------------------------------
! compute ab_cross_z on staggered grid
!----------------------------------------------------------------------------------------------------------
! ax_pcc_xpencil to ax_ccp_zpencil
    apcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => fl%wk1
    apcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk2
    if(trim(str) == 'ub_cross') then
      call transpose_y_to_z (fl%qx, apcc_zpencil, dm%dpcc)
      call Get_z_midp_C2P_3D(apcc_zpencil, apcp_zpencil, dm, iacc, dm%ibcz_qx(:), dm%fbcz_qx(:, :, :))
    else
      call transpose_y_to_z (mh%jx, apcc_zpencil, dm%dpcc)
      call Get_z_midp_C2P_3D(apcc_zpencil, apcp_zpencil, dm, iacc, mh%ibcz_jx(:), mh%fbcz_jx(:, :, :))
    end if
    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3)) => fl%wk1
    call transpose_z_to_y (apcp_zpencil, apcp_ypencil, dm%dpcp)
    apcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk2
    call transpose_y_to_x (apcp_ypencil, apcp_xpencil, dm%dpcp)
    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk1
    if(trim(str) == 'ub_cross') then
      call Get_x_midp_P2C_3D(apcp_xpencil, accp_xpencil, dm, iacc, dm%ibcx_qx(:))
    else
      call Get_x_midp_P2C_3D(apcp_xpencil, accp_xpencil, dm, iacc, mh%ibcx_jx(:))
    end if
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk2
    call transpose_x_to_y (accp_xpencil, accp_ypencil, dm%dccp)
    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk1
    call transpose_y_to_z (accp_ypencil, accp_zpencil, dm%dccp)
    ax_ccp_zpencil => accp_zpencil

! by_cpc_xpencil to by_ccp_zpencil
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
    call transpose_x_to_y (mh%by, acpc_ypencil, dm%dcpc)
    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => fl%wk2
    call transpose_y_to_z (acpc_ypencil, acpc_zpencil, dm%dcpc) 
    acpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk3
    call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, iacc, mh%ibcz_by(:), mh%fbcz_by(:, :, :))
    acpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk2
    call transpose_z_to_y (acpp_zpencil, acpp_ypencil, dm%dcpp)
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk3
    call Get_y_midp_P2C_3D(acpp_ypencil, accp_ypencil, dm, iacc, mh%ibcy_by(:))
    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk2
    call transpose_y_to_z (accp_ypencil, accp_zpencil, dm%dccp)
    by_ccp_zpencil => accp_zpencil

! ay_cpc_xpencil to ay_ccp_zpencil
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk3
    acpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3)) => fl%wk4
    acpp_zpencil(1:ncppz(1),1:ncppz(2),1:ncppz(3)) => fl%wk5
    if(trim(str) == 'ub_cross') then
      call transpose_x_to_y (fl%qy, acpc_ypencil, dm%dcpc) 
      call transpose_y_to_z (acpc_ypencil, acpc_zpencil, dm%dcpc) 
      call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, iacc, dm%ibcz_qy(:), dm%fbcz_qy(:, :, :))
    else
      call transpose_x_to_y (mh%jy, acpc_ypencil, dm%dcpc) 
      call transpose_y_to_z (acpc_ypencil, acpc_zpencil, dm%dcpc) 
      call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, iacc, mh%ibcz_jy(:), mh%fbcz_jy(:, :, :))
    end if
    acpp_ypencil(1:ncppy(1),1:ncppy(2),1:ncppy(3)) => fl%wk3
    call transpose_z_to_y (acpp_zpencil, acpp_ypencil, dm%dcpp)
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk4
    if(trim(str) == 'ub_cross') then
      call Get_y_midp_P2C_3D(acpp_ypencil, accp_ypencil, dm, iacc, dm%ibcy_qy(:))
    else
      call Get_y_midp_P2C_3D(acpp_ypencil, accp_ypencil, dm, iacc, mh%ibcy_jy(:))
    end if
    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk3
    call transpose_y_to_z (accp_ypencil, accp_zpencil, dm%dccp)
    ay_ccp_zpencil => accp_zpencil

! bx_pcc_xpencil to bx_ccp_zpencil
    apcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk5
    call transpose_x_to_y (mh%bx, apcc_ypencil, dm%dpcc)
    apcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => fl%wk4
    call transpose_y_to_z (apcc_ypencil, apcc_zpencil, dm%dpcc)
    apcp_zpencil(1:npcpz(1),1:npcpz(2),1:npcpz(3)) => fl%wk5
    call Get_z_midp_C2P_3D(apcc_zpencil, apcp_zpencil, dm, iacc, mh%ibcz_bx(:), mh%fbcz_bx(:, :, :))
    apcp_ypencil(1:npcpy(1),1:npcpy(2),1:npcpy(3)) => fl%wk4
    call transpose_z_to_y (apcp_zpencil, apcp_ypencil, dm%dpcp)
    apcp_xpencil(1:npcpx(1),1:npcpx(2),1:npcpx(3)) => fl%wk5
    call transpose_y_to_x (apcp_ypencil, apcp_xpencil, dm%dpcp)
    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk4
    call Get_x_midp_P2C_3D(apcp_xpencil, accp_xpencil, dm, iacc, mh%ibcx_bx(:))
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk5
    call transpose_x_to_y (accp_xpencil, accp_ypencil, dm%dccp)
    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk4
    call transpose_y_to_z (accp_ypencil, accp_zpencil, dm%dccp)
    bx_ccp_zpencil => accp_zpencil

    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk5
    nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz; do j=1,ny; do i=1,nx
      accp_zpencil(i,j,k) = ax_ccp_zpencil(i,j,k) * by_ccp_zpencil(i,j,k) &
                          - ay_ccp_zpencil(i,j,k) * bx_ccp_zpencil(i,j,k)
    end do; end do; end do
    !$acc end parallel loop
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk4
    call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x(accp_ypencil, ab_cross_z,   dm%dccp)

    return
  end subroutine
!==========================================================================================================
  subroutine compute_Lorentz_force(fl, mh, dm)
    use udf_type_mod
    use operations
    use decomp_2d
    use continuity_eq_mod
    use poisson_interface_mod
    use visualisation_field_mod
    implicit none
!----------------------------------------------------------------------------------------------------------
! calculate the Lozrentz-force based on a static magnetic field B, B is time-independent
!----------------------------------------------------------------------------------------------------------
    type(t_flow), intent(inout) :: fl
    type(t_mhd),  intent(inout) :: mh
    type(t_domain), intent(in)  :: dm

    integer :: i, j, k
    integer :: nx1, ny1, nz1
    integer, dimension(3) ::  ncccy, ncccz, npccx, ncpcx, ncpcy, nccpx, nccpy, nccpz

    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)) :: ub_cross_x
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)) :: ub_cross_y
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)) :: ub_cross_z

    real(WP), pointer, dimension(:,:,:) :: apcc_xpencil, &
                                           acpc_xpencil, &
                                           accp_xpencil, &
                                           accc_ypencil, &
                                           acpc_ypencil, &
                                           accp_ypencil, &
                                           accc_zpencil, &
                                           accp_zpencil

!----------------------------------------------------------------------------------------------------------
! calculate vector u cross-product vector b in x-pencil
!----------------------------------------------------------------------------------------------------------
    mh%iteration = fl%iteration

    ncccy = dm%dccc%ysz
    ncccz = dm%dccc%zsz
    npccx = dm%dpcc%xsz
    ncpcx = dm%dcpc%xsz
    ncpcy = dm%dcpc%ysz
    nccpx = dm%dccp%xsz
    nccpy = dm%dccp%ysz
    nccpz = dm%dccp%zsz

    !$acc data create(ub_cross_x, ub_cross_y, ub_cross_z)

    call cross_production_mhd(fl, mh, ub_cross_x, ub_cross_y, ub_cross_z, 'ub_cross', dm)

#ifdef DEBUG_STEPS
    call write_visu_any3darray(ub_cross_x, 'ub_cross_x', 'debug', dm%dpcc, dm, fl%iteration)
    call write_visu_any3darray(ub_cross_y, 'ub_cross_y', 'debug', dm%dcpc, dm, fl%iteration)
    call write_visu_any3darray(ub_cross_z, 'ub_cross_z', 'debug', dm%dccp, dm, fl%iteration)
#endif
!----------------------------------------------------------------------------------------------------------
! calculate div(ub_cross) in x-pencil
!----------------------------------------------------------------------------------------------------------
!!  FIXME: boundary condition issue?
    call Get_divergence_vector(ub_cross_x, ub_cross_y, ub_cross_z, mh%ep, fl, dm)

#ifdef DEBUG_STEPS
    call write_visu_any3darray(mh%ep, 'ep1', 'debug', dm%dccc, dm, fl%iteration)
#endif
!----------------------------------------------------------------------------------------------------------
! solving the Poisson equation for the electric potential
!----------------------------------------------------------------------------------------------------------
    call solve_fft_poisson(mh%ep, fl, dm)   ! TODO: uncomment after testing

#ifdef DEBUG_STEPS
    call write_visu_any3darray(mh%ep, 'ep2', 'debug', dm%dccc, dm, fl%iteration)
#endif
!----------------------------------------------------------------------------------------------------------
! calculate the current density jx, jy, jz (a vector)
!----------------------------------------------------------------------------------------------------------
    apcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk1
    call Get_x_1der_C2P_3D(mh%ep, apcc_xpencil, dm, dm%iAccuracy, mh%ibcx_ep, mh%fbcx_ep)
    nx1 = npccx(1); ny1 = npccx(2); nz1 = npccx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz1; do j=1,ny1; do i=1,nx1
      mh%jx(i,j,k) = - apcc_xpencil(i,j,k) + ub_cross_x(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    accc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk1
    call transpose_x_to_y (mh%ep, accc_ypencil, dm%dccc)
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
    call Get_y_1der_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mh%ibcy_ep, mh%fbcy_ep)
    acpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3)) => fl%wk3
    call transpose_y_to_x (acpc_ypencil, acpc_xpencil, dm%dcpc)
    nx1 = ncpcx(1); ny1 = ncpcx(2); nz1 = ncpcx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz1; do j=1,ny1; do i=1,nx1
      mh%jy(i,j,k) = - acpc_xpencil(i,j,k) + ub_cross_y(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

    accc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk2
    call transpose_y_to_z (accc_ypencil, accc_zpencil, dm%dccc)
    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk1
    call Get_z_1der_C2P_3D(accc_zpencil, accp_zpencil, dm, dm%iAccuracy, mh%ibcz_ep, mh%fbcz_ep)
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk2
    call transpose_z_to_y (accp_zpencil, accp_ypencil, dm%dccp)
    accp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3)) => fl%wk1
    call transpose_y_to_x (accp_ypencil, accp_xpencil, dm%dccp)
    nx1 = nccpx(1); ny1 = nccpx(2); nz1 = nccpx(3)
    !$acc parallel loop collapse(3) default(present)
    do k=1,nz1; do j=1,ny1; do i=1,nx1
      mh%jz(i,j,k) = - accp_xpencil(i,j,k) + ub_cross_z(i,j,k)
    end do; end do; end do
    !$acc end parallel loop

#ifdef DEBUG_STEPS
    call write_visu_any3darray(mh%jx, 'jx', 'debug', dm%dpcc, dm, fl%iteration)
    call write_visu_any3darray(mh%jy, 'jy', 'debug', dm%dcpc, dm, fl%iteration)
    call write_visu_any3darray(mh%jz, 'jz', 'debug', dm%dccp, dm, fl%iteration)
#endif
!----------------------------------------------------------------------------------------------------------
! calculate the Lorentz force lrfx, lrfy, lrfz (a vector)
!----------------------------------------------------------------------------------------------------------
    call cross_production_mhd(fl, mh, fl%lrfx, fl%lrfy, fl%lrfz, 'jb_cross', dm)
!----------------------------------------------------------------------------------------------------------
! un-dimensionlise the value
!----------------------------------------------------------------------------------------------------------
    !$acc kernels default(present)
    fl%lrfx(:,:,:) = fl%lrfx(:,:,:) * mh%Nstuart
    fl%lrfy(:,:,:) = fl%lrfy(:,:,:) * mh%Nstuart
    fl%lrfz(:,:,:) = fl%lrfz(:,:,:) * mh%Nstuart
    !$acc end kernels

#ifdef DEBUG_STEPS
    call write_visu_any3darray(fl%lrfx, 'lrfx', 'debug', dm%dpcc, dm, fl%iteration)
    call write_visu_any3darray(fl%lrfy, 'lrfy', 'debug', dm%dcpc, dm, fl%iteration)
    call write_visu_any3darray(fl%lrfz, 'lrfz', 'debug', dm%dccp, dm, fl%iteration)
#endif

  !$acc end data

  return
  end subroutine

!==========================================================================================================
  subroutine check_current_conservation(mh, fl, dm)
    use find_max_min_ave_mod
    use udf_type_mod
    use decomp_2d
    use wtformat_mod
    use continuity_eq_mod
    implicit none
    type(t_mhd),    intent(in)    :: mh
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in)    :: dm
    real(WP) :: intg_m, intg_fbcx(2), intg_fbcy(2), intg_fbcz(2), crrt_imbalance(8)
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div

    if(.not. dm%is_mhd) return
    !-----------------------------------------------------------------
    ! divergence-free current density flux
    !-----------------------------------------------------------------
    !$acc data create(div)
    call Get_divergence_vector(mh%jx, mh%jy, mh%jz, div, fl, dm)
    call Find_max_min_3d(div, opt_calc='MAXI', opt_name="elementary div(j_vec) =")

    call Get_volumetric_average_3d(dm, dm%dccc, div, intg_m, SPACE_INTEGRAL, 'Ivol')
    !$acc end data
    !-----------------------------------------------------------------
    ! current density flux through b.c. = integral_surface 
    !-----------------------------------------------------------------
    ! x-bc
    intg_fbcx = ZERO
    if(dm%ibcx_qx(1)/=IBC_PERIODIC)then
      call Get_area_average_2d_for_fbcx(dm, dm%dpcc, mh%fbcx_jx, intg_fbcx, SPACE_INTEGRAL, 'fbcx')
    end if
    ! y-bc
    intg_fbcy = ZERO
    if(dm%ibcy_qy(1)/=IBC_PERIODIC)then
      call Get_area_average_2d_for_fbcy(dm, dm%dcpc, mh%fbcy_jy, intg_fbcy, SPACE_INTEGRAL, 'fbcy', is_rf=.true.)
    end if
    ! z-bc
    intg_fbcz = ZERO
    if(dm%ibcz_qz(1)/=IBC_PERIODIC)then
      call Get_area_average_2d_for_fbcz(dm, dm%dccp, mh%fbcz_jz, intg_fbcz, SPACE_INTEGRAL, 'fbcz')
    end if

    ! current change rate
    crrt_imbalance(1:2) = intg_fbcx(1:2)
    crrt_imbalance(3:4) = intg_fbcy(1:2)
    crrt_imbalance(5:6) = intg_fbcz(1:2)
    crrt_imbalance(7)   = intg_m
    crrt_imbalance(8)   = intg_m + &
                          intg_fbcx(1) - intg_fbcx(2) + &
                          intg_fbcy(1) - intg_fbcy(2) + &
                          intg_fbcz(1) - intg_fbcz(2) 

    if (nrank == 0) then
        write (*, wrtfmt1el) 'global electric current imbalance = ', crrt_imbalance(8)
    end if

  end subroutine

end module
