module convert_primary_conservative_mod
  public :: convert_primary_conservative
 contains 

  subroutine convert_primary_conservative(fl, dm, dens, itag, iloc, qx, qy, qz, gx, gy, gz)
    use udf_type_mod
    use operations
    use decomp_2d
    use parameters_constant_mod
    use cylindrical_rn_mod
    implicit none
    type(t_domain), intent(inout)   :: dm
    type(t_flow  ), intent(inout) :: fl
    integer, intent(in) :: itag
    integer, intent(in) :: iloc
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(in) :: dens
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(inout), optional ::  qx, gx
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(inout), optional ::  qy, gy
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(inout), optional ::  qz, gz

    integer  :: i, j, k
    integer  :: nx, ny, nz
    integer, dimension(3) ::  ncccx, ncccy, ncccz, &
                              npccx, npccy, npccz, &
                              ncpcx, ncpcy, ncpcz, &
                              nccpx, nccpy, nccpz

    real(WP), pointer, dimension(:,:,:) :: d_ccc_ypencil, &
                                           d_ccc_zpencil, &
                                           d_pcc_xpencil, &
                                           d_pcc_ypencil, &
                                           d_pcc_zpencil, &
                                           d_cpc_xpencil, &
                                           d_cpc_zpencil, &
                                           d_cpc_ypencil, &
                                           d_ccp_xpencil, &
                                           d_ccp_ypencil, &
                                           d_ccp_zpencil, &

                                           acpc_ypencil,  &
                                           accp_ypencil,  &
                                           accp_zpencil,  &

                                           fbcx_4cc,      &
                                           fbcy_c4c,      &
                                           fbcz_cc4

    if(.not. dm%is_thermo) return

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

    if(iloc == IBLK .or. iloc == IALL) then
      if(.not. present(qx)) call Print_error_msg(' Lack of primary variables')
      if(.not. present(qy)) call Print_error_msg(' Lack of primary variables')
      if(.not. present(qz)) call Print_error_msg(' Lack of primary variables')
      if(.not. present(gx)) call Print_error_msg(' Lack of conservative variables')
      if(.not. present(gy)) call Print_error_msg(' Lack of conservative variables')
      if(.not. present(gz)) call Print_error_msg(' Lack of conservative variables')
    end if
! wk1, wk2 used as workspace buffers, wk3, wk4 and wk5 used as storage
!----------------------------------------------------------------------------------------------------------
! x-pencil : u1 -> g1 = u1_pcc * d_pcc
!----------------------------------------------------------------------------------------------------------
    fbcx_4cc(1:4,1:npccx(2),1:npccx(3)) => fl%wkbc1
    !$acc kernels default(present)
    fbcx_4cc(:, :, :) = dm%fbcx_ftp(:, :, :)%d
    !$acc end kernels
    ! save d_pcc_xpencil to wk3 for future use
    d_pcc_xpencil(1:npccx(1),1:npccx(2),1:npccx(3)) => fl%wk3
    call Get_x_midp_C2P_3D (dens, d_pcc_xpencil, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4cc)
    if(iloc == IBLK .or. iloc == IALL) then
      if(itag == IQ2G) then
        nx = npccx(1); ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          gx(i,j,k) = qx(i,j,k) * d_pcc_xpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
      if(itag == IG2Q) then
        nx = npccx(1); ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          qx(i,j,k) = gx(i,j,k) / d_pcc_xpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
      end if
    end if
!----------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2 = u2_cpc * d_cpc
!----------------------------------------------------------------------------------------------------------
    d_ccc_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3))  => fl%wk1
    call transpose_x_to_y(dens, d_ccc_ypencil, dm%dccc)
    fbcy_c4c(1:ncpcy(1),1:4,1:ncpcy(3)) => fl%wkbc1
    !$acc kernels default(present)
    fbcy_c4c(:, :, :) = dm%fbcy_ftp(:, :, :)%d
    !$acc end kernels
    ! save d_cpc_ypencil to wk4 for future use
    d_cpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3))  => fl%wk4
    call Get_y_midp_C2P_3D (d_ccc_ypencil, d_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)
    call axis_estimating_radial_xpx(d_cpc_ypencil, dm%dcpc, IPENCIL(2), dm, IDIM(1))
    acpc_ypencil(1:ncpcy(1),1:ncpcy(2),1:ncpcy(3)) => fl%wk2
    if(iloc == IBLK .or. iloc == IALL) then
      if(itag == IQ2G) then
        call transpose_x_to_y(qy, acpc_ypencil, dm%dcpc)
        nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpc_ypencil(i,j,k) = acpc_ypencil(i,j,k) * d_cpc_ypencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
        call transpose_y_to_x(acpc_ypencil, gy, dm%dcpc)
      else if(itag == IG2Q) then
        call transpose_x_to_y(gy, acpc_ypencil, dm%dcpc)
        nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          acpc_ypencil(i,j,k) = acpc_ypencil(i,j,k) / d_cpc_ypencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
        call transpose_y_to_x(acpc_ypencil, qy, dm%dcpc)
      else
      end if
    end if
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3 = u3_ccp * d_ccp
!----------------------------------------------------------------------------------------------------------
    d_ccc_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3))  => fl%wk2
    call transpose_y_to_z(d_ccc_ypencil, d_ccc_zpencil, dm%dccc)
    fbcz_cc4(1:nccpz(1),1:nccpz(2),1:4) => fl%wkbc1
    !$acc kernels default(present)
    fbcz_cc4(:, :, :) = dm%fbcz_ftp(:, :, :)%d
    !$acc end kernels
    ! save d_ccp_zpencil to wk5 for future use
    d_ccp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3))  => fl%wk5
    call Get_z_midp_C2P_3D (d_ccc_zpencil, d_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp, fbcz_cc4)
    accp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3)) => fl%wk1
    accp_zpencil(1:nccpz(1),1:nccpz(2),1:nccpz(3)) => fl%wk2
    if(iloc == IBLK .or. iloc == IALL) then
      if(itag == IQ2G) then
        call transpose_x_to_y(qz, accp_ypencil, dm%dccp)
        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
        nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accp_zpencil(i,j,k) = accp_zpencil(i,j,k) * d_ccp_zpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
        call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccp)
        call transpose_y_to_x(accp_ypencil, gz, dm%dccp)
      else if(itag == IG2Q) then
        call transpose_x_to_y(gz, accp_ypencil, dm%dccp)
        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
        nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,nx
          accp_zpencil(i,j,k) = accp_zpencil(i,j,k) / d_ccp_zpencil(i,j,k)
        end do; end do; end do
        !$acc end parallel loop
        call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccp)
        call transpose_y_to_x(accp_ypencil, qz, dm%dccp)
      else
      end if
    end if

!----------------------------------------------------------------------------------------------------------
! BC: - x pencil
!----------------------------------------------------------------------------------------------------------
    if(iloc == IBND .or. iloc == IALL) then
    if(dm%ibcx_qx(1) == IBC_DIRICHLET .or. dm%ibcx_qx(2) == IBC_DIRICHLET) then 
      if(itag == IQ2G) then
        ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,4
          dm%fbcx_gx(i,j,k) = dm%fbcx_qx(i,j,k) * dm%fbcx_ftp(i,j,k)%d
        end do; end do; end do
        !$acc end parallel loop
      end if
      if(itag == IG2Q) then
        ny = npccx(2); nz = npccx(3)
        !$acc parallel loop collapse(3) default(present)
        do k=1,nz; do j=1,ny; do i=1,4
          dm%fbcx_qx(i,j,k) = dm%fbcx_gx(i,j,k) / dm%fbcx_ftp(i,j,k)%d
        end do; end do; end do
        !$acc end parallel loop
      end if
    end if

    if(dm%ibcx_qy(1) == IBC_DIRICHLET .or. dm%ibcx_qy(2) == IBC_DIRICHLET) then
      d_cpc_xpencil(1:ncpcx(1),1:ncpcx(2),1:ncpcx(3))  => fl%wk1
      call transpose_y_to_x(d_cpc_ypencil, d_cpc_xpencil, dm%dcpc)
      if(itag == IQ2G) then
        nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do j=1,ny
          dm%fbcx_gy(1, j, k) = dm%fbcx_qy(1, j, k) * d_cpc_xpencil(1,  j, k)
          dm%fbcx_gy(2, j, k) = dm%fbcx_qy(2, j, k) * d_cpc_xpencil(nx, j, k)
          dm%fbcx_gy(3, j, k) = dm%fbcx_gy(1, j, k)
          dm%fbcx_gy(4, j, k) = dm%fbcx_gy(2, j, k)
        end do; end do
        !$acc end parallel loop
      else  if(itag == IG2Q)then
        nx = ncpcx(1); ny = ncpcx(2); nz = ncpcx(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do j=1,ny
          dm%fbcx_qy(1, j, k) = dm%fbcx_gy(1, j, k) / d_cpc_xpencil(1,  j, k)
          dm%fbcx_qy(2, j, k) = dm%fbcx_gy(2, j, k) / d_cpc_xpencil(nx, j, k)
          dm%fbcx_qy(3, j, k) = dm%fbcx_qy(1, j, k)
          dm%fbcx_qy(4, j, k) = dm%fbcx_qy(2, j, k)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if

    if(dm%ibcx_qz(1) == IBC_DIRICHLET .or. dm%ibcx_qz(2) == IBC_DIRICHLET) then
      d_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))  => fl%wk1
      d_ccp_xpencil(1:nccpx(1),1:nccpx(2),1:nccpx(3))  => fl%wk2
      call transpose_z_to_y(d_ccp_zpencil, d_ccp_ypencil, dm%dccp)
      call transpose_y_to_x(d_ccp_ypencil, d_ccp_xpencil, dm%dccp)
      if(itag == IQ2G) then
        nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do j=1,ny
          dm%fbcx_gz(1, j, k) = dm%fbcx_qz(1, j, k) * d_ccp_xpencil(1,  j, k)
          dm%fbcx_gz(2, j, k) = dm%fbcx_qz(2, j, k) * d_ccp_xpencil(nx, j, k)
          dm%fbcx_gz(3, j, k) = dm%fbcx_gz(1, j, k)
          dm%fbcx_gz(4, j, k) = dm%fbcx_gz(2, j, k)
        end do; end do
        !$acc end parallel loop
      else  if(itag == IG2Q)then
        nx = nccpx(1); ny = nccpx(2); nz = nccpx(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do j=1,ny
          dm%fbcx_qz(1, j, k) = dm%fbcx_gz(1, j, k) / d_ccp_xpencil(1,  j, k)
          dm%fbcx_qz(2, j, k) = dm%fbcx_gz(2, j, k) / d_ccp_xpencil(nx, j, k)
          dm%fbcx_qz(3, j, k) = dm%fbcx_qz(1, j, k)
          dm%fbcx_qz(4, j, k) = dm%fbcx_qz(2, j, k)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if
!----------------------------------------------------------------------------------------------------------
! BC: - y pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qx(1) == IBC_DIRICHLET .or. dm%ibcy_qx(2) == IBC_DIRICHLET) then
      d_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk1
      call transpose_x_to_y(d_pcc_xpencil, d_pcc_ypencil, dm%dpcc)
      if(itag == IQ2G) then
        nx = npccy(1); ny = npccy(2); nz = npccy(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do i=1,nx
          dm%fbcy_gx(i, 1, k) = dm%fbcy_qx(i, 1, k) * d_pcc_ypencil(i,  1, k)
          dm%fbcy_gx(i, 2, k) = dm%fbcy_qx(i, 2, k) * d_pcc_ypencil(i, ny, k)
          dm%fbcy_gx(i, 3, k) = dm%fbcy_gx(i, 1, k)
          dm%fbcy_gx(i, 4, k) = dm%fbcy_gx(i, 2, k)
        end do; end do
        !$acc end parallel loop
      else if(itag == IG2Q) then
        nx = npccy(1); ny = npccy(2); nz = npccy(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do i=1,nx
          dm%fbcy_qx(i, 1, k) = dm%fbcy_gx(i, 1, k) / d_pcc_ypencil(i,  1, k)
          dm%fbcy_qx(i, 2, k) = dm%fbcy_gx(i, 2, k) / d_pcc_ypencil(i, ny, k)
          dm%fbcy_qx(i, 3, k) = dm%fbcy_qx(i, 1, k)
          dm%fbcy_qx(i, 4, k) = dm%fbcy_qx(i, 2, k)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if

    if(dm%ibcy_qy(1) == IBC_DIRICHLET .or. dm%ibcy_qy(2) == IBC_DIRICHLET) then
      if(itag == IQ2G) then
        nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do i=1,nx
          dm%fbcy_gy(i, 1, k) = dm%fbcy_qy(i, 1, k) * d_cpc_ypencil(i,  1, k)
          dm%fbcy_gy(i, 2, k) = dm%fbcy_qy(i, 2, k) * d_cpc_ypencil(i, ny, k)
          dm%fbcy_gy(i, 3, k) = dm%fbcy_gy(i, 1, k)
          dm%fbcy_gy(i, 4, k) = dm%fbcy_gy(i, 2, k)
        end do; end do
        !$acc end parallel loop
      else if(itag == IG2Q) then
        nx = ncpcy(1); ny = ncpcy(2); nz = ncpcy(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do i=1,nx
          dm%fbcy_qy(i, 1, k) = dm%fbcy_gy(i, 1, k) / d_cpc_ypencil(i,  1, k)
          dm%fbcy_qy(i, 2, k) = dm%fbcy_gy(i, 2, k) / d_cpc_ypencil(i, ny, k)
          dm%fbcy_qy(i, 3, k) = dm%fbcy_qy(i, 1, k)
          dm%fbcy_qy(i, 4, k) = dm%fbcy_qy(i, 2, k)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if

    if(dm%ibcy_qz(1) == IBC_DIRICHLET .or. dm%ibcy_qz(2) == IBC_DIRICHLET) then
      d_ccp_ypencil(1:nccpy(1),1:nccpy(2),1:nccpy(3))  => fl%wk1
      call transpose_z_to_y(d_ccp_zpencil, d_ccp_ypencil, dm%dccp)
      if(itag == IQ2G) then
        nx = nccpy(1); ny = nccpy(2); nz = nccpy(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do i=1,nx
          dm%fbcy_gz(i, 1, k) = dm%fbcy_qz(i, 1, k) * d_ccp_ypencil(i,  1, k)
          dm%fbcy_gz(i, 2, k) = dm%fbcy_qz(i, 2, k) * d_ccp_ypencil(i, ny, k)
          dm%fbcy_gz(i, 3, k) = dm%fbcy_gz(i, 1, k)
          dm%fbcy_gz(i, 4, k) = dm%fbcy_gz(i, 2, k)
        end do; end do
        !$acc end parallel loop
      else if(itag == IG2Q) then
        nx = nccpy(1); ny = nccpy(2); nz = nccpy(3)
        !$acc parallel loop collapse(2) default(present)
        do k=1,nz; do i=1,nx
          dm%fbcy_qz(i, 1, k) = dm%fbcy_gz(i, 1, k) / d_ccp_ypencil(i,  1, k)
          dm%fbcy_qz(i, 2, k) = dm%fbcy_gz(i, 2, k) / d_ccp_ypencil(i, ny, k)
          dm%fbcy_qz(i, 3, k) = dm%fbcy_qz(i, 1, k)
          dm%fbcy_qz(i, 4, k) = dm%fbcy_qz(i, 2, k)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if
!----------------------------------------------------------------------------------------------------------
! BC: - z pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcz_qx(1) == IBC_DIRICHLET .or. dm%ibcz_qx(2) == IBC_DIRICHLET) then
      d_pcc_ypencil(1:npccy(1),1:npccy(2),1:npccy(3)) => fl%wk1
      d_pcc_zpencil(1:npccz(1),1:npccz(2),1:npccz(3)) => fl%wk2
      call transpose_x_to_y(d_pcc_xpencil, d_pcc_ypencil)
      call transpose_y_to_z(d_pcc_ypencil, d_pcc_zpencil, dm%dpcc)
      if(itag == IQ2G) then
        nx = npccz(1); ny = npccz(2); nz = npccz(3)
        !$acc parallel loop collapse(2) default(present)
        do j=1,ny; do i=1,nx
          dm%fbcz_gx(i, j, 1) = dm%fbcz_qx(i, j, 1) * d_pcc_zpencil(i, j, 1)
          dm%fbcz_gx(i, j, 2) = dm%fbcz_qx(i, j, 2) * d_pcc_zpencil(i, j, nz)
          dm%fbcz_gx(i, j, 3) = dm%fbcz_gx(i, j, 1)
          dm%fbcz_gx(i, j, 4) = dm%fbcz_gx(i, j, 2)
        end do; end do
        !$acc end parallel loop
      else if(itag == IG2Q) then
        nx = npccz(1); ny = npccz(2); nz = npccz(3)
        !$acc parallel loop collapse(2) default(present)
        do j=1,ny; do i=1,nx
          dm%fbcz_qx(i, j, 1) = dm%fbcz_gx(i, j, 1) / d_pcc_zpencil(i, j, 1)
          dm%fbcz_qx(i, j, 2) = dm%fbcz_gx(i, j, 2) / d_pcc_zpencil(i, j, nz)
          dm%fbcz_qx(i, j, 3) = dm%fbcz_qx(i, j, 1)
          dm%fbcz_qx(i, j, 4) = dm%fbcz_qx(i, j, 2)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if

    if(dm%ibcz_qy(1) == IBC_DIRICHLET .or. dm%ibcz_qy(2) == IBC_DIRICHLET) then
      d_cpc_zpencil(1:ncpcz(1),1:ncpcz(2),1:ncpcz(3))  => fl%wk1
      call transpose_y_to_x(d_cpc_ypencil, d_cpc_zpencil, dm%dcpc)
      if(itag == IQ2G) then
        nx = ncpcz(1); ny = ncpcz(2); nz = ncpcz(3)
        !$acc parallel loop collapse(2) default(present)
        do j=1,ny; do i=1,nx
          dm%fbcz_gy(i, j, 1) = dm%fbcz_qy(i, j, 1) * d_cpc_zpencil(i, j, 1)
          dm%fbcz_gy(i, j, 2) = dm%fbcz_qy(i, j, 2) * d_cpc_zpencil(i, j, nz)
          dm%fbcz_gy(i, j, 3) = dm%fbcz_gy(i, j, 1)
          dm%fbcz_gy(i, j, 4) = dm%fbcz_gy(i, j, 2)
        end do; end do
        !$acc end parallel loop
      else if(itag == IG2Q) then
        nx = ncpcz(1); ny = ncpcz(2); nz = ncpcz(3)
        !$acc parallel loop collapse(2) default(present)
        do j=1,ny; do i=1,nx
          dm%fbcz_qy(i, j, 1) = dm%fbcz_gy(i, j, 1) / d_cpc_zpencil(i, j, 1)
          dm%fbcz_qy(i, j, 2) = dm%fbcz_gy(i, j, 2) / d_cpc_zpencil(i, j, nz)
          dm%fbcz_qy(i, j, 3) = dm%fbcz_qy(i, j, 1)
          dm%fbcz_qy(i, j, 4) = dm%fbcz_qy(i, j, 2)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if

    if(dm%ibcz_qz(1) == IBC_DIRICHLET .or. dm%ibcz_qz(2) == IBC_DIRICHLET) then 
      if(itag == IQ2G) then
        nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
        !$acc parallel loop collapse(2) default(present)
        do j=1,ny; do i=1,nx
          dm%fbcz_gz(i, j, 1) = dm%fbcz_qz(i, j, 1) * d_ccp_zpencil(i, j, 1)
          dm%fbcz_gz(i, j, 2) = dm%fbcz_qz(i, j, 2) * d_ccp_zpencil(i, j, nz)
          dm%fbcz_gz(i, j, 3) = dm%fbcz_gz(i, j, 1)
          dm%fbcz_gz(i, j, 4) = dm%fbcz_gz(i, j, 2)
        end do; end do
        !$acc end parallel loop
      else if(itag == IG2Q) then
        nx = nccpz(1); ny = nccpz(2); nz = nccpz(3)
        !$acc parallel loop collapse(2) default(present)
        do j=1,ny; do i=1,nx
          dm%fbcz_qz(i, j, 1) = dm%fbcz_gz(i, j, 1) / d_ccp_zpencil(i, j, 1)
          dm%fbcz_qz(i, j, 2) = dm%fbcz_gz(i, j, 2) / d_ccp_zpencil(i, j, nz)
          dm%fbcz_qz(i, j, 3) = dm%fbcz_qz(i, j, 1)
          dm%fbcz_qz(i, j, 4) = dm%fbcz_qz(i, j, 2)
        end do; end do
        !$acc end parallel loop
      else
      end if
    end if
    end if

    return
  end subroutine

end module
