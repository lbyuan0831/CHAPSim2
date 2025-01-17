module statistics_mod

  character(13), parameter :: io_name = "statistics-io"
  integer, allocatable :: ncl_stat(:, :)

  private :: write_statistics_array
  private :: read_statistics_array

  public  :: init_statistics_flow
  public  :: update_statistics_flow
  public  :: write_statistics_flow

  public  :: init_statistics_thermo
  public  :: update_statistics_thermo
  public  :: write_statistics_thermo

contains
!==========================================================================================================
!==========================================================================================================
  subroutine init_statistics_flow(fl, dm)
    use udf_type_mod
    use parameters_constant_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    !integer :: i

    allocate (ncl_stat(3, nxdomain))
    ncl_stat = 0

    ! do i = 1, 3
    !   if(dm%is_periodic(i)) then 
    !     ncl_stat(i, dm%idom) = xszS(i)
    !   else 
    !     ncl_stat(i, dm%idom) = MAX(xszS(i) - 1, 1)
    !   end if
    ! end do
    ncl_stat(1, dm%idom) = dm%dccc%xsz(1) ! default skip is 1.
    ncl_stat(2, dm%idom) = dm%dccc%xsz(2) ! default skip is 1.
    ncl_stat(3, dm%idom) = dm%dccc%xsz(3) ! default skip is 1.


    allocate ( fl%pr_mean        (ncl_stat(1, dm%idom), ncl_stat(2, dm%idom), ncl_stat(3, dm%idom)   ) )
    allocate ( fl%u_vector_mean  (ncl_stat(1, dm%idom), ncl_stat(2, dm%idom), ncl_stat(3, dm%idom), 3) )
    allocate ( fl%uu_tensor6_mean(ncl_stat(1, dm%idom), ncl_stat(2, dm%idom), ncl_stat(3, dm%idom), 6) )
    fl%u_vector_mean   = ZERO
    fl%pr_mean         = ZERO
    fl%uu_tensor6_mean = ZERO

    if(fl%inittype == INIT_RESTART .and. fl%iteration > dm%stat_istart) then
      call read_statistics_array(fl%pr_mean,                     'mean_pr', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%u_vector_mean  (:, :, :, 1), 'mean_ux', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%u_vector_mean  (:, :, :, 2), 'mean_uy', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%u_vector_mean  (:, :, :, 3), 'mean_uz', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%uu_tensor6_mean(:, :, :, 1), 'mean_uu', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%uu_tensor6_mean(:, :, :, 2), 'mean_vv', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%uu_tensor6_mean(:, :, :, 3), 'mean_ww', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%uu_tensor6_mean(:, :, :, 4), 'mean_uv', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%uu_tensor6_mean(:, :, :, 5), 'mean_uw', dm%idom, fl%iterfrom, dm%dccc)
      call read_statistics_array(fl%uu_tensor6_mean(:, :, :, 6), 'mean_vw', dm%idom, fl%iterfrom, dm%dccc)
    end if

    return
  end subroutine
  !==========================================================================================================
!==========================================================================================================
  subroutine init_statistics_thermo(tm, dm)
    use udf_type_mod
    use parameters_constant_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(inout) :: tm

    if(.not. dm%is_thermo) return

    allocate ( tm%t_mean  (ncl_stat(1, dm%idom), ncl_stat(2, dm%idom), ncl_stat(3, dm%idom)) )
    allocate ( tm%tt_mean (ncl_stat(1, dm%idom), ncl_stat(2, dm%idom), ncl_stat(3, dm%idom)) )
    
    tm%t_mean  = ZERO
    tm%tt_mean = ZERO
    if(tm%inittype == INIT_RESTART .and. tm%iteration > dm%stat_istart) then
      call read_statistics_array(tm%t_mean,  'mean_t',  dm%idom, tm%iterfrom, dm%dccc)
      call read_statistics_array(tm%tt_mean, 'mean_tt', dm%idom, tm%iterfrom, dm%dccc)
    end if

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================

  subroutine update_statistics_flow(fl, dm)
    use udf_type_mod
    use parameters_constant_mod
    use operations
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc1
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc2
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc3
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil

    real(WP) :: ac, am
    integer :: nstat

! to do: coarse mesh does not work right now due to the above allocation.
    
!----------------------------------------------------------------------------------------------------------
!   this is for a asymptotic averaging ... 
!----------------------------------------------------------------------------------------------------------
    nstat = fl%iteration - dm%stat_istart + 1
    ac = ONE / real(nstat, WP)
    am = real(nstat - 1, WP) / real(nstat, WP)
!----------------------------------------------------------------------------------------------------------
!   pressure, stored in cell centre
!----------------------------------------------------------------------------------------------------------
    fl%pr_mean(:, :, :) = am * fl%pr_mean(:, :, :) + ac * fl%pres(:, :, :)
!----------------------------------------------------------------------------------------------------------
!   ux
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(fl%qx, accc1, dm, dm%ibcx(:, 1), dm%fbcx_var(:, :, :, 1) )
    fl%u_vector_mean(:, :, :, 1) = am * fl%u_vector_mean(:, :, :, 1) + ac * accc1(:, :, :)
!----------------------------------------------------------------------------------------------------------
!   uy
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy_var(:, :, :, 2))
    call transpose_y_to_x(accc_ypencil, accc2, dm%dccc)
    fl%u_vector_mean(:, :, :, 2) = am * fl%u_vector_mean(:, :, :, 2) + ac * accc2(:, :, :)
!----------------------------------------------------------------------------------------------------------
!   uz
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qz, accp_ypencil, dm%dccp)
    call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
    call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_var(:, :, :, 3) )
    call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
    call transpose_y_to_x(accc_ypencil, accc3, dm%dccc)
    fl%u_vector_mean(:, :, :, 3) = am * fl%u_vector_mean(:, :, :, 3) + ac * accc3(:, :, :)
!----------------------------------------------------------------------------------------------------------
!   tensor, uu, vv, ww, uv, uw, vw, x-pencil, stored in cell centre
!----------------------------------------------------------------------------------------------------------
    fl%uu_tensor6_mean(:, :, :, 1) = am * fl%uu_tensor6_mean(:, :, :, 1) + ac * accc1(:, :, :) * accc1(:, :, :)
    fl%uu_tensor6_mean(:, :, :, 2) = am * fl%uu_tensor6_mean(:, :, :, 2) + ac * accc2(:, :, :) * accc2(:, :, :)
    fl%uu_tensor6_mean(:, :, :, 3) = am * fl%uu_tensor6_mean(:, :, :, 3) + ac * accc3(:, :, :) * accc3(:, :, :)
    fl%uu_tensor6_mean(:, :, :, 4) = am * fl%uu_tensor6_mean(:, :, :, 4) + ac * accc1(:, :, :) * accc2(:, :, :)
    fl%uu_tensor6_mean(:, :, :, 5) = am * fl%uu_tensor6_mean(:, :, :, 5) + ac * accc1(:, :, :) * accc3(:, :, :)
    fl%uu_tensor6_mean(:, :, :, 6) = am * fl%uu_tensor6_mean(:, :, :, 6) + ac * accc2(:, :, :) * accc3(:, :, :)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_statistics_flow(fl, dm)
    use udf_type_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl

    call write_statistics_array(fl%pr_mean,                     'mean_pr', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%u_vector_mean  (:, :, :, 1), 'mean_ux', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%u_vector_mean  (:, :, :, 2), 'mean_uy', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%u_vector_mean  (:, :, :, 3), 'mean_uz', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%uu_tensor6_mean(:, :, :, 1), 'mean_uu', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%uu_tensor6_mean(:, :, :, 2), 'mean_vv', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%uu_tensor6_mean(:, :, :, 3), 'mean_ww', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%uu_tensor6_mean(:, :, :, 4), 'mean_uv', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%uu_tensor6_mean(:, :, :, 5), 'mean_uw', dm%idom, fl%iteration, dm%dccc)
    call write_statistics_array(fl%uu_tensor6_mean(:, :, :, 6), 'mean_vw', dm%idom, fl%iteration, dm%dccc)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine update_statistics_thermo(tm, dm)
    use udf_type_mod
    use parameters_constant_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(inout) :: tm
    real(WP) :: ac, am
    integer :: nstat

! to do: coarse mesh does not work right now due to the above allocation.
    
!----------------------------------------------------------------------------------------------------------
!   this is for a asymptotic averaging ... 
!----------------------------------------------------------------------------------------------------------
    nstat = tm%iteration - dm%stat_istart + 1
    ac = ONE / real(nstat, WP)
    am = real(nstat - 1, WP) / real(nstat, WP)
!----------------------------------------------------------------------------------------------------------
!   temperature
!----------------------------------------------------------------------------------------------------------
    tm%t_mean (:, :, :) = am * tm%t_mean(:, :, :)  + ac * tm%tTemp(:, :, :)
    tm%tt_mean(:, :, :) = am * tm%tt_mean(:, :, :) + ac * tm%tTemp(:, :, :) * tm%tTemp(:, :, :)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_statistics_thermo(tm, dm)
    use udf_type_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    call write_statistics_array(tm%t_mean,  'mean_t', dm%idom, tm%iteration, dm%dccc)
    call write_statistics_array(tm%tt_mean, 'mean_tt', dm%idom, tm%iteration, dm%dccc)


    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_statistics_array(var, keyword, idm, iter, dtmp)
    use udf_type_mod
    use files_io_mod
    use io_tools_mod
    use decomp_2d_io
    implicit none 
    real(WP), intent(in) :: var( :, :, :)
    type(DECOMP_INFO), intent(in) :: dtmp
    character(*), intent(in) :: keyword
    integer, intent(in) :: idm
    integer, intent(in) :: iter

    character(120):: data_flname_path

    call generate_pathfile_name(data_flname_path, idm, trim(keyword), dir_data, 'bin', iter)
    call decomp_2d_write_one(X_PENCIL, var, trim(data_flname_path), dtmp)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine read_statistics_array(var, keyword, idm, iter, dtmp)
    use precision_mod
    !use files_io_mod, only: dir_data
    use io_tools_mod
    use decomp_2d_io
    implicit none 
    
    real(WP), dimension(:, :, :), intent(in) :: var
    character(*), intent(in) :: keyword
    integer, intent(in) :: idm
    integer, intent(in) :: iter
    type(DECOMP_INFO), intent(in) :: dtmp

    character(120):: data_flname

    call generate_file_name(data_flname, idm, trim(keyword), 'bin', iter)
    !call decomp_2d_read_one(X_PENCIL, var, trim(dir_data), trim(data_flname), io_name, dtmp, reduce_prec=.false.)
    !to do, check why the above line does not work.

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================


end module