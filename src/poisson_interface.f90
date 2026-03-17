
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================================================================================================
module poisson_interface_mod
  use parameters_constant_mod
  use fft2decomp_interface_mod
  use decomp_2d_poisson
  use transpose_extended_mod
  use fishpack_fft
  implicit none

  public :: initialise_fft
  public :: solve_fft_poisson

contains
!==========================================================================================================
!==========================================================================================================
  subroutine initialise_fft(dm)
    use udf_type_mod
    implicit none 
    type(t_domain), intent(in) :: dm

    if(nrank == 0 ) call Print_debug_start_msg("Initialising the Poisson solver ...")
    
    if(dm%ifft_lib == FFT_2DECOMP_3DFFT ) then 
      call build_up_fft2decomp_interface(dm)
      call decomp_2d_poisson_init()
    else if(dm%ifft_lib == FFT_FISHPACK_2DFFT) then 
      call fishpack_fft_init(dm)
    else 
      call Print_error_msg('Error in selecting FFT libs')
    end if

    if(nrank == 0 ) call Print_debug_end_msg()

    return 
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine solve_fft_poisson(rhs_xpencil, fl, dm)
    use udf_type_mod
    use decomp_extended_mod
    implicit none 
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl
    integer :: i, j, k
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(INOUT) :: rhs_xpencil

    real(WP), pointer, dimension(:,:,:) :: rhs_ypencil
    real(WP), pointer, dimension(:,:,:) :: rhs_zpencil
    real(WP), dimension(dm%dccc%zst(1) : dm%dccc%zen(1), &
                        dm%dccc%zst(2) : dm%dccc%zen(2), &
                        dm%dccc%zst(3) : dm%dccc%zen(3)) :: rhs_zpencil_ggg
    integer, dimension(3) :: ncccy, ncccz

    ncccy = dm%dccc%ysz
    ncccz = dm%dccc%zsz
    if(dm%ifft_lib == FFT_2DECOMP_3DFFT ) then
      rhs_ypencil(1:ncccy(1),1:ncccy(2),1:ncccy(3)) => fl%wk1
      rhs_zpencil(1:ncccz(1),1:ncccz(2),1:ncccz(3)) => fl%wk2
      call transpose_x_to_y (rhs_xpencil, rhs_ypencil, dm%dccc)
      call transpose_y_to_z (rhs_ypencil, rhs_zpencil, dm%dccc)
      !$acc data create(rhs_zpencil_ggg)
      call zpencil_index_llg2ggg(rhs_zpencil, rhs_zpencil_ggg, dm%dccc)
      call poisson(rhs_zpencil_ggg, fl)
      call zpencil_index_ggg2llg(rhs_zpencil_ggg, rhs_zpencil, dm%dccc)
      call transpose_z_to_y (rhs_zpencil, rhs_ypencil, dm%dccc)
      call transpose_y_to_x (rhs_ypencil, rhs_xpencil, dm%dccc)
      !$acc end data
    else if(dm%ifft_lib == FFT_FISHPACK_2DFFT) then
#ifdef USE_GPU
      call Print_error_msg('The current fishpack_fft does not support GPUs')
#else
      call fishpack_fft_simple(rhs_xpencil, dm)
#endif
    else
      call Print_error_msg('Error in selecting FFT libs')
    end if

    return
  end subroutine


end module 
