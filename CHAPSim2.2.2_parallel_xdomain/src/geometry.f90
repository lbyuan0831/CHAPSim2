!-------------------------------------------------------------------------------
!                      CHAPSim version 2.0.0
!                      --------------------------
! This file is part of CHAPSim, a general-purpose CFD tool.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 3 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------
!===============================================================================
!> \file geometry.f90
!>
!> \brief Building up the geometry and mesh inwrtfmt1ion.
!>
!===============================================================================
module geometry_mod
  use var_dft_mod, only : domain
  implicit none

  !private
  private :: Buildup_grid_mapping_1D
  private :: Buildup_neibour_index
  public  :: Buildup_geometry_mesh_info
  
contains
!===============================================================================
!===============================================================================
!> \brief Building up the mesh mapping relation between physical domain and mesh
!>  to a computational domain and mesh.   
!>
!> This subroutine is used locally for 1D only.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str          string to indicate mapping of cell centre or nodes
!> \param[in]     n            number of mapping points
!> \param[out]    y            the physical coordinate array
!> \param[out]    mp           the mapping relations for 1st and 2nd deriviatives
!_______________________________________________________________________________
  subroutine Buildup_grid_mapping_1D (str, n, y, mp)
!===============================================================================
! Module files
!===============================================================================
    use math_mod
    use input_general_mod
    use parameters_constant_mod
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    character(len = *), intent(in) :: str
    integer(4), intent( in )       :: n
    real(WP), intent( out )        :: y(n)
    real(WP), intent( out )        :: mp(n, 3)
!===============================================================================
! Local Arguments
!===============================================================================    
    integer(4) :: j
    real(WP) :: eta_shift
    real(WP) :: eta_delta
    real(WP) :: alpha, beta, gamma, delta, cc, dd, ee, st1, st2, mm
    real(WP), dimension(n) :: eta

    eta_shift = ZERO
    eta_delta = ONE
    if ( trim( str ) == 'nd' ) then
      eta_shift = ZERO
      eta_delta = ONE / real( n - 1, WP )
    else if ( trim( str ) == 'cl' ) then
      eta_shift = ONE / ( real(n, WP) ) * HALF
      eta_delta = ONE / real( n, WP )
    else 
      call Print_error_msg('Grid stretching location not defined in Subroutine: '// &
      "Buildup_grid_mapping_1D")
    end if

    ! to build up the computational domain \eta \in [0, 1] uniform mesh
    eta(1) = ZERO + eta_shift

    do j = 2, n
      eta(j) = eta(1) + real(j - 1, WP) * eta_delta
    end do

    ! to build up the physical domain y stretching grids based on Eq(53) of Leizet2009JCP
    ! and to build up the derivates based on Eq(53) and (47) in Leizet2009JCP
    gamma = ONE
    delta = ZERO
    if (istret == ISTRET_NO) then
      y(:) = eta(:)
      y(:) = y(:) * (lyt - lyb) + lyb
      mp(:, 1) = ONE
      mp(:, 2) = ONE
      mp(:, 3) = ONE
      return
    else if (istret == ISTRET_CENTRE) then
      gamma = ONE
      delta = ZERO
    else if (istret == ISTRET_2SIDES) then
      gamma = ONE
      delta = HALF
    else if (istret == ISTRET_BOTTOM) then
      gamma = HALF
      delta = HALF
    else if (istret == ISTRET_TOP) then
      gamma = HALF
      delta = ZERO
    else
      call Print_error_msg('Grid stretching flag is not valid in Subroutine: '// &
      "Buildup_grid_mapping_1D")
    end if

    beta = rstret
    alpha =  ( -ONE + sqrt_wp( ONE + FOUR * PI * PI * beta * beta ) ) / beta * HALF

    cc = sqrt_wp( alpha * beta + ONE ) / sqrt_wp( beta )
    dd = cc / sqrt_wp( alpha )
    ee = cc * sqrt_wp( alpha )

    st1 = (ONE   - TWO * delta) / gamma * HALF
    st2 = (THREE - TWO * delta) / gamma * HALF

    do j = 1, n
      mm = PI * (gamma * eta(j) + delta)

      ! y \in [0, 1]
      y(j) = atan_wp ( dd * tan_wp( mm ) ) - &
            atan_wp ( dd * tan_wp( PI * delta) ) + &
            PI * ( heaviside_step( eta(j) - st1 ) + heaviside_step( eta(j) - st2 ) )
      y(j) = ONE / (gamma * ee) * y(j)
      ! y \in [lyb, lyt]
      y(j) = y(j) * (lyt - lyb) + lyb

      ! 1/h'
      mp(j, 1) = (alpha / PI + sin_wp(mm) * sin_wp(mm) / PI / beta)  / (lyt - lyb)

      ! (1/h')^2
      mp(j, 2) = mp(j, 1) * mp(j, 1)

      ! -h"/(h'^3) = 1/h' * [ d(1/h') / d\eta]
      mp(j, 3) = gamma / (lyt - lyb) / beta * sin_wp(TWO * mm) * mp(j, 1)

    end do

    return
  end subroutine Buildup_grid_mapping_1D
!===============================================================================
!===============================================================================
!> \brief Building up the neibouring index of a given index array.   
!>
!> This subroutine is used locally for the bulk part of the grids. The two points
!> near the boundary are not considered except periodic b.c.
!> The neibouring index reduces the repeated calculation of index increase
!> /decrease for a 5-point stencil. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n            number of the given index range
!> \param[in]     is_peri      whether the given index array is of periodic b.c.
!> \param[out]    nbr          the neibouring index in order of -2, -1, +1, +2
!_______________________________________________________________________________
  subroutine Buildup_npneibour_index(n, ibc, nbr)
    use parameters_constant_mod
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    integer(4), intent(in)  :: n
    integer(4), intent(in)  :: ibc
    integer(4), intent(inout) :: nbr(4, 4)
!===============================================================================
! Local arguments
!===============================================================================
!-------------------------------------------------------------------------------
! nbr(1, i)  nbr(2, i)      nbr(2, i)   nbr(2, i)
!  i -2       i - 1     i    i + 1      i + 2
!-------------------------------------------------------------------------------
    nbr(:, :) = huge(i)

    nbr(1, 3) = 1 + 1
    nbr(1, 4) = 1 + 2
    nbr(2, 2) = 2 - 1
    nbr(2, 3) = 2 + 1 
    nbr(2, 4) = 2 + 2
    if(ibc(1) == IBC_PERIODIC) then
      nbr(1, 1) = n - 1 ! 1'-2' = -1'
      nbr(1, 2) = n     ! 1'-1' = 0'
      nbr(2, 1) = n     ! 2'-2' = 0'
    else if (ibc(1) == IBC_SYMMETRIC .or. ibc(1) == IBC_ASYMMETRIC) then
      nbr(1, 1) = 3 ! 1'-2' = -1'
      nbr(1, 2) = 2     ! 1'-1' = 0'
      nbr(2, 1) = 2     ! 2'-2' = 0'
    else
      ! do nothing, should not be used. 
    end if

    nbr(n, 1) = n - 2
    nbr(n, 2) = n - 1
    nbr(n - 1, 1) = n - 1 - 2
    nbr(n - 1, 2) = n - 1 - 1
    nbr(n - 1, 3) = n - 1 + 1

    if(ibc(2) == IBC_PERIODIC) then 
      nbr(n, 3) = 1   ! n' + 1 
      nbr(n, 4) = 2   ! n' + 2 
      nbr(n - 1, 4) = 1 ! n' - 1 + 2 = n' + 1 
    else if (ibc(2) == IBC_SYMMETRIC .or. ibc(2) == IBC_ASYMMETRIC) then
      nbr(n, 3) = n - 1   ! n' + 1 
      nbr(n, 4) = n - 2   ! n' + 2 
      nbr(n - 1, 4) = n - 1 ! n' - 1 + 2 = n' + 1 
    else
      ! do nothing, should not be used. 
    end if

    return
  end subroutine

!_______________________________________________________________________________
  subroutine Buildup_ncneibour_index(n, ibc, nbr)
    use parameters_constant_mod
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    integer(4), intent(in)  :: n
    integer(4), intent(in)  :: ibc
    integer(4), intent(inout) :: nbr(4, 4)
!===============================================================================
! Local arguments
!===============================================================================
!-------------------------------------------------------------------------------
! nbr(1, i)  nbr(2, i)      nbr(2, i)   nbr(2, i)
!  i -2       i - 1     i    i + 1      i + 2
!-------------------------------------------------------------------------------
    nbr(:, :) = huge(i)

    nbr(1, 3) = 1 + 1
    nbr(1, 4) = 1 + 2
    nbr(2, 2) = 2 - 1
    nbr(2, 3) = 2 + 1 
    nbr(2, 4) = 2 + 2
    if(ibc(1) == IBC_PERIODIC) then
      nbr(1, 1) = n - 1 ! 1-2 = -1
      nbr(1, 2) = n     ! 1-1 = 0
      nbr(2, 1) = n     ! 2-2 = 0
    else if (ibc(1) == IBC_SYMMETRIC .or. ibc(1) == IBC_ASYMMETRIC) then
      nbr(1, 1) = 2     ! 1-2 = -1
      nbr(1, 2) = 1     ! 1-1 = 0
      nbr(2, 1) = 1     ! 2-2 = 0
    else
      ! do nothing, should not be used. 
    end if

    nbr(n, 1) = n - 2
    nbr(n, 2) = n - 1
    nbr(n - 1, 1) = n - 1 - 2
    nbr(n - 1, 2) = n - 1 - 1
    nbr(n - 1, 3) = n - 1 + 1

    if(ibc(2) == IBC_PERIODIC) then 
      nbr(n, 3) = 1   ! n + 1 
      nbr(n, 4) = 2   ! n + 2 
      nbr(n - 1, 4) = 1 ! n - 1 + 2 = n + 1 
    else if (ibc(2) == IBC_SYMMETRIC .or. ibc(2) == IBC_ASYMMETRIC) then
      nbr(n, 3) = n   ! n + 1 
      nbr(n, 4) = n - 1   ! n + 2 
      nbr(n - 1, 4) = n ! n - 1 + 2 = n + 1 
    else
      ! do nothing, should not be used. 
    end if

    return
  end subroutine
!===============================================================================
!===============================================================================
  subroutine Buildup_geometry_mesh_info (dm)
    use mpi_mod
    use input_general_mod
    use math_mod
    use parameters_constant_mod, only : ONE, HALF, ZERO, MAXP, MINP, TRUNCERR
    implicit none
    type(t_domain), intent(inout) :: dm
    integer(4) :: i
    logical    :: dbg = .false.

    if(nrank == 0) call Print_debug_start_msg("Initializing domain geometric inwrtfmt1ion ...")
    ! Build up domain info

    dm%is_periodic(:) == .false.
    if(dm%ibcx(1, 1) == IBC_PERIODIC) dm%is_periodic(1) = .true.
    if(dm%ibcy(1, 1) == IBC_PERIODIC) dm%is_periodic(2) = .true.
    if(dm%ibcz(1, 1) == IBC_PERIODIC) dm%is_periodic(3) = .true.

    dm%np_geo(1) = dm%nc(1) + 1 
    dm%np_geo(2) = dm%nc(2) + 1
    dm%np_geo(3) = dm%nc(3) + 1

    do i = 1, 3
      if ( dm%is_periodic(i) ) then
        dm%np(i) = dm%nc(i)
      else 
        dm%np(i) = dm%np_geo(i)
      end if
    end do

    dm%is_stretching(:) = .false.
    if (dm%istret /= ISTRET_NO) dm%is_stretching(2) = .true.
    
    if(dm%is_stretching(2)) then
      dm%h(2) = ONE / real(dm%nc(2), WP)
    else 
      dm%h(2) = (dm%lyt - dm%lyb) / real(dm%nc(2), WP) ! mean dy
    end if
    dm%h(1) = dm%lxx / real(dm%nc(1), WP)
    dm%h(3) = dm%lzz / real(dm%nc(3), WP)
    dm%h2r(:) = ONE / dm%h(:) / dm%h(:)
    dm%h1r(:) = ONE / dm%h(:)

    !build up index sequence for boundary part

    call Buildup_npneibour_index (dm%np(1), dm%ibcx(1,1:2), dm%ipnbr(:, :) )
    call Buildup_npneibour_index (dm%np(2), dm%ibcy(1,1:2), dm%jpnbr(:, :) )
    call Buildup_npneibour_index (dm%np(3), dm%ibcz(1,1:2), dm%kpnbr(:, :) )

    call Buildup_ncneibour_index (dm%nc(1), dm%ibcx(1,1:2), dm%icnbr(:, :) )
    call Buildup_ncneibour_index (dm%nc(2), dm%ibcy(1,1:2), dm%jcnbr(:, :) )
    call Buildup_ncneibour_index (dm%nc(3), dm%ibcz(1,1:2), dm%kcnbr(:, :) )

    ! allocate  variables for mapping physical domain to computational domain
    allocate ( dm%yp( dm%np_geo(2) ) ); dm%yp(:) = ZERO
    allocate ( dm%yc( dm%nc(2) ) ); dm%yc(:) = ZERO

    allocate ( dm%yMappingpt( dm%np_geo(2), 3 ) ); dm%yMappingpt(:, :) = ONE
    allocate ( dm%yMappingcc( dm%nc(2),     3 ) ); dm%yMappingcc(:, :) = ONE

    call Buildup_grid_mapping_1D ('nd', dm%np_geo(2), dm%yp(:), dm%yMappingPt(:, :))
    call Buildup_grid_mapping_1D ('cl', dm%nc(2),     dm%yc(:), dm%yMappingcc(:, :))

    ! print out for debugging
    if(dbg) then
      do i = 1, dm%np_geo(2)
        write (OUTPUT_UNIT, '(I5, 1F8.4)') i, dm%yp(i)
      end do
    end if

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine  Buildup_geometry_mesh_info
end module geometry_mod
