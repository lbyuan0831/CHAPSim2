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
! This program is disatributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------
!===============================================================================
!> \file operations.f90
!>
!> \brief A general operation of derivative and interpolation in 1D.
!>
!===============================================================================
module operations
  use precision_mod, only : WP
  implicit none

  private
!-------------------------------------------------------------------------------
! basic coefficients for TDMA of 1st deriviative  
! to store coefficients for TDMA
! eg, d1fC2C(5, 3, 5)
!     First column: 1:2 for one side b.c.
!                   4:5 for the other side b.c.
!                   3   for interior
!     Second column: 1 for coefficients of f^(1)_{i-1}
!                    2 for coefficients of f^(1)_{i}
!                    3 for coefficients of f^(1)_{i+1}
!     Third column:  for b.c. flags
!                 integer, parameter :: IBC_INTERIOR    = 0, &
!                 IBC_PERIODIC    = 1, &
!                 IBC_SYMMETRIC   = 2, &
!                 IBC_ASYMMETRIC  = 3, &
!                 IBC_DIRICHLET   = 4, &
!                 IBC_NEUMANN     = 5, &
!                 all others      = 6
!     d1fC2C vs d1rC2C :
!       f : coefficients in the LHS, unknown side.
!       r : coefficients in the RHS, known side. 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! for 1st derivative
!-------------------------------------------------------------------------------
  ! collocated C2C
  real(WP), save, public :: d1fC2C(5, 3, 6)
  real(WP), save, public :: d1rC2C(5, 4, 6)
  
  ! collocated P2P
  real(WP) :: d1fP2P(5, 3, 6)
  real(WP) :: d1rP2P(5, 4, 6)

  ! staggered C2P
  real(WP) :: d1fC2P(5, 3, 6)
  real(WP) :: d1rC2P(5, 4, 6)

  ! staggered P2C
  real(WP) :: d1fP2C(5, 3, 6)
  real(WP) :: d1rP2C(5, 4, 6)
!-------------------------------------------------------------------------------
! for 2nd derivative
!-------------------------------------------------------------------------------
  ! collocated C2C
  real(WP), public :: d2fC2C(5, 3, 6)
  real(WP), public :: d2rC2C(5, 4, 6) ! one more value used. 
  
  ! collocated P2P
  real(WP) :: d2fP2P(5, 3, 6)
  real(WP) :: d2rP2P(5, 4, 6)
!-------------------------------------------------------------------------------
! for iterpolation
!-------------------------------------------------------------------------------
  ! interpolation P2C
  real(WP) :: m1fP2C(5, 3, 6)
  real(WP) :: m1rP2C(5, 4, 6)

  ! interpolation C2P
  real(WP) :: m1fC2P(5, 3, 6)
  real(WP) :: m1rC2P(5, 4, 6)

!-------------------------------------------------------------------------------
! coefficients array for TDMA of 1st deriviative  
! to store coefficients array for TDMA
!-------------------------------------------------------------------------------
  type t_xtdma_lhs
!-------------------------------------------------------------------------------
!   x : pre-processed TDMA LHS Matrix for 1st deriviative
!-------------------------------------------------------------------------------
    real(WP), allocatable :: ad1x_P2P(:, :, :)
    real(WP), allocatable :: bd1x_P2P(:, :, :)
    real(WP), allocatable :: cd1x_P2P(:, :, :)
    real(WP), allocatable :: dd1x_P2P(:, :, :)

    real(WP), allocatable :: ad1x_C2C(:, :, :)
    real(WP), allocatable :: bd1x_C2C(:, :, :)
    real(WP), allocatable :: cd1x_C2C(:, :, :)
    real(WP), allocatable :: dd1x_C2C(:, :, :)

    real(WP), allocatable :: ad1x_P2C(:, :, :)
    real(WP), allocatable :: bd1x_P2C(:, :, :)
    real(WP), allocatable :: cd1x_P2C(:, :, :)
    real(WP), allocatable :: dd1x_P2C(:, :, :)

    real(WP), allocatable :: ad1x_C2P(:, :, :)
    real(WP), allocatable :: bd1x_C2P(:, :, :)
    real(WP), allocatable :: cd1x_C2P(:, :, :)
    real(WP), allocatable :: dd1x_C2P(:, :, :)
!-------------------------------------------------------------------------------
!   x : pre-processed TDMA LHS Matrix for 2nd deriviative
!-------------------------------------------------------------------------------
    real(WP), allocatable :: ad2x_P2P(:, :, :)
    real(WP), allocatable :: bd2x_P2P(:, :, :)
    real(WP), allocatable :: cd2x_P2P(:, :, :)
    real(WP), allocatable :: dd2x_P2P(:, :, :)

    real(WP), allocatable :: ad2x_C2C(:, :, :)
    real(WP), allocatable :: bd2x_C2C(:, :, :)
    real(WP), allocatable :: cd2x_C2C(:, :, :)
    real(WP), allocatable :: dd2x_C2C(:, :, :)
!-------------------------------------------------------------------------------
!   x : pre-processed TDMA LHS Matrix for mid-point interpolation
!-------------------------------------------------------------------------------
    real(WP), allocatable :: am1x_P2C(:, :, :)
    real(WP), allocatable :: bm1x_P2C(:, :, :)
    real(WP), allocatable :: cm1x_P2C(:, :, :)
    real(WP), allocatable :: dm1x_P2C(:, :, :)

    real(WP), allocatable :: am1x_C2P(:, :, :)
    real(WP), allocatable :: bm1x_C2P(:, :, :)
    real(WP), allocatable :: cm1x_C2P(:, :, :)
    real(WP), allocatable :: dm1x_C2P(:, :, :)
  end type t_xtdma_lhs

  type(t_xtdma_lhs), allocatable :: xtdma_lhs(:) 

!-------------------------------------------------------------------------------
! y : pre-processed TDMA LHS Matrix for 1st deriviative
!-------------------------------------------------------------------------------
  real(WP), allocatable :: ad1y_P2P(:, :, :)
  real(WP), allocatable :: bd1y_P2P(:, :, :)
  real(WP), allocatable :: cd1y_P2P(:, :, :)
  real(WP), allocatable :: dd1y_P2P(:, :, :)

  real(WP), allocatable :: ad1y_C2C(:, :, :)
  real(WP), allocatable :: bd1y_C2C(:, :, :)
  real(WP), allocatable :: cd1y_C2C(:, :, :)
  real(WP), allocatable :: dd1y_C2C(:, :, :)

  real(WP), allocatable :: ad1y_P2C(:, :, :)
  real(WP), allocatable :: bd1y_P2C(:, :, :)
  real(WP), allocatable :: cd1y_P2C(:, :, :)
  real(WP), allocatable :: dd1y_P2C(:, :, :)

  real(WP), allocatable :: ad1y_C2P(:, :, :)
  real(WP), allocatable :: bd1y_C2P(:, :, :)
  real(WP), allocatable :: cd1y_C2P(:, :, :)
  real(WP), allocatable :: dd1y_C2P(:, :, :)
!-------------------------------------------------------------------------------
! y : pre-processed TDMA LHS Matrix for 2nd deriviative
!-------------------------------------------------------------------------------
  real(WP), allocatable :: ad2y_P2P(:, :, :)
  real(WP), allocatable :: bd2y_P2P(:, :, :)
  real(WP), allocatable :: cd2y_P2P(:, :, :)
  real(WP), allocatable :: dd2y_P2P(:, :, :)

  real(WP), allocatable :: ad2y_C2C(:, :, :)
  real(WP), allocatable :: bd2y_C2C(:, :, :)
  real(WP), allocatable :: cd2y_C2C(:, :, :)
  real(WP), allocatable :: dd2y_C2C(:, :, :)
!-------------------------------------------------------------------------------
! y : pre-processed TDMA LHS Matrix for mid-point interpolation
!-------------------------------------------------------------------------------
  real(WP), allocatable :: am1y_P2C(:, :, :)
  real(WP), allocatable :: bm1y_P2C(:, :, :)
  real(WP), allocatable :: cm1y_P2C(:, :, :)
  real(WP), allocatable :: dm1y_P2C(:, :, :)

  real(WP), allocatable :: am1y_C2P(:, :, :)
  real(WP), allocatable :: bm1y_C2P(:, :, :)
  real(WP), allocatable :: cm1y_C2P(:, :, :)
  real(WP), allocatable :: dm1y_C2P(:, :, :)

!-------------------------------------------------------------------------------
! z : pre-processed TDMA LHS Matrix for 1st deriviative
!-------------------------------------------------------------------------------
  real(WP), allocatable :: ad1z_P2P(:, :, :)
  real(WP), allocatable :: bd1z_P2P(:, :, :)
  real(WP), allocatable :: cd1z_P2P(:, :, :)
  real(WP), allocatable :: dd1z_P2P(:, :, :)

  real(WP), allocatable :: ad1z_C2C(:, :, :)
  real(WP), allocatable :: bd1z_C2C(:, :, :)
  real(WP), allocatable :: cd1z_C2C(:, :, :)
  real(WP), allocatable :: dd1z_C2C(:, :, :)

  real(WP), allocatable :: ad1z_P2C(:, :, :)
  real(WP), allocatable :: bd1z_P2C(:, :, :)
  real(WP), allocatable :: cd1z_P2C(:, :, :)
  real(WP), allocatable :: dd1z_P2C(:, :, :)

  real(WP), allocatable :: ad1z_C2P(:, :, :)
  real(WP), allocatable :: bd1z_C2P(:, :, :)
  real(WP), allocatable :: cd1z_C2P(:, :, :)
  real(WP), allocatable :: dd1z_C2P(:, :, :)
!-------------------------------------------------------------------------------
! z : pre-processed TDMA LHS Matrix for 2nd deriviative
!-------------------------------------------------------------------------------
  real(WP), allocatable :: ad2z_P2P(:, :, :)
  real(WP), allocatable :: bd2z_P2P(:, :, :)
  real(WP), allocatable :: cd2z_P2P(:, :, :)
  real(WP), allocatable :: dd2z_P2P(:, :, :)

  real(WP), allocatable :: ad2z_C2C(:, :, :)
  real(WP), allocatable :: bd2z_C2C(:, :, :)
  real(WP), allocatable :: cd2z_C2C(:, :, :)
  real(WP), allocatable :: dd2z_C2C(:, :, :)
!-------------------------------------------------------------------------------
! z : pre-processed TDMA LHS Matrix for mid-point interpolation
!-------------------------------------------------------------------------------
  real(WP), allocatable :: am1z_P2C(:, :, :)
  real(WP), allocatable :: bm1z_P2C(:, :, :)
  real(WP), allocatable :: cm1z_P2C(:, :, :)
  real(WP), allocatable :: dm1z_P2C(:, :, :)

  real(WP), allocatable :: am1z_C2P(:, :, :)
  real(WP), allocatable :: bm1z_C2P(:, :, :)
  real(WP), allocatable :: cm1z_C2P(:, :, :)
  real(WP), allocatable :: dm1z_C2P(:, :, :)
  
!-------------------------------------------------------------------------------
! processures
!-------------------------------------------------------------------------------
  private :: Prepare_compact_coefficients
  private :: Buildup_TDMA_LHS_array

  private :: Prepare_TDMA_interp_P2C_RHS_array
  private :: Prepare_TDMA_interp_C2P_RHS_array ! need fbc for Dirichlet

  private :: Prepare_TDMA_1deri_C2C_RHS_array
  private :: Prepare_TDMA_1deri_P2P_RHS_array ! need fbc for Neumann
  private :: Prepare_TDMA_1deri_C2P_RHS_array ! need fbc for Neumann
  private :: Prepare_TDMA_1deri_P2C_RHS_array

  private :: Prepare_TDMA_2deri_C2C_RHS_array
  private :: Prepare_TDMA_2deri_P2P_RHS_array ! need fbc for Neumann

  public  :: Prepare_LHS_coeffs_for_operations

  public  :: Get_x_midp_C2P_1D
  public  :: Get_y_midp_C2P_1D
  public  :: Get_z_midp_C2P_1D
  public  :: Get_x_midp_P2C_1D
  public  :: Get_y_midp_P2C_1D 
  public  :: Get_z_midp_P2C_1D

  public  :: Get_x_midp_C2P_3D
  public  :: Get_y_midp_C2P_3D
  public  :: Get_z_midp_C2P_3D
  public  :: Get_x_midp_P2C_3D
  public  :: Get_y_midp_P2C_3D
  public  :: Get_z_midp_P2C_3D

  

  public  :: Get_x_1st_derivative_P2P_1D
  public  :: Get_y_1st_derivative_P2P_1D
  public  :: Get_z_1st_derivative_P2P_1D

  public  :: Get_x_1st_derivative_C2P_1D
  public  :: Get_y_1st_derivative_C2P_1D
  public  :: Get_z_1st_derivative_C2P_1D

  public  :: Get_x_1st_derivative_C2C_1D
  public  :: Get_y_1st_derivative_C2C_1D
  public  :: Get_z_1st_derivative_C2C_1D

  public  :: Get_x_1st_derivative_P2C_1D
  public  :: Get_y_1st_derivative_P2C_1D
  public  :: Get_z_1st_derivative_P2C_1D

  public  :: Get_x_1st_derivative_P2P_3D
  public  :: Get_y_1st_derivative_P2P_3D
  public  :: Get_z_1st_derivative_P2P_3D

  public  :: Get_x_1st_derivative_C2P_3D
  public  :: Get_y_1st_derivative_C2P_3D
  public  :: Get_z_1st_derivative_C2P_3D

  public  :: Get_x_1st_derivative_C2C_3D
  public  :: Get_y_1st_derivative_C2C_3D
  public  :: Get_z_1st_derivative_C2C_3D
  
  public  :: Get_x_1st_derivative_P2C_3D
  public  :: Get_y_1st_derivative_P2C_3D
  public  :: Get_z_1st_derivative_P2C_3D

  public  :: Get_x_2nd_derivative_C2C_1D
  public  :: Get_y_2nd_derivative_C2C_1D
  public  :: Get_z_2nd_derivative_C2C_1D

  public  :: Get_x_2nd_derivative_P2P_1D
  public  :: Get_y_2nd_derivative_P2P_1D
  public  :: Get_z_2nd_derivative_P2P_1D

  public  :: Get_x_2nd_derivative_C2C_3D
  public  :: Get_y_2nd_derivative_C2C_3D
  public  :: Get_z_2nd_derivative_C2C_3D

  public  :: Get_x_2nd_derivative_P2P_3D
  public  :: Get_y_2nd_derivative_P2P_3D
  public  :: Get_z_2nd_derivative_P2P_3D

contains
!===============================================================================
!> \brief Assigned the cooefficients for the compact schemes     
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!>
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           !
!-------------------------------------------------------------------------------
!> \param[in]     iaccu         the accuracy given by user
!===============================================================================
  subroutine Prepare_compact_coefficients(iaccu)
    use parameters_constant_mod
    use input_general_mod
    use mpi_mod
    implicit none

    integer, intent(in) :: iaccu

    real(WP) :: alpha,  a,  b,  c,  d
    real(WP) :: alpha1, a1, b1, c1, d1
    real(WP) :: alpha2, a2, b2, c2, d2

    if(nrank == 0) then
       call Print_debug_start_msg &
         ("Assigning coefficient matrix for the compact FD ...")
       write(*, *) "The given numerical accuracy =", iaccu
    end if

!-------------------------------------------------------------------------------
!   initialisation
!-------------------------------------------------------------------------------
    d1fC2C(:, :, :) = ZERO
    d1rC2C(:, :, :) = ZERO
    d1fP2P(:, :, :) = ZERO
    d1rP2P(:, :, :) = ZERO

    d1fC2P(:, :, :) = ZERO
    d1rC2P(:, :, :) = ZERO
    d1fP2C(:, :, :) = ZERO
    d1rP2C(:, :, :) = ZERO

    d2fC2C(:, :, :) = ZERO
    d2rC2C(:, :, :) = ZERO
    d2fP2P(:, :, :) = ZERO
    d2rP2P(:, :, :) = ZERO

    m1fC2P(:, :, :) = ZERO
    m1rC2P(:, :, :) = ZERO
    m1fP2C(:, :, :) = ZERO
    m1rP2C(:, :, :) = ZERO
!===============================================================================
! Set 1 : P2P, C2P, periodic & symmetric & asymmetric
!         1st derivative on collocated grids, C2C/P2P bulk coefficients
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!===============================================================================
    alpha = ZERO
        a = ZERO
        b = ZERO
        c = ZERO
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
          a = ONE
          b = ZERO
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
          a = FOUR / THREE
          b = - ONE / THREE
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / FOUR
          a = THREE / TWO
          b = ZERO
    else if (iaccu == IACCU_CP6) then
      alpha = ONE / THREE
          a = FOURTEEN / NINE
          b = ONE / NINE
    else ! default 2nd CD
      alpha = ZERO
          a = ONE
          b = ZERO
    end if
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2C : periodic b.c.
! d1fC2C : "d1"=first deriviative, "f"=f'  side, "C2C"= center 2 centre 
! d1rC2C : "d1"=first deriviative, "r"=rhs side, "C2C"= center 2 centre 
! [ 1    alpha                   alpha][f'_1]=[a/2 * (f_{2}   - f_{n})/h   + b/4 * (f_{3}   - f_{n-1})/h]
! [      alpha 1     alpha            ][f'_2] [a/2 * (f_{3}   - f_{1})/h   + b/4 * (f_{4}   - f_{n})/h  ]
! [            alpha 1     alpha      ][f'_i] [a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                  alpha 1     alpha][f'_4] [a/2 * (f_{n}   - f_{n-2})/h + b/4 * (f_{1}   - f_{n-3})/h]
! [alpha                   alpha 1    ][f'_5] [a/2 * (f_{1}   - f_{n-1})/h + b/4 * (f_{2}   - f_{n-2})/h]
!-------------------------------------------------------------------------------
    d1fC2C(1:5, 1, IBC_PERIODIC) = alpha
    d1fC2C(1:5, 2, IBC_PERIODIC) = ONE
    d1fC2C(1:5, 3, IBC_PERIODIC) = alpha

    d1rC2C(1:5, 1, IBC_PERIODIC) = a / TWO  ! a/2
    d1rC2C(1:5, 2, IBC_PERIODIC) = b / FOUR ! b/4
    d1rC2C(1:5, 3, IBC_PERIODIC) = c        ! not used
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2P : periodic b.c.  Same as C2C
!-------------------------------------------------------------------------------
    d1fP2P(:, :, IBC_PERIODIC) = d1fC2C(:, :, IBC_PERIODIC)
    d1rP2P(:, :, IBC_PERIODIC) = d1rC2C(:, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2C : symmetric b.c.
! [ 1-alpha  alpha                          ][f'_1]=[a/2 * (f_{2}   - f_{1})/h   + b/4 * (f_{3}   - f_{2})/h  ]
! [          alpha 1     alpha              ][f'_2] [a/2 * (f_{3}   - f_{1})/h   + b/4 * (f_{4}   - f_{1})/h  ]
! [                alpha 1     alpha        ][f'_i] [a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                      alpha 1     alpha  ][f'_4] [a/2 * (f_{n}   - f_{n-2})/h + b/4 * (f_{n}   - f_{n-3})/h]
! [                            alpha 1-alpha][f'_5] [a/2 * (f_{n}   - f_{n-1})/h + b/4 * (f_{n-1} - f_{n-2})/h]
!-------------------------------------------------------------------------------
    d1fC2C(1,   1, IBC_SYMMETRIC) = ZERO        ! not used
    d1fC2C(1,   2, IBC_SYMMETRIC) = ONE - alpha
    d1fC2C(1,   3, IBC_SYMMETRIC) = alpha
    d1fC2C(5,   1, IBC_SYMMETRIC) = d1fC2C(1,   3, IBC_SYMMETRIC)
    d1fC2C(5,   2, IBC_SYMMETRIC) = d1fC2C(1,   2, IBC_SYMMETRIC)
    d1fC2C(5,   3, IBC_SYMMETRIC) = d1fC2C(1,   1, IBC_SYMMETRIC)

    d1fC2C(2:4, :, IBC_SYMMETRIC) = d1fC2C(2:4, :, IBC_PERIODIC )

    d1rC2C(:,   :, IBC_SYMMETRIC) = d1rC2C(:,   :, IBC_PERIODIC )
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2C : asymmetric b.c.
! [ 1+alpha  alpha                          ][f'_1]=[a/2 * (f_{2}   + f_{1})/h   + b/4 * (f_{3}   + f_{2})/h  ]
! [          alpha 1     alpha              ][f'_2] [a/2 * (f_{3}   - f_{1})/h   + b/4 * (f_{4}   + f_{1})/h  ]
! [                alpha 1     alpha        ][f'_i] [a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                      alpha 1     alpha  ][f'_4] [a/2 * (f_{n}   - f_{n-2})/h + b/4 * (-f_{n}   - f_{n-3})/h]
! [                            alpha 1+alpha][f'_5] [a/2 * (-f_{n}   - f_{n-1})/h + b/4 * (-f_{n-1} - f_{n-2})/h]
!-------------------------------------------------------------------------------
    d1fC2C(1,   1, IBC_ASYMMETRIC) = ZERO        ! not used
    d1fC2C(1,   2, IBC_ASYMMETRIC) = ONE + alpha
    d1fC2C(1,   3, IBC_ASYMMETRIC) = alpha

    d1fC2C(5,   1, IBC_ASYMMETRIC) = d1fC2C(1,   3, IBC_ASYMMETRIC)
    d1fC2C(5,   2, IBC_ASYMMETRIC) = d1fC2C(1,   2, IBC_ASYMMETRIC)
    d1fC2C(5,   3, IBC_ASYMMETRIC) = d1fC2C(1,   1, IBC_ASYMMETRIC)

    d1fC2C(2:4, :, IBC_ASYMMETRIC) = d1fC2C(2:4, :, IBC_PERIODIC  )

    d1rC2C(:,   :, IBC_ASYMMETRIC) = d1rC2C(:,   :, IBC_PERIODIC  )
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2P : symmetric b.c.
! [ 1  0                              ][f'_{1'}]=[a/2 * (f_{2'}   - f_{2'})/h   + b/4 * (f_{3'}   - f_{3'})/h  ]
! [    alpha 1     alpha              ][f'_{2'}] [a/2 * (f_{3'}   - f_{1'})/h   + b/4 * (f_{4'}   - f_{2'})/h  ]
! [          alpha 1     alpha        ][f'_{i'}] [a/2 * (f_{i'+1} - f_{i'-1})/h + b/4 * (f_{i'+2} - f_{i'-2})/h]
! [                alpha 1     alpha  ][f'_{4'}] [a/2 * (f_{n'}   - f_{n'-2})/h + b/4 * (f_{n'-1} - f_{n'-3})/h]
! [                      0     1      ][f'_{5'}] [a/2 * (f_{n'-1} - f_{n'-1})/h + b/4 * (f_{n'-2} - f_{n'-2})/h]
!-------------------------------------------------------------------------------
    d1fP2P(1,   1, IBC_SYMMETRIC) = ZERO ! not used
    d1fP2P(1,   2, IBC_SYMMETRIC) = ONE
    d1fP2P(1,   3, IBC_SYMMETRIC) = ZERO

    d1fP2P(5,   1, IBC_SYMMETRIC) = d1fP2P(1,   3, IBC_SYMMETRIC)
    d1fP2P(5,   2, IBC_SYMMETRIC) = d1fP2P(1,   2, IBC_SYMMETRIC)
    d1fP2P(5,   3, IBC_SYMMETRIC) = d1fP2P(1,   1, IBC_SYMMETRIC)

    d1fP2P(2:4, :, IBC_SYMMETRIC) = d1fP2P(2:4, :, IBC_PERIODIC )

    d1rP2P(:,   :, IBC_SYMMETRIC) = d1rP2P(:,   :, IBC_PERIODIC )
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2P : asymmetric b.c.
! [ 1  2alpha                          ][f'_{1'}]=[a/2 * (f_{2'}   + f_{2'})/h   + b/4 * (f_{3'}   + f_{3'})/h  ]
! [    alpha 1     alpha               ][f'_{2'}] [a/2 * (f_{3'}   - f_{1'})/h   + b/4 * (f_{4'}   + f_{2'})/h  ]
! [          alpha 1     alpha         ][f'_{i'}] [a/2 * (f_{i'+1} - f_{i'-1})/h + b/4 * (f_{i'+2} - f_{i'-2})/h]
! [                alpha 1      alpha  ][f'_{4'}] [a/2 * (f_{n'}   - f_{n'-2})/h + b/4 * (-f_{n'-1} - f_{n'-3})/h]
! [                      2alpha 1      ][f'_{5'}] [a/2 * (-f_{n'-1} - f_{n'-1})/h + b/4 * (-f_{n'-2} - f_{n'-2})/h]
!-------------------------------------------------------------------------------
    d1fP2P(1,   1, IBC_ASYMMETRIC) = ZERO ! not used
    d1fP2P(1,   2, IBC_ASYMMETRIC) = ONE
    d1fP2P(1,   3, IBC_ASYMMETRIC) = TWO * alpha

    d1fP2P(5,   1, IBC_ASYMMETRIC) = d1fP2P(1,   3, IBC_ASYMMETRIC)
    d1fP2P(5,   2, IBC_ASYMMETRIC) = d1fP2P(1,   2, IBC_ASYMMETRIC)
    d1fP2P(5,   3, IBC_ASYMMETRIC) = d1fP2P(1,   1, IBC_ASYMMETRIC)

    d1fP2P(2:4, :, IBC_ASYMMETRIC) = d1fP2P(2:4, :, IBC_PERIODIC)

    d1rP2P(:,   :, IBC_ASYMMETRIC) = d1rP2P(:,   :, IBC_PERIODIC)
!===============================================================================
! Set 2: no bc required
!       C2C : no specified, Neumann
!       P2P : no specified, Dirichlet B.C.
! 1st derivative on collocated grids, C2C/P2P coefficients : Dirichlet B.C.
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!===============================================================================
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
    alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then ! degrade to 2nd CD
      alpha1 = ZERO
          a1 = -THREE / TWO
          b1 = TWO
          c1 = -ONE / TWO
      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6) then ! degrade to 3rd CP
      alpha1 = TWO
          a1 = -FIVE / TWO
          b1 = TWO
          c1 = ONE / TWO
      alpha2 = ONE / FOUR
          a2 = THREE / TWO
          b2 = ZERO
    else ! default 2nd CD
      alpha1 = ZERO
          a1 = -THREE / TWO
          b1 = TWO
          c1 = -ONE / TWO
     alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
    end if
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2C : no specified, Neumann
! P2P : no specified, Dirichlet B.C.
! [ 1     alpha1                            ][f'_1]=[a1 * f_{1}/h  + b1 * f_{2}/h + c1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2] [a2/2 * (f_{3} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i] [ a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4] [a2/2 * (f_{n} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5] [-a1 * f_{n}/h  - b1 * f_{n-1}/h - c1 * f_{n-2}/h]
!-------------------------------------------------------------------------------
    d1fC2C(1, 1,   IBC_INTRPL) = ZERO ! not used
    d1fC2C(1, 2,   IBC_INTRPL) = ONE
    d1fC2C(1, 3,   IBC_INTRPL) = alpha1
    d1rC2C(1, 1,   IBC_INTRPL) = a1
    d1rC2C(1, 2,   IBC_INTRPL) = b1
    d1rC2C(1, 3,   IBC_INTRPL) = c1

    d1fC2C(5, 1,   IBC_INTRPL) =   d1fC2C(1, 3, IBC_INTRPL)
    d1fC2C(5, 2,   IBC_INTRPL) =   d1fC2C(1, 2, IBC_INTRPL)
    d1fC2C(5, 3,   IBC_INTRPL) =   d1fC2C(1, 1, IBC_INTRPL)
    d1rC2C(5, 1,   IBC_INTRPL) = - d1rC2C(1, 1, IBC_INTRPL)
    d1rC2C(5, 2,   IBC_INTRPL) = - d1rC2C(1, 2, IBC_INTRPL)
    d1rC2C(5, 3,   IBC_INTRPL) = - d1rC2C(1, 3, IBC_INTRPL)
  
    d1fC2C(2, 1,   IBC_INTRPL) = alpha2
    d1fC2C(2, 2,   IBC_INTRPL) = ONE
    d1fC2C(2, 3,   IBC_INTRPL) = alpha2
    d1rC2C(2, 1,   IBC_INTRPL) = a2 / TWO
    d1rC2C(2, 2,   IBC_INTRPL) = b2 / FOUR
    d1rC2C(2, 3,   IBC_INTRPL) = c2  ! not used  

    d1fC2C(4, 1,   IBC_INTRPL) = d1fC2C(2, 3, IBC_INTRPL)
    d1fC2C(4, 2,   IBC_INTRPL) = d1fC2C(2, 2, IBC_INTRPL)
    d1fC2C(4, 3,   IBC_INTRPL) = d1fC2C(2, 1, IBC_INTRPL)
    d1rC2C(4, 1,   IBC_INTRPL) = d1rC2C(2, 1, IBC_INTRPL)
    d1rC2C(4, 2,   IBC_INTRPL) = d1rC2C(2, 2, IBC_INTRPL)
    d1rC2C(4, 3,   IBC_INTRPL) = d1rC2C(2, 3, IBC_INTRPL)

    d1fC2C(3, 1:3, IBC_INTRPL) = d1fC2C(3, 1:3, IBC_PERIODIC)
    d1rC2C(3, 1:3, IBC_INTRPL) = d1rC2C(3, 1:3, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2C : neumann
!-------------------------------------------------------------------------------
    d1fC2C(:, :, IBC_NEUMANN) = d1fC2C(:, :, IBC_INTRPL)
    d1rC2C(:, :, IBC_NEUMANN) = d1rC2C(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2P : no specified
!-------------------------------------------------------------------------------
    d1fP2P(:, :, IBC_INTRPL) = d1fC2C(:, :, IBC_INTRPL)
    d1rP2P(:, :, IBC_INTRPL) = d1rC2C(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2P : Dirichlet
!-------------------------------------------------------------------------------
    d1fP2P(:, :, IBC_DIRICHLET) = d1fP2P(:, :, IBC_INTRPL)
    d1rP2P(:, :, IBC_DIRICHLET) = d1rP2P(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2P : NEUMANN
! [ 1     0                                 ][f'_1]=[known]
! [alpha2 1      alpha2                     ][f'_2] [a2/2 * (f_{3} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i] [ a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4] [a2/2 * (f_{n} - f_{n-2})/h]
! [                            0      1     ][f'_5] [known]
!-------------------------------------------------------------------------------
    d1fP2P(1, 1,   IBC_NEUMANN) = ZERO ! not used
    d1fP2P(1, 2,   IBC_NEUMANN) = ONE
    d1fP2P(1, 3,   IBC_NEUMANN) = ZERO
    d1rP2P(1, 1,   IBC_NEUMANN) = ZERO
    d1rP2P(1, 2,   IBC_NEUMANN) = ZERO
    d1rP2P(1, 3,   IBC_NEUMANN) = ZERO

    d1fP2P(5, 1,   IBC_NEUMANN) =   d1fP2P(1, 3, IBC_NEUMANN)
    d1fP2P(5, 2,   IBC_NEUMANN) =   d1fP2P(1, 2, IBC_NEUMANN)
    d1fP2P(5, 3,   IBC_NEUMANN) =   d1fP2P(1, 1, IBC_NEUMANN)
    d1rP2P(5, 1,   IBC_NEUMANN) = - d1rP2P(1, 1, IBC_NEUMANN)
    d1rP2P(5, 2,   IBC_NEUMANN) = - d1rP2P(1, 2, IBC_NEUMANN)
    d1rP2P(5, 3,   IBC_NEUMANN) = - d1rP2P(1, 3, IBC_NEUMANN)

    d1fP2P(2, 1,   IBC_NEUMANN) = alpha2
    d1fP2P(2, 2,   IBC_NEUMANN) = ONE
    d1fP2P(2, 3,   IBC_NEUMANN) = alpha2
    d1rP2P(2, 1,   IBC_NEUMANN) = a2 / TWO
    d1rP2P(2, 2,   IBC_NEUMANN) = b2 / FOUR ! not used
    d1rP2P(2, 3,   IBC_NEUMANN) = c2    

    d1fP2P(4, 1,   IBC_NEUMANN) = d1fP2P(2, 3, IBC_NEUMANN)
    d1fP2P(4, 2,   IBC_NEUMANN) = d1fP2P(2, 2, IBC_NEUMANN)
    d1fP2P(4, 3,   IBC_NEUMANN) = d1fP2P(2, 1, IBC_NEUMANN)
    d1rP2P(4, 1,   IBC_NEUMANN) = d1rP2P(2, 1, IBC_NEUMANN)
    d1rP2P(4, 2,   IBC_NEUMANN) = d1rP2P(2, 2, IBC_NEUMANN)
    d1rP2P(4, 3,   IBC_NEUMANN) = d1rP2P(2, 3, IBC_NEUMANN)

    d1fP2P(3, 1:3, IBC_NEUMANN) = d1fP2P(3, 1:3, IBC_PERIODIC)
    d1rP2P(3, 1:3, IBC_NEUMANN) = d1rP2P(3, 1:3, IBC_PERIODIC)
!===============================================================================
! Set 3: Dirichlet for C2C (unique)
! 1st derivative on collocated grids, C2C/P2P coefficients : Dirichlet B.C.
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!===============================================================================
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
        d1 = ZERO
    alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then ! degrade to 2nd CD
      alpha1 = ZERO
          a1 = -FOUR / THREE
          b1 = ONE
          c1 = ONE / THREE
          d1 = ZERO
      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP
      alpha1 = ONE / THREE
          a1 = - EIGHT / NINE
          b1 = ZERO
          c1 = EIGHT / NINE
          d1 = ZERO
      alpha2 = ONE / FOUR
          a2 = THREE / TWO
          b2 = ZERO
    else if (iaccu == IACCU_CP6) then ! degrade to 4th CP
      alpha1 = TWO / THREE
          a1 = - THIRTYTWO / FOURTYFIVE
          b1 = - SEVEN / FIVE
          c1 = TEN / NINE
          d1 = ONE / TEN
      alpha2 = ONE / FOUR
          a2 = THREE / TWO
          b2 = ZERO
    else ! default 2nd CD
      alpha1 = ZERO
          a1 = -FOUR / THREE
          b1 = ONE
          c1 = ONE / THREE
          d1 = ZERO
      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
    end if
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2C : Dirchlet
! [ 1     alpha1                            ][f'_1]=[a1 * f_{1'}/h  + b1 * f_{1}/h + c1 * f_{2}/h + d1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2] [a2/2 * (f_{3} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i] [ a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4] [a2/2 * (f_{n} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5] [-a1 * f_{n'+1}/h -b1 * f_{n}/h  - c1 * f_{n-1}/h - d1 * f_{n-2}/h]
!-------------------------------------------------------------------------------
    d1fC2C(1, 1,   IBC_DIRICHLET) = ZERO ! not used
    d1fC2C(1, 2,   IBC_DIRICHLET) = ONE
    d1fC2C(1, 3,   IBC_DIRICHLET) = alpha1
    d1rC2C(1, 1,   IBC_DIRICHLET) = a1
    d1rC2C(1, 2,   IBC_DIRICHLET) = b1
    d1rC2C(1, 3,   IBC_DIRICHLET) = c1
    d1rC2C(1, 4,   IBC_DIRICHLET) = d1

    d1fC2C(5, 1,   IBC_DIRICHLET) =   d1fC2C(1, 3, IBC_DIRICHLET)
    d1fC2C(5, 2,   IBC_DIRICHLET) =   d1fC2C(1, 2, IBC_DIRICHLET)
    d1fC2C(5, 3,   IBC_DIRICHLET) =   d1fC2C(1, 1, IBC_DIRICHLET)
    d1rC2C(5, 1,   IBC_DIRICHLET) = - d1rC2C(1, 1, IBC_DIRICHLET)
    d1rC2C(5, 2,   IBC_DIRICHLET) = - d1rC2C(1, 2, IBC_DIRICHLET)
    d1rC2C(5, 3,   IBC_DIRICHLET) = - d1rC2C(1, 3, IBC_DIRICHLET)
    d1rC2C(5, 4,   IBC_DIRICHLET) = - d1rC2C(1, 4, IBC_DIRICHLET)

    d1fC2C(2, 1,   IBC_DIRICHLET) = alpha2
    d1fC2C(2, 2,   IBC_DIRICHLET) = ONE
    d1fC2C(2, 3,   IBC_DIRICHLET) = alpha2
    d1rC2C(2, 1,   IBC_DIRICHLET) = a2 / TWO
    d1rC2C(2, 2,   IBC_DIRICHLET) = b2 / FOUR ! not used
    d1rC2C(2, 3,   IBC_DIRICHLET) = c2    

    d1fC2C(4, 1,   IBC_DIRICHLET) = d1fC2C(2, 3, IBC_DIRICHLET)
    d1fC2C(4, 2,   IBC_DIRICHLET) = d1fC2C(2, 2, IBC_DIRICHLET)
    d1fC2C(4, 3,   IBC_DIRICHLET) = d1fC2C(2, 1, IBC_DIRICHLET)
    d1rC2C(4, 1,   IBC_DIRICHLET) = d1rC2C(2, 1, IBC_DIRICHLET)
    d1rC2C(4, 2,   IBC_DIRICHLET) = d1rC2C(2, 2, IBC_DIRICHLET)
    d1rC2C(4, 3,   IBC_DIRICHLET) = d1rC2C(2, 3, IBC_DIRICHLET)

    d1fC2C(3, 1:3, IBC_DIRICHLET) = d1fC2C(3, 1:3, IBC_PERIODIC)
    d1rC2C(3, 1:3, IBC_DIRICHLET) = d1rC2C(3, 1:3, IBC_PERIODIC)
!===============================================================================
! 1st derivative on staggered grids P2C and C2P : Periodic or Symmetric B.C.
! P2C ==>
! alpha * f'_{i-1} +  f'_i +  alpha * f'_{i+1}  = a/(h ) * ( f_{i'+1} - f_{i'} ) + &
!                                                 b/(3h) * ( f_{i'+2} - f_{i'-1} )
! C2P ==>
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(h ) * ( f_{i}   - f_{i-1} ) + &
!                                                 b/(3h) * ( f_{i+1} - f_{i-2} )
!===============================================================================
    alpha = ZERO
        a = ZERO
        b = ZERO
        c = ZERO
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
          a = ONE
          b = ZERO
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
          a = NINE / EIGHT
          b = -ONE / EIGHT
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / TWENTYTWO
          a = TWELVE / ELEVEN
          b = ZERO
    else if (iaccu == IACCU_CP6) then
      alpha = NINE / SIXTYTWO
          a = SIXTYTHREE / SIXTYTWO
          b = SEVENTEEN / SIXTYTWO
    else  ! default 2nd CD
      alpha = ZERO
          a = ONE
          b = ZERO
    end if
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2P, P2C: periodic b.c. : staggered 
! [ 1    alpha                   alpha][f'_1']=[a * (f_{1}   - f_{n})/h   + b/3 * (f_{2}   - f_{n-1})/h]
! [      alpha 1     alpha            ][f'_2'] [a * (f_{1}   - f_{1})/h   + b/3 * (f_{3}   - f_{n})/h  ]
! [            alpha 1     alpha      ][f'_i'] [a * (f_{i}   - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h]
! [                  alpha 1     alpha][f'_4'] [a * (f_{n-1} - f_{n-2})/h + b/3 * (f_{n}   - f_{n-3})/h]
! [alpha                   alpha 1    ][f'_5'] [a * (f_{n}   - f_{n-1})/h + b/3 * (f_{1}   - f_{n-2})/h]
!-------------------------------------------------------------------------------
    d1fC2P(1:5, 1, IBC_PERIODIC) = alpha
    d1fC2P(1:5, 2, IBC_PERIODIC) = ONE
    d1fC2P(1:5, 3, IBC_PERIODIC) = alpha
    d1rC2P(1:5, 1, IBC_PERIODIC) = a         ! a
    d1rC2P(1:5, 2, IBC_PERIODIC) = b / THREE ! b/3
    d1rC2P(1:5, 3, IBC_PERIODIC) = c         ! not used
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2C : periodic b.c.  Same as C2P
!-------------------------------------------------------------------------------
    d1fP2C(:, :, IBC_PERIODIC) = d1fC2P(:, :, IBC_PERIODIC)
    d1rP2C(:, :, IBC_PERIODIC) = d1rC2P(:, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2P : symmetric b.c.
! [ 1     0                        ][f'_1]=[a * (f_{1}   - f_{1})/h   + b/3 * (f_{2}   - f_{2})/h   ]
! [ alpha 1     alpha              ][f'_2] [a * (f_{2}   - f_{1})/h   + b/3 * (f_{3}   - f_{1})/h   ]
! [       alpha 1     alpha        ][f'_i] [a * (f_{i}   - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h ]
! [             alpha 1     alpha  ][f'_4] [a * (f_{n-1} - f_{n-2})/h + b/3 * (f_{n-1} - f_{n-3})/h ]
! [                   0     1      ][f'_5] [a * (f_{n-1} - f_{n-1})/h + b/3 * (f_{n-2} - f_{n-2})/h ]
!-------------------------------------------------------------------------------
    d1fC2P(1,   1, IBC_SYMMETRIC) = ZERO ! not used
    d1fC2P(1,   2, IBC_SYMMETRIC) = ONE
    d1fC2P(1,   3, IBC_SYMMETRIC) = ZERO

    d1fC2P(5,   1, IBC_SYMMETRIC) = d1fC2P(1,   3, IBC_SYMMETRIC)
    d1fC2P(5,   2, IBC_SYMMETRIC) = d1fC2P(1,   2, IBC_SYMMETRIC)
    d1fC2P(5,   3, IBC_SYMMETRIC) = d1fC2P(1,   1, IBC_SYMMETRIC)

    d1fC2P(2:4, :, IBC_SYMMETRIC) = d1fC2P(2:4, :, IBC_PERIODIC)

    d1rC2P(:,   :, IBC_SYMMETRIC) = d1rC2P(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2P : asymmetric b.c.
! [ 1     2alpha                   ][f'_1]=[a * (f_{1}   + f_{1})/h   + b/3 * (f_{2}   + f_{2})/h   ]
! [ alpha 1     alpha              ][f'_2] [a * (f_{2}   - f_{1})/h   + b/3 * (f_{3}   + f_{1})/h   ]
! [       alpha 1     alpha        ][f'_i] [a * (f_{i}   - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h ]
! [             alpha 1     alpha  ][f'_4] [a * (f_{n-1} - f_{n-2})/h + b/3 * (-f_{n-1} - f_{n-3})/h ]
! [                   2alpha     1 ][f'_5] [a * (-f_{n-1} - f_{n-1})/h + b/3 * (-f_{n-2} - f_{n-2})/h ]
!-------------------------------------------------------------------------------
    d1fC2P(1,   1, IBC_ASYMMETRIC) = ZERO ! not used
    d1fC2P(1,   2, IBC_ASYMMETRIC) = ONE
    d1fC2P(1,   3, IBC_ASYMMETRIC) = TWO * alpha

    d1fC2P(5,   1, IBC_ASYMMETRIC) = d1fC2P(1,   3, IBC_ASYMMETRIC)
    d1fC2P(5,   2, IBC_ASYMMETRIC) = d1fC2P(1,   2, IBC_ASYMMETRIC)
    d1fC2P(5,   3, IBC_ASYMMETRIC) = d1fC2P(1,   1, IBC_ASYMMETRIC)

    d1fC2P(2:4, :, IBC_ASYMMETRIC) = d1fC2P(2:4, :, IBC_PERIODIC)

    d1rC2P(:,   :, IBC_ASYMMETRIC) = d1rC2P(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2C : symmetric b.c.
! [ 1-alpha  alpha                          ][f_1]=[a * (f_{2'}   - f_{1'})/h   + b/3 * (f_{3'}   - f_{2'})/h    ]
! [          alpha 1     alpha              ][f_2] [a * (f_{3'}   - f_{2'})/h   + b/3 * (f_{4'}   - f_{1'})/h    ]
! [                alpha 1     alpha        ][f_i] [a * (f_{i+1}  - f_{i-1})/h  + b/3 * (f_{i+2}  - f_{i-2})/h   ]
! [                      alpha 1     alpha  ][f_4] [a * (f_{n'}   - f_{n'-1})/h + b/3 * (f_{n'+1} - f_{n'-2})/h  ]
! [                            alpha 1-alpha][f_5] [a * (f_{n'+1} - f_{n'})/h   + b/3 * (f_{n'}   - f_{n'-1'})/h ]
!-------------------------------------------------------------------------------
    d1fP2C(1,   1, IBC_SYMMETRIC) = ZERO ! not used
    d1fP2C(1,   2, IBC_SYMMETRIC) = ONE - alpha
    d1fP2C(1,   3, IBC_SYMMETRIC) = alpha

    d1fP2C(5,   1, IBC_SYMMETRIC) = d1fP2C(1,   3, IBC_SYMMETRIC)
    d1fP2C(5,   2, IBC_SYMMETRIC) = d1fP2C(1,   2, IBC_SYMMETRIC)
    d1fP2C(5,   3, IBC_SYMMETRIC) = d1fP2C(1,   1, IBC_SYMMETRIC)

    d1fP2C(2:4, :, IBC_SYMMETRIC) = d1fP2C(2:4, :, IBC_PERIODIC)

    d1rP2C(:,   :, IBC_SYMMETRIC) = d1rP2C(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2C : asymmetric b.c.
! [ 1+alpha  alpha                          ][f_1]=[a * (f_{2'}   - f_{1'})/h   + b/3 * (f_{3'}   + f_{2'})/h    ]
! [          alpha 1     alpha              ][f_2] [a * (f_{3'}   - f_{2'})/h   + b/3 * (f_{4'}   - f_{1'})/h    ]
! [                alpha 1     alpha        ][f_i] [a * (f_{i+1}  - f_{i-1})/h  + b/3 * (f_{i+2}  - f_{i-2})/h   ]
! [                      alpha 1     alpha  ][f_4] [a * (f_{n'}   - f_{n'-1})/h + b/3 * (f_{n'+1} - f_{n'-2})/h  ]
! [                            alpha 1+alpha][f_5] [a * (f_{n'+1} - f_{n'})/h   + b/3 * (-f_{n'}   - f_{n'-1'})/h ]
!-------------------------------------------------------------------------------
    d1fP2C(1,   1, IBC_ASYMMETRIC) = ZERO ! not used
    d1fP2C(1,   2, IBC_ASYMMETRIC) = ONE + alpha
    d1fP2C(1,   3, IBC_ASYMMETRIC) = alpha

    d1fP2C(5,   1, IBC_ASYMMETRIC) = d1fP2C(1,   3, IBC_ASYMMETRIC)
    d1fP2C(5,   2, IBC_ASYMMETRIC) = d1fP2C(1,   2, IBC_ASYMMETRIC)
    d1fP2C(5,   3, IBC_ASYMMETRIC) = d1fP2C(1,   1, IBC_ASYMMETRIC)

    d1fP2C(2:4, :, IBC_ASYMMETRIC) = d1fP2C(2:4, :, IBC_PERIODIC)

    d1rP2C(:,   :, IBC_ASYMMETRIC) = d1rP2C(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2P : no specified, interpolation
! [ 1     alpha1                            ][f'_1']=[a1 * f_{1}/h  + b1 * f_{2}/h + c1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2'] [a2 * (f_{2} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i'] [a *  (f_{i} - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4'] [a2 * (f_{n-1} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5'] [-a1 * f_{n-1}/h  - b1 * f_{n-2}/h - c1 * f_{n-3}/h]
!-------------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
    alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
          a1 = -TWO
          b1 = THREE
          c1 = -ONE

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = TWENTYTHREE
          a1 = -TWENTYFIVE
          b1 = TWENTYSIX
          c1 = -ONE

      alpha2 = ONE / TWENTYTWO
          a2 = TWELVE / ELEVEN
          b2 = ZERO ! not used
          c2 = ZERO ! not used
      
    else  ! default 2nd CD
     alpha1 = ZERO
          a1 = -TWO
          b1 = THREE
          c1 = -ONE

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    end if

    d1fC2P(1, 1, IBC_INTRPL) = ZERO ! not used
    d1fC2P(1, 2, IBC_INTRPL) = ONE
    d1fC2P(1, 3, IBC_INTRPL) = alpha1
    d1rC2P(1, 1, IBC_INTRPL) = a1
    d1rC2P(1, 2, IBC_INTRPL) = b1
    d1rC2P(1, 3, IBC_INTRPL) = c1

    d1fC2P(5, 1, IBC_INTRPL) = d1fC2P(1, 3, IBC_INTRPL)
    d1fC2P(5, 2, IBC_INTRPL) = d1fC2P(1, 2, IBC_INTRPL)
    d1fC2P(5, 3, IBC_INTRPL) = d1fC2P(1, 1, IBC_INTRPL)
    d1rC2P(5, 1, IBC_INTRPL) = - d1rC2P(1, 1, IBC_INTRPL)
    d1rC2P(5, 2, IBC_INTRPL) = - d1rC2P(1, 2, IBC_INTRPL)
    d1rC2P(5, 3, IBC_INTRPL) = - d1rC2P(1, 3, IBC_INTRPL)

    d1fC2P(2, 1, IBC_INTRPL) = alpha2
    d1fC2P(2, 2, IBC_INTRPL) = ONE
    d1fC2P(2, 3, IBC_INTRPL) = alpha2
    d1rC2P(2, 1, IBC_INTRPL) = a2
    d1rC2P(2, 2, IBC_INTRPL) = b2 / THREE ! not used
    d1rC2P(2, 3, IBC_INTRPL) = c2 ! not used

    d1fC2P(4, 1, IBC_INTRPL) = d1fC2P(2, 1, IBC_INTRPL)
    d1fC2P(4, 2, IBC_INTRPL) = d1fC2P(2, 2, IBC_INTRPL)
    d1fC2P(4, 3, IBC_INTRPL) = d1fC2P(2, 3, IBC_INTRPL)
    d1rC2P(4, 1, IBC_INTRPL) = d1rC2P(2, 1, IBC_INTRPL)
    d1rC2P(4, 2, IBC_INTRPL) = d1rC2P(2, 2, IBC_INTRPL)
    d1rC2P(4, 3, IBC_INTRPL) = d1rC2P(2, 3, IBC_INTRPL)

    d1fC2P(3, :, IBC_INTRPL) = d1fC2P(3, :, IBC_PERIODIC)
    d1rC2P(3, :, IBC_INTRPL) = d1rC2P(3, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2P : neumann
! [ 1     0                                 ][f'_1']=[known]
! [alpha2 1      alpha2                     ][f'_2'] [a2 * (f_{2} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i'] [a *  (f_{i} - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4'] [a2 * (f_{n-1} - f_{n-2})/h]
! [                            0      1     ][f'_5'] [known]
!-------------------------------------------------------------------------------
    d1fC2P(1, 1, IBC_NEUMANN) = ZERO ! not used
    d1fC2P(1, 2, IBC_NEUMANN) = ONE
    d1fC2P(1, 3, IBC_NEUMANN) = ZERO
    d1rC2P(1, 1, IBC_NEUMANN) = ZERO ! not used
    d1rC2P(1, 2, IBC_NEUMANN) = ZERO ! not used
    d1rC2P(1, 3, IBC_NEUMANN) = ZERO ! not used

    d1fC2P(5, 1, IBC_NEUMANN) = d1fC2P(1, 3, IBC_NEUMANN)
    d1fC2P(5, 2, IBC_NEUMANN) = d1fC2P(1, 2, IBC_NEUMANN)
    d1fC2P(5, 3, IBC_NEUMANN) = d1fC2P(1, 1, IBC_NEUMANN)
    d1rC2P(5, 1, IBC_NEUMANN) = d1rC2P(1, 1, IBC_NEUMANN)
    d1rC2P(5, 2, IBC_NEUMANN) = d1rC2P(1, 2, IBC_NEUMANN)
    d1rC2P(5, 3, IBC_NEUMANN) = d1rC2P(1, 3, IBC_NEUMANN)

    d1fC2P(2, 1, IBC_NEUMANN) = alpha2
    d1fC2P(2, 2, IBC_NEUMANN) = ONE
    d1fC2P(2, 3, IBC_NEUMANN) = alpha2
    d1rC2P(2, 1, IBC_NEUMANN) = a2
    d1rC2P(2, 2, IBC_NEUMANN) = b2 / THREE ! not used
    d1rC2P(2, 3, IBC_NEUMANN) = c2 ! not used

    d1fC2P(4, 1, IBC_NEUMANN) = d1fC2P(2, 1, IBC_NEUMANN)
    d1fC2P(4, 2, IBC_NEUMANN) = d1fC2P(2, 2, IBC_NEUMANN)
    d1fC2P(4, 3, IBC_NEUMANN) = d1fC2P(2, 3, IBC_NEUMANN)
    d1rC2P(4, 1, IBC_NEUMANN) = d1rC2P(2, 1, IBC_NEUMANN)
    d1rC2P(4, 2, IBC_NEUMANN) = d1rC2P(2, 2, IBC_NEUMANN)
    d1rC2P(4, 3, IBC_NEUMANN) = d1rC2P(2, 3, IBC_NEUMANN)

    d1fC2P(3, :, IBC_NEUMANN) = d1fC2P(3, :, IBC_PERIODIC)
    d1rC2P(3, :, IBC_NEUMANN) = d1rC2P(3, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! C2P : Dirichlet
! [ 1     alpha1                            ][f'_1']=[a1 * f_{1'}/h + b1 * f_{1}/h  + c1 * f_{2}/h + d1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2'] [a2 * (f_{2} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i'] [a *  (f_{i} - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4'] [a2 * (f_{n-1} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5'] [-a1 * f_{n'}/h - b1 * f_{n-1}/h  - c1 * f_{n-2}/h - d1 * f_{n-3}/h]
!-------------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
        d1 = ZERO
    alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
          a1 = - EIGHT / THREE
          b1 = THREE
          c1 = - ONE / THREE
          d1 = ZERO
      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP

      alpha1 = THREE
          a1 = - EIGHT / THREE
          b1 = ZERO
          c1 = EIGHT / THREE
          d1 = ZERO
      alpha2 = ONE / TWENTYTWO
          a2 = TWELVE / ELEVEN
          b2 = ZERO ! not used
          c2 = ZERO ! not used
    else if (iaccu == IACCU_CP6) then ! degrade to 4th CP

      alpha1 = FIFTEEN
          a1 = - SIXTEEN / FIFTEEN
          b1 = - FIFTEEN
          c1 = FIFTY / THREE
          d1 = - THREE / FIVE

      alpha2 = ONE / TWENTYTWO
          a2 = TWELVE / ELEVEN
          b2 = ZERO ! not used
          c2 = ZERO ! not used
      
    else  ! default 2nd CD
     alpha1 = ZERO
          a1 = -TWO
          b1 = THREE
          c1 = -ONE

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    end if
  
    d1fC2P(1,   1, IBC_DIRICHLET) = ZERO ! not used
    d1fC2P(1,   2, IBC_DIRICHLET) = ONE
    d1fC2P(1,   3, IBC_DIRICHLET) = alpha1
    d1rC2P(1,   1, IBC_DIRICHLET) = a1
    d1rC2P(1,   2, IBC_DIRICHLET) = b1
    d1rC2P(1,   3, IBC_DIRICHLET) = c1
    d1rC2P(1,   4, IBC_DIRICHLET) = d1

    d1fC2P(5,   1, IBC_DIRICHLET) =   d1fC2P(1, 3, IBC_DIRICHLET)
    d1fC2P(5,   2, IBC_DIRICHLET) =   d1fC2P(1, 2, IBC_DIRICHLET)
    d1fC2P(5,   3, IBC_DIRICHLET) =   d1fC2P(1, 1, IBC_DIRICHLET)
    d1rC2P(5,   1, IBC_DIRICHLET) = - d1rC2P(1, 1, IBC_DIRICHLET)
    d1rC2P(5,   2, IBC_DIRICHLET) = - d1rC2P(1, 2, IBC_DIRICHLET)
    d1rC2P(5,   3, IBC_DIRICHLET) = - d1rC2P(1, 3, IBC_DIRICHLET)
    d1rC2P(5,   4, IBC_DIRICHLET) = - d1rC2P(1, 4, IBC_DIRICHLET)

    d1fC2P(2,   1, IBC_DIRICHLET) = alpha2
    d1fC2P(2,   2, IBC_DIRICHLET) = ONE
    d1fC2P(2,   3, IBC_DIRICHLET) = alpha2
    d1rC2P(2,   1, IBC_DIRICHLET) = a2
    d1rC2P(2,   2, IBC_DIRICHLET) = b2 / THREE ! not used
    d1rC2P(2,   3, IBC_DIRICHLET) = c2 ! not used

    d1fC2P(4,   1, IBC_DIRICHLET) = d1fC2P(2, 1, IBC_DIRICHLET)
    d1fC2P(4,   2, IBC_DIRICHLET) = d1fC2P(2, 2, IBC_DIRICHLET)
    d1fC2P(4,   3, IBC_DIRICHLET) = d1fC2P(2, 3, IBC_DIRICHLET)
    d1rC2P(4,   1, IBC_DIRICHLET) = d1rC2P(2, 1, IBC_DIRICHLET)
    d1rC2P(4,   2, IBC_DIRICHLET) = d1rC2P(2, 2, IBC_DIRICHLET)
    d1rC2P(4,   3, IBC_DIRICHLET) = d1rC2P(2, 3, IBC_DIRICHLET)

    d1fC2P(3,   :, IBC_DIRICHLET) = d1fC2P(3, :, IBC_PERIODIC)
    d1rC2P(3,   :, IBC_DIRICHLET) = d1rC2P(3, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2C : no specified = Dirichlet B.C.
! [ 1     alpha1                            ][f'_1]=[a1 * f_{1'}/h  + b1 * f_{2'}/h + c1 * f_{3'}/h  ]
! [alpha2 1      alpha2                     ][f'_2] [a2 * (f_{3'} - f_{2'})/h  ]
! [       alpha  1      alpha               ][f'_i] [a *  (f_{i'+1} - f_{i'})/h + b/3 * (f_{i'+2} - f_{i'-1})/h]
! [                     alpha2 1      alpha2][f'_4] [a2 * (f_{n'} - f_{n'-1})/h]
! [                            alpha1 1     ][f'_5] [-a1 * f_{n'+1}/h  - b1 * f_{n'}/h - c1 * f_{n'-1}/h]
!-------------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
      alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then
      alpha1 = ZERO
          a1 = -ONE
          b1 = ONE
          c1 = ZERO

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = -ONE
          a1 = -ONE
          b1 = TWO
          c1 = -ONE

      alpha2 = ONE / TWENTYTWO
          a2 = TWELVE / ELEVEN
          b2 = ZERO ! not used
          c2 = ZERO ! not used
      
    else  ! default 2nd CD
      alpha1 = ZERO
          a1 = -ONE
          b1 = ONE
          c1 = ZERO

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used
    end if

    d1fP2C(1, 1, IBC_INTRPL) = ZERO ! not used
    d1fP2C(1, 2, IBC_INTRPL) = ONE
    d1fP2C(1, 3, IBC_INTRPL) = alpha1
    d1rP2C(1, 1, IBC_INTRPL) = a1
    d1rP2C(1, 2, IBC_INTRPL) = b1
    d1rP2C(1, 3, IBC_INTRPL) = c1

    d1fP2C(5, 1, IBC_INTRPL) =   d1fP2C(1, 3, IBC_INTRPL)
    d1fP2C(5, 2, IBC_INTRPL) =   d1fP2C(1, 2, IBC_INTRPL)
    d1fP2C(5, 3, IBC_INTRPL) =   d1fP2C(1, 1, IBC_INTRPL)
    d1rP2C(5, 1, IBC_INTRPL) = - d1rP2C(1, 1, IBC_INTRPL)
    d1rP2C(5, 2, IBC_INTRPL) = - d1rP2C(1, 2, IBC_INTRPL)
    d1rP2C(5, 3, IBC_INTRPL) = - d1rP2C(1, 3, IBC_INTRPL)

    d1fP2C(2, 1, IBC_INTRPL) = alpha2
    d1fP2C(2, 2, IBC_INTRPL) = ONE
    d1fP2C(2, 3, IBC_INTRPL) = alpha2
    d1rP2C(2, 1, IBC_INTRPL) = a2
    d1rP2C(2, 2, IBC_INTRPL) = b2 / THREE ! not used
    d1rP2C(2, 3, IBC_INTRPL) = c2 ! not used

    d1fP2C(4, 1, IBC_INTRPL) = d1fP2C(2, 1, IBC_INTRPL)
    d1fP2C(4, 2, IBC_INTRPL) = d1fP2C(2, 2, IBC_INTRPL)
    d1fP2C(4, 3, IBC_INTRPL) = d1fP2C(2, 3, IBC_INTRPL)
    d1rP2C(4, 1, IBC_INTRPL) = d1rP2C(2, 1, IBC_INTRPL)
    d1rP2C(4, 2, IBC_INTRPL) = d1rP2C(2, 2, IBC_INTRPL)
    d1rP2C(4, 3, IBC_INTRPL) = d1rP2C(2, 3, IBC_INTRPL)

    d1fP2C(3, :, IBC_INTRPL) = d1fP2C(3, :, IBC_PERIODIC)
    d1rP2C(3, :, IBC_INTRPL) = d1rP2C(3, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 1st-derivative : 
! P2C : no specified = Dirichlet B.C. = Neumann
!-------------------------------------------------------------------------------
    d1fP2C(:, :, IBC_DIRICHLET) = d1fP2C(:, :, IBC_INTRPL)
    d1rP2C(:, :, IBC_DIRICHLET) = d1rP2C(:, :, IBC_INTRPL)

    d1fP2C(:, :, IBC_NEUMANN  ) = d1fP2C(:, :, IBC_INTRPL)
    d1rP2C(:, :, IBC_NEUMANN  ) = d1rP2C(:, :, IBC_INTRPL)
!===============================================================================
!interpolation. P2C and C2P Periodic or Symmetric B.C.
! P2C : i_max = nc
! alpha * f_{i-1} + f_i + alpha * f_{i+1} =    a/2 * ( f_{i'}   + f_{i'+1} ) + &
!                                              b/2 * ( f_{i'+2} + f_{i'-1} )
! C2P : i'_max = np
! alpha * f_{i'-1} + f_i' + alpha * f_{i'+1} = a/2 * ( f_{i}   + f_{i-1} ) + &
!                                              b/2 * ( f_{i+1} + f_{i-2} )
!===============================================================================
    alpha = ZERO
        a = ONE
        b = ZERO
        c = ZERO
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
          a = ONE
          b = ZERO
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
          a = NINE / EIGHT
          b = -ONE / EIGHT
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / SIX
          a = FOUR / THREE
          b = ZERO
    else if (iaccu == IACCU_CP6) then
      alpha = THREE / TEN
          a = THREE / TWO
          b = ONE / TEN
    else  ! default 2nd CD
      alpha = ZERO
          a = ONE
          b = ZERO
    end if
!-------------------------------------------------------------------------------
! interpolation. C2P: periodic & symmetric 
! [ 1    alpha                   alpha][f_1']=[a/2 * (f_{1}   + f_{n})   + b/2 * (f_{2}   + f_{n-1})]
! [      alpha 1     alpha            ][f_2'] [a/2 * (f_{3}   + f_{n})   + b/2 * (f_{2}   + f_{1})  ]
! [            alpha 1     alpha      ][f_i'] [a/2 * (f_{i}   + f_{i-1}) + b/2 * (f_{i+1} + f_{i-2})]
! [                  alpha 1     alpha][f_4'] [a/2 * (f_{n}   + f_{n-3}) + b/2 * (f_{n-1} + f_{n-2})]
! [alpha                   alpha 1    ][f_5'] [a/2 * (f_{1}   + f_{n-2}) + b/2 * (f_{n}   + f_{n-1})]
!-----------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!interpolation. C2P for periodic b.c.
!-------------------------------------------------------------------------------
    m1fC2P(1:5, 1, IBC_PERIODIC) = alpha
    m1fC2P(1:5, 2, IBC_PERIODIC) = ONE
    m1fC2P(1:5, 3, IBC_PERIODIC) = alpha
    m1rC2P(1:5, 1, IBC_PERIODIC) = a / TWO
    m1rC2P(1:5, 2, IBC_PERIODIC) = b / TWO
    m1rC2P(1:5, 3, IBC_PERIODIC) = c ! not used
!-------------------------------------------------------------------------------
!interpolation. P2C for periodic b.c.
!-------------------------------------------------------------------------------
    m1fP2C(:, :, IBC_PERIODIC) = m1fC2P(:, :, IBC_PERIODIC)
    m1rP2C(:, :, IBC_PERIODIC) = m1rC2P(:, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
!interpolation. C2P. symmetric, orthogonal, eg. u in y direction.
!-------------------------------------------------------------------------------
    m1fC2P(1,   1, IBC_SYMMETRIC) = ZERO ! not used
    m1fC2P(1,   2, IBC_SYMMETRIC) = ONE
    m1fC2P(1,   3, IBC_SYMMETRIC) = alpha + alpha
    m1fC2P(5,   1, IBC_SYMMETRIC) = m1fC2P(1,   3, IBC_SYMMETRIC)
    m1fC2P(5,   2, IBC_SYMMETRIC) = m1fC2P(1,   2, IBC_SYMMETRIC)
    m1fC2P(5,   3, IBC_SYMMETRIC) = m1fC2P(1,   1, IBC_SYMMETRIC)
    m1fC2P(2:4, :, IBC_SYMMETRIC) = m1fC2P(2:4, :, IBC_PERIODIC)

    m1rC2P(:,   :, IBC_SYMMETRIC) = m1rC2P(:, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
!interpolation. C2P. asymmetric, orthogonal, eg. v in y direction.
!-------------------------------------------------------------------------------
    m1fC2P(1,   1, IBC_ASYMMETRIC) = ZERO ! not used
    m1fC2P(1,   2, IBC_ASYMMETRIC) = ONE
    m1fC2P(1,   3, IBC_ASYMMETRIC) = ZERO
    m1fC2P(5,   1, IBC_ASYMMETRIC) = m1fC2P(1,   3, IBC_ASYMMETRIC)
    m1fC2P(5,   2, IBC_ASYMMETRIC) = m1fC2P(1,   2, IBC_ASYMMETRIC)
    m1fC2P(5,   3, IBC_ASYMMETRIC) = m1fC2P(1,   1, IBC_ASYMMETRIC)
    m1fC2P(2:4, :, IBC_ASYMMETRIC) = m1fC2P(2:4, :, IBC_PERIODIC)

    m1rC2P(:,   :, IBC_ASYMMETRIC) = m1rC2P(:, :, IBC_PERIODIC)
    m1rC2P(1,   :, IBC_ASYMMETRIC) = ZERO ! double safe, not necessary
    m1rC2P(5,   :, IBC_ASYMMETRIC) = ZERO
!-------------------------------------------------------------------------------
!interpolation. P2C. symmetric, orthogonal, eg. u in y direction.
!-------------------------------------------------------------------------------
    m1fP2C(1,   1, IBC_SYMMETRIC) = ZERO ! not used
    m1fP2C(1,   2, IBC_SYMMETRIC) = ONE + alpha
    m1fP2C(1,   3, IBC_SYMMETRIC) = alpha
    m1fP2C(5,   1, IBC_SYMMETRIC) = m1fP2C(1,   3, IBC_SYMMETRIC)
    m1fP2C(5,   2, IBC_SYMMETRIC) = m1fP2C(1,   2, IBC_SYMMETRIC)
    m1fP2C(5,   3, IBC_SYMMETRIC) = m1fP2C(1,   1, IBC_SYMMETRIC)
    m1fP2C(2:4, :, IBC_SYMMETRIC) = m1fP2C(2:4, :, IBC_PERIODIC)

    m1rP2C(:,   :, IBC_SYMMETRIC) = m1rP2C(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
!interpolation. P2C. asymmetric, orthogonal, eg. v in y direction.
!-------------------------------------------------------------------------------
    m1fP2C(1,   1, IBC_ASYMMETRIC) = ZERO ! not used
    m1fP2C(1,   2, IBC_ASYMMETRIC) = ONE - alpha
    m1fP2C(1,   3, IBC_ASYMMETRIC) = alpha
    m1fP2C(5,   1, IBC_ASYMMETRIC) = m1fP2C(1,   3, IBC_ASYMMETRIC)
    m1fP2C(5,   2, IBC_ASYMMETRIC) = m1fP2C(1,   2, IBC_ASYMMETRIC)
    m1fP2C(5,   3, IBC_ASYMMETRIC) = m1fP2C(1,   1, IBC_ASYMMETRIC)
    m1fP2C(2:4, :, IBC_ASYMMETRIC) = m1fP2C(2:4, :, IBC_PERIODIC)

    m1rP2C(:,   :, IBC_ASYMMETRIC) = m1rP2C(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! interpolation. P2C: Dirichlet = no specified = Neumann
! [ 1    alpha1                          ][f_1']=[a1 * f_{1'} + b1 * f_{2'} + c1 * f_{3'}  ]
! [      alpha2 1     alpha2             ][f_2'] [a2/2 * (f_{2'}   + f_{3'})]
! [             alpha 1      alpha       ][f_i'] [a/2  * (f_{i'}   + f_{i'+1}) + b/2 * (f_{i'+2} + f_{i'-1})]
! [                   alpha2 1     alpha2][f_4'] [a2/2 * (f_{n'-1}   + f_{n'})]
! [                          alpha1 1    ][f_5'] [a1   * f_{n'+1} + b1 * f_{n'} + c1 * f_{n'-1}]
!-----------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
    alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
          a1 = THREE / EIGHT
          b1 = THREE / FOUR
          c1 = -ONE / EIGHT

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = ONE
          a1 = ONE / FOUR
          b1 = THREE / TWO
          c1 = ONE / FOUR

      alpha2 = ONE / SIX
          a2 = FOUR / THREE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else  ! default 2nd CD
      alpha1 = ZERO
          a1 = THREE / EIGHT
          b1 = THREE / FOUR
          c1 = -ONE / EIGHT

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used
    end if
    !P2C
    m1fP2C(1, 1, IBC_INTRPL) = ZERO ! not used
    m1fP2C(1, 2, IBC_INTRPL) = ONE
    m1fP2C(1, 3, IBC_INTRPL) = alpha1
    m1rP2C(1, 1, IBC_INTRPL) = a1
    m1rP2C(1, 2, IBC_INTRPL) = b1
    m1rP2C(1, 3, IBC_INTRPL) = c1

    m1fP2C(5, 1, IBC_INTRPL) = m1fP2C(1, 3, IBC_INTRPL)
    m1fP2C(5, 2, IBC_INTRPL) = m1fP2C(1, 2, IBC_INTRPL)
    m1fP2C(5, 3, IBC_INTRPL) = m1fP2C(1, 1, IBC_INTRPL)
    m1rP2C(5, 1, IBC_INTRPL) = m1rP2C(1, 1, IBC_INTRPL)
    m1rP2C(5, 2, IBC_INTRPL) = m1rP2C(1, 2, IBC_INTRPL)
    m1rP2C(5, 3, IBC_INTRPL) = m1rP2C(1, 3, IBC_INTRPL)

    m1fP2C(2:4, :, IBC_INTRPL) = m1fP2C(2:4, :, IBC_PERIODIC)
    m1rP2C(2:4, :, IBC_INTRPL) = m1rP2C(2:4, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! interpolation. P2C: Dirichlet = no specified = Neumann
!-------------------------------------------------------------------------------
    m1fP2C(:, :, IBC_DIRICHLET) = m1fP2C(:, :, IBC_INTRPL)
    m1rP2C(:, :, IBC_DIRICHLET) = m1rP2C(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! interpolation. P2C: Dirichlet = no specified = Neumann
!-------------------------------------------------------------------------------
    m1fP2C(:, :, IBC_NEUMANN) = m1fP2C(:, :, IBC_INTRPL)
    m1rP2C(:, :, IBC_NEUMANN) = m1rP2C(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! interpolation. C2P: No specified = neumann
! [ 1    alpha1                          ][f_1']=[a1 * f_{1} + b1 * f_{2} + c1 * f_{3}  ]
! [      alpha2 1     alpha2             ][f_2'] [a2/2 * (f_{2}   + f_{1})]
! [             alpha 1      alpha       ][f_i'] [a/2 * (f_{i}   + f_{i-1}) + b/2 * (f_{i+1} + f_{i-2})]
! [                   alpha2 1     alpha2][f_4'] [a2/2 * (f_{n-1}   + f_{n-2})]
! [                          alpha1 1    ][f_5'] [a1 * f_{n-1} + b1 * f_{n-2} + c1 * f_{n-3}]
!-----------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
      alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
          a1 = FIFTEEN / EIGHT
          b1 = - FIVE / FOUR
          c1 = THREE / EIGHT

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = FIVE
          a1 = FIFTEEN / FOUR
          b1 = FIVE / TWO
          c1 = -ONE / FOUR

      alpha2 = ONE / SIX 
          a2 = FOUR / THREE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else  ! default 2nd CD

      alpha1 = ZERO
          a1 = FIFTEEN / EIGHT
          b1 = - FIVE / FOUR
          c1 = THREE / EIGHT

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used
     
    end if

    m1fC2P(1, 1, IBC_INTRPL) = ZERO ! not used
    m1fC2P(1, 2, IBC_INTRPL) = ONE
    m1fC2P(1, 3, IBC_INTRPL) = alpha1
    m1rC2P(1, 1, IBC_INTRPL) = a1
    m1rC2P(1, 2, IBC_INTRPL) = b1
    m1rC2P(1, 3, IBC_INTRPL) = c1

    m1fC2P(5, 1, IBC_INTRPL) = m1fC2P(1, 3, IBC_INTRPL)
    m1fC2P(5, 2, IBC_INTRPL) = m1fC2P(1, 2, IBC_INTRPL)
    m1fC2P(5, 3, IBC_INTRPL) = m1fC2P(1, 1, IBC_INTRPL)
    m1rC2P(5, 1, IBC_INTRPL) = m1rC2P(1, 1, IBC_INTRPL)
    m1rC2P(5, 2, IBC_INTRPL) = m1rC2P(1, 2, IBC_INTRPL)
    m1rC2P(5, 3, IBC_INTRPL) = m1rC2P(1, 3, IBC_INTRPL)

    m1fC2P(2, 1, IBC_INTRPL) = alpha2
    m1fC2P(2, 2, IBC_INTRPL) = ONE
    m1fC2P(2, 3, IBC_INTRPL) = alpha2
    m1rC2P(2, 1, IBC_INTRPL) = a2 / TWO
    m1rC2P(2, 2, IBC_INTRPL) = ZERO ! not used
    m1rC2P(2, 3, IBC_INTRPL) = ZERO ! not used

    m1fC2P(4, 1, IBC_INTRPL) = m1fC2P(2, 1, IBC_INTRPL)
    m1fC2P(4, 2, IBC_INTRPL) = m1fC2P(2, 2, IBC_INTRPL)
    m1fC2P(4, 3, IBC_INTRPL) = m1fC2P(2, 3, IBC_INTRPL)
    m1rC2P(4, 1, IBC_INTRPL) = m1rC2P(2, 1, IBC_INTRPL)
    m1rC2P(4, 2, IBC_INTRPL) = m1rC2P(2, 2, IBC_INTRPL)
    m1rC2P(4, 3, IBC_INTRPL) = m1rC2P(2, 3, IBC_INTRPL)

    m1fC2P(3, :, IBC_INTRPL) = m1fC2P(3, :, IBC_PERIODIC)
    m1rC2P(3, :, IBC_INTRPL) = m1rC2P(3, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! interpolation. C2P: No specified = neumann
!-------------------------------------------------------------------------------
    m1fC2P(:, :, IBC_NEUMANN) = m1fC2P(:, :, IBC_INTRPL)
    m1rC2P(:, :, IBC_NEUMANN) = m1rC2P(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! interpolation. C2P: Dirichlet
! [ 1    0                              ][f_1']=known
! [      alpha2 1     alpha2             ][f_2'] [a2/2 * (f_{2}   + f_{1})]
! [             alpha 1      alpha       ][f_i'] [a/2 * (f_{i}   + f_{i-1}) + b/2 * (f_{i+1} + f_{i-2})]
! [                   alpha2 1     alpha2][f_4'] [a2/2 * (f_{n-1}   + f_{n-2})]
! [                          0     1    ][f_5'] known
!-----------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
      alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
          a1 = ZERO ! not used
          b1 = ZERO ! not used
          c1 = ZERO ! not used

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = ZERO
          a1 = ZERO ! not used
          b1 = ZERO ! not used
          c1 = ZERO ! not used

      alpha2 = ONE / SIX 
          a2 = FOUR / THREE
          b2 = ZERO ! not used
          c2 = ZERO ! not used

    else  ! default 2nd CD

      alpha1 = ZERO
          a1 = ZERO ! not used
          b1 = ZERO ! not used
          c1 = ZERO ! not used

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO ! not used
          c2 = ZERO ! not used
     
    end if

    m1fC2P(1, 1, IBC_DIRICHLET) = ZERO ! not used
    m1fC2P(1, 2, IBC_DIRICHLET) = ONE
    m1fC2P(1, 3, IBC_DIRICHLET) = alpha1
    m1rC2P(1, 1, IBC_DIRICHLET) = a1
    m1rC2P(1, 2, IBC_DIRICHLET) = b1
    m1rC2P(1, 3, IBC_DIRICHLET) = c1

    m1fC2P(5, 1, IBC_DIRICHLET) = m1fC2P(1, 3, IBC_DIRICHLET)
    m1fC2P(5, 2, IBC_DIRICHLET) = m1fC2P(1, 2, IBC_DIRICHLET)
    m1fC2P(5, 3, IBC_DIRICHLET) = m1fC2P(1, 1, IBC_DIRICHLET)
    m1rC2P(5, 1, IBC_DIRICHLET) = m1rC2P(1, 1, IBC_DIRICHLET)
    m1rC2P(5, 2, IBC_DIRICHLET) = m1rC2P(1, 2, IBC_DIRICHLET)
    m1rC2P(5, 3, IBC_DIRICHLET) = m1rC2P(1, 3, IBC_DIRICHLET)

    m1fC2P(2, 1, IBC_DIRICHLET) = alpha2
    m1fC2P(2, 2, IBC_DIRICHLET) = ONE
    m1fC2P(2, 3, IBC_DIRICHLET) = alpha2
    m1rC2P(2, 1, IBC_DIRICHLET) = a2 / TWO
    m1rC2P(2, 2, IBC_DIRICHLET) = ZERO ! not used
    m1rC2P(2, 3, IBC_DIRICHLET) = ZERO ! not used

    m1fC2P(4, 1, IBC_DIRICHLET) = m1fC2P(2, 1, IBC_DIRICHLET)
    m1fC2P(4, 2, IBC_DIRICHLET) = m1fC2P(2, 2, IBC_DIRICHLET)
    m1fC2P(4, 3, IBC_DIRICHLET) = m1fC2P(2, 3, IBC_DIRICHLET)
    m1rC2P(4, 1, IBC_DIRICHLET) = m1rC2P(2, 1, IBC_DIRICHLET)
    m1rC2P(4, 2, IBC_DIRICHLET) = m1rC2P(2, 2, IBC_DIRICHLET)
    m1rC2P(4, 3, IBC_DIRICHLET) = m1rC2P(2, 3, IBC_DIRICHLET)

    m1fC2P(3, :, IBC_DIRICHLET) = m1fC2P(3, :, IBC_PERIODIC)
    m1rC2P(3, :, IBC_DIRICHLET) = m1rC2P(3, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 2nd diriviative P2P and C2C, periodic, symmetric, asymmetric
!-------------------------------------------------------------------------------
    alpha = ZERO
        a = ONE
        b = ZERO
        c = ZERO ! not used
        d = ZERO ! not used
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
          a = ONE
          b = ZERO
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
          a = FOUR / THREE
          b = -ONE / THREE
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / TEN
          a = SIX / FIVE
          b = ZERO
    else if (iaccu == IACCU_CP6) then
      alpha = TWO / ELEVEN
          a = TWELVE / ELEVEN
          b = THREE / ELEVEN
    else  ! default 2nd CD
      alpha = ZERO
          a = ONE
          b = ZERO
    end if
!-------------------------------------------------------------------------------
! 2nd diriviative C2C, periodic
!-------------------------------------------------------------------------------
    d2fC2C(1:5, 1, IBC_PERIODIC) = alpha
    d2fC2C(1:5, 2, IBC_PERIODIC) = ONE
    d2fC2C(1:5, 3, IBC_PERIODIC) = alpha

    d2rC2C(1:5, 1, IBC_PERIODIC) = a / ONE
    d2rC2C(1:5, 2, IBC_PERIODIC) = b / FOUR
    d2rC2C(1:5, 3, IBC_PERIODIC) = c ! not used
    d2rC2C(1:5, 4, IBC_PERIODIC) = d ! not used
!-------------------------------------------------------------------------------
! 2nd diriviative P2P , periodic
!-------------------------------------------------------------------------------
    d2fP2P(:, :, IBC_PERIODIC) = d2fC2C(:, :, IBC_PERIODIC)
    d2rP2P(:, :, IBC_PERIODIC) = d2rC2C(:, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 2nd diriviative C2C, symmetric
!-------------------------------------------------------------------------------
    d2fC2C(1,   1, IBC_SYMMETRIC) = ZERO ! not used
    d2fC2C(1,   2, IBC_SYMMETRIC) = ONE + alpha
    d2fC2C(1,   3, IBC_SYMMETRIC) = alpha
    d2fC2C(5,   1, IBC_SYMMETRIC) = d2fC2C(1,   3, IBC_SYMMETRIC)
    d2fC2C(5,   2, IBC_SYMMETRIC) = d2fC2C(1,   2, IBC_SYMMETRIC)
    d2fC2C(5,   3, IBC_SYMMETRIC) = d2fC2C(1,   1, IBC_SYMMETRIC)
    d2fC2C(2:4, :, IBC_SYMMETRIC) = d2fC2C(2:4, :, IBC_PERIODIC)

    d2rC2C(:,   :, IBC_SYMMETRIC) = d2rC2C(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 2nd diriviative C2C, asymmetric
!-------------------------------------------------------------------------------
    d2fC2C(1,   1, IBC_ASYMMETRIC) = ZERO ! not used
    d2fC2C(1,   2, IBC_ASYMMETRIC) = ONE - alpha
    d2fC2C(1,   3, IBC_ASYMMETRIC) = alpha
    d2fC2C(5,   1, IBC_ASYMMETRIC) = d2fC2C(1,   3, IBC_ASYMMETRIC)
    d2fC2C(5,   2, IBC_ASYMMETRIC) = d2fC2C(1,   2, IBC_ASYMMETRIC)
    d2fC2C(5,   3, IBC_ASYMMETRIC) = d2fC2C(1,   1, IBC_ASYMMETRIC)
    d2fC2C(2:4, :, IBC_ASYMMETRIC) = d2fC2C(2:4, :, IBC_PERIODIC)

    d2rC2C(:,   :, IBC_ASYMMETRIC) = d2rC2C(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 2nd diriviative P2P, symmetric
!-------------------------------------------------------------------------------
    d2fP2P(1,   1, IBC_SYMMETRIC) = ZERO ! not used
    d2fP2P(1,   2, IBC_SYMMETRIC) = ONE
    d2fP2P(1,   3, IBC_SYMMETRIC) = alpha + alpha
    d2fP2P(5,   1, IBC_SYMMETRIC) = d2fP2P(1,   3, IBC_SYMMETRIC)
    d2fP2P(5,   2, IBC_SYMMETRIC) = d2fP2P(1,   2, IBC_SYMMETRIC)
    d2fP2P(5,   3, IBC_SYMMETRIC) = d2fP2P(1,   1, IBC_SYMMETRIC)
    d2fP2P(2:4, :, IBC_SYMMETRIC) = d2fP2P(2:4, :, IBC_PERIODIC)
    d2rP2P(:,   :, IBC_SYMMETRIC) = d2rP2P(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 2nd diriviative P2P, asymmetric
!-------------------------------------------------------------------------------
    d2fP2P(1,   1, IBC_ASYMMETRIC) = ZERO ! not used
    d2fP2P(1,   2, IBC_ASYMMETRIC) = ONE
    d2fP2P(1,   3, IBC_ASYMMETRIC) = ZERO
    d2fP2P(5,   1, IBC_ASYMMETRIC) = d2fP2P(1,   3, IBC_ASYMMETRIC)
    d2fP2P(5,   2, IBC_ASYMMETRIC) = d2fP2P(1,   2, IBC_ASYMMETRIC)
    d2fP2P(5,   3, IBC_ASYMMETRIC) = d2fP2P(1,   1, IBC_ASYMMETRIC)
    d2fP2P(2:4, :, IBC_ASYMMETRIC) = d2fP2P(2:4, :, IBC_PERIODIC)
    d2rP2P(:,   :, IBC_ASYMMETRIC) = d2rP2P(:,   :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 2nd diriviative P2P, no specified = dirichlet = neumann 
!-------------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
        d1 = ZERO

    alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
        d2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then
      alpha1 = ZERO
          a1 = TWO
          b1 = -FIVE
          c1 = FOUR
          d1 = -ONE

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
          c2 = ZERO ! not used
          d2 = ZERO ! not used
    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6 ) then ! degrade to 3rd CP
      alpha1 = ELEVEN
          a1 = THIRTEEN
          b1 = -TWENTYSEVEN
          c1 = FIFTEEN
          d1 = -ONE

      alpha2 = ONE / TEN
          a2 = SIX / FIVE
          b2 = ZERO
          c2 = ZERO ! not used
          d2 = ZERO ! not used
    else  ! default 2nd CD
      alpha1 = ZERO
          a1 = TWO
          b1 = -FIVE
          c1 = FOUR
          d1 = -ONE

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
          c2 = ZERO ! not used
          d2 = ZERO ! not used
    end if

    d2fP2P(1, 1, IBC_INTRPL) = ZERO ! not used
    d2fP2P(1, 2, IBC_INTRPL) = ONE
    d2fP2P(1, 3, IBC_INTRPL) = alpha1
    d2rP2P(1, 1, IBC_INTRPL) = a1
    d2rP2P(1, 2, IBC_INTRPL) = b1
    d2rP2P(1, 3, IBC_INTRPL) = c1
    d2rP2P(1, 4, IBC_INTRPL) = d1

    d2fP2P(5, 1, IBC_INTRPL) = d2fP2P(1, 3, IBC_INTRPL)
    d2fP2P(5, 2, IBC_INTRPL) = d2fP2P(1, 2, IBC_INTRPL)
    d2fP2P(5, 3, IBC_INTRPL) = d2fP2P(1, 1, IBC_INTRPL)
    d2rP2P(5, 1, IBC_INTRPL) = d2rP2P(1, 1, IBC_INTRPL)
    d2rP2P(5, 2, IBC_INTRPL) = d2rP2P(1, 2, IBC_INTRPL)
    d2rP2P(5, 3, IBC_INTRPL) = d2rP2P(1, 3, IBC_INTRPL)
    d2rP2P(5, 4, IBC_INTRPL) = d2rP2P(1, 4, IBC_INTRPL)

    d2fP2P(2, 1, IBC_INTRPL) = alpha2
    d2fP2P(2, 2, IBC_INTRPL) = ONE
    d2fP2P(2, 3, IBC_INTRPL) = alpha2
    d2rP2P(2, 1, IBC_INTRPL) = a2
    d2rP2P(2, 2, IBC_INTRPL) = ZERO ! not used
    d2rP2P(2, 3, IBC_INTRPL) = ZERO ! not used
    d2rP2P(2, 4, IBC_INTRPL) = ZERO ! not used

    d2fP2P(4, 1, IBC_INTRPL) = d2fP2P(2, 1, IBC_INTRPL)
    d2fP2P(4, 2, IBC_INTRPL) = d2fP2P(2, 2, IBC_INTRPL)
    d2fP2P(4, 3, IBC_INTRPL) = d2fP2P(2, 3, IBC_INTRPL)
    d2rP2P(4, 1, IBC_INTRPL) = d2rP2P(2, 1, IBC_INTRPL)
    d2rP2P(4, 2, IBC_INTRPL) = d2rP2P(2, 2, IBC_INTRPL)
    d2rP2P(4, 3, IBC_INTRPL) = d2rP2P(2, 3, IBC_INTRPL)
    d2rP2P(4, 4, IBC_INTRPL) = d2rP2P(2, 4, IBC_INTRPL)

    d2fP2P(3, :, IBC_INTRPL) = d2fP2P(3, :, IBC_PERIODIC)
    d2rP2P(3, :, IBC_INTRPL) = d2rP2P(3, :, IBC_PERIODIC)
!-------------------------------------------------------------------------------
! 2nd diriviative P2P, no specified = dirichlet 
!-------------------------------------------------------------------------------
    d2fP2P(:, :, IBC_DIRICHLET) = d2fP2P(:, :, IBC_INTRPL)
    d2rP2P(:, :, IBC_DIRICHLET) = d2rP2P(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! 2nd diriviative P2P, no specified = neumann 
!-------------------------------------------------------------------------------
    d2fP2P(:, :, IBC_NEUMANN) = d2fP2P(:, :, IBC_INTRPL)
    d2rP2P(:, :, IBC_NEUMANN) = d2rP2P(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! 2nd diriviative C2C, no specified =  neumann = P2C unspecified
!-------------------------------------------------------------------------------
    d2fC2C(:, :, IBC_INTRPL) = d2fP2P(:, :, IBC_INTRPL)
    d2rC2C(:, :, IBC_INTRPL) = d2rP2P(:, :, IBC_INTRPL)

    d2fC2C(:, :, IBC_NEUMANN) = d2fC2C(:, :, IBC_INTRPL)
    d2rC2C(:, :, IBC_NEUMANN) = d2rC2C(:, :, IBC_INTRPL)
!-------------------------------------------------------------------------------
! 2nd diriviative C2C, Dirchilet
!-------------------------------------------------------------------------------
    alpha1 = ZERO
        a1 = ZERO
        b1 = ZERO
        c1 = ZERO
        d1 = ZERO

    alpha2 = ZERO
        a2 = ZERO
        b2 = ZERO
        c2 = ZERO
        d2 = ZERO
    if (iaccu == IACCU_CD2 .or. iaccu == IACCU_CD4) then
      alpha1 = ZERO
          a1 = SIXTEEN / SEVEN
          b1 = -TWENTYFIVE / SEVEN
          c1 = TEN / SEVEN
          d1 = -ONE / SEVEN

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
          c2 = ZERO ! not used
          d2 = ZERO ! not used
    else if (iaccu == IACCU_CP4 .or. iaccu == IACCU_CP6 ) then ! degrade to 3rd CP
      alpha1 = FIVE / FOURTEEN
          a1 = SIXTEEN / SEVEN
          b1 = -FOURTYFIVE / FOURTEEN
          c1 = FIVE / SEVEN
          d1 = THREE / FOURTEEN

      alpha2 = ONE / TEN
          a2 = SIX / FIVE
          b2 = ZERO
          c2 = ZERO ! not used
          d2 = ZERO ! not used
    else  ! default 2nd CD
      alpha1 = ZERO
          a1 = SIXTEEN / SEVEN
          b1 = -TWENTYFIVE / SEVEN
          c1 = TEN / SEVEN
          d1 = -ONE / SEVEN

      alpha2 = ZERO
          a2 = ONE
          b2 = ZERO
          c2 = ZERO ! not used
          d2 = ZERO ! not used
    end if

    d2fC2C(1, 1, IBC_DIRICHLET) = ZERO ! not used
    d2fC2C(1, 2, IBC_DIRICHLET) = ONE
    d2fC2C(1, 3, IBC_DIRICHLET) = alpha1
    d2rC2C(1, 1, IBC_DIRICHLET) = a1
    d2rC2C(1, 2, IBC_DIRICHLET) = b1
    d2rC2C(1, 3, IBC_DIRICHLET) = c1
    d2rC2C(1, 4, IBC_DIRICHLET) = d1

    d2fC2C(5, 1, IBC_DIRICHLET) = d2fC2C(1, 3, IBC_DIRICHLET)
    d2fC2C(5, 2, IBC_DIRICHLET) = d2fC2C(1, 2, IBC_DIRICHLET)
    d2fC2C(5, 3, IBC_DIRICHLET) = d2fC2C(1, 1, IBC_DIRICHLET)
    d2rC2C(5, 1, IBC_DIRICHLET) = d2rC2C(1, 1, IBC_DIRICHLET)
    d2rC2C(5, 2, IBC_DIRICHLET) = d2rC2C(1, 2, IBC_DIRICHLET)
    d2rC2C(5, 3, IBC_DIRICHLET) = d2rC2C(1, 3, IBC_DIRICHLET)
    d2rC2C(5, 4, IBC_DIRICHLET) = d2rC2C(1, 4, IBC_DIRICHLET)

    d2fC2C(2, 1, IBC_DIRICHLET) = alpha2
    d2fC2C(2, 2, IBC_DIRICHLET) = ONE
    d2fC2C(2, 3, IBC_DIRICHLET) = alpha2
    d2rC2C(2, 1, IBC_DIRICHLET) = a2
    d2rC2C(2, 2, IBC_DIRICHLET) = ZERO ! not used
    d2rC2C(2, 3, IBC_DIRICHLET) = ZERO ! not used
    d2rC2C(2, 4, IBC_DIRICHLET) = ZERO ! not used

    d2fC2C(4, 1, IBC_DIRICHLET) = d2fC2C(2, 1, IBC_DIRICHLET)
    d2fC2C(4, 2, IBC_DIRICHLET) = d2fC2C(2, 2, IBC_DIRICHLET)
    d2fC2C(4, 3, IBC_DIRICHLET) = d2fC2C(2, 3, IBC_DIRICHLET)
    d2rC2C(4, 1, IBC_DIRICHLET) = d2rC2C(2, 1, IBC_DIRICHLET)
    d2rC2C(4, 2, IBC_DIRICHLET) = d2rC2C(2, 2, IBC_DIRICHLET)
    d2rC2C(4, 3, IBC_DIRICHLET) = d2rC2C(2, 3, IBC_DIRICHLET)
    d2rC2C(4, 4, IBC_DIRICHLET) = d2rC2C(2, 4, IBC_DIRICHLET)
    
    d2fC2C(3, :, IBC_DIRICHLET) = d2fC2C(3, :, IBC_PERIODIC)
    d2rC2C(3, :, IBC_DIRICHLET) = d2rC2C(3, :, IBC_PERIODIC)


    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine Prepare_compact_coefficients
!===============================================================================
!> \brief Assigning the sparse matrix in the LHS of the compact scheme, and
!> calculating the geometry-only dependent variables for the TDMA scheme.
!>
!> This subroutine is called once locally.
!>
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           !
!-------------------------------------------------------------------------------
!> \param[in]     n             the number of unknown array
!> \param[in]     bc            the boundary condition at two ends of the unknown
!> \param[in]     coeff         the basic TDMA coefficients defined above.
!> \param[out]    a             the coefficients for TDMA
!> \param[out]    b             a_i * x_(i-1) + b_i * x_(i) + c_i * x_(i+1)
!> \param[out]    c             = RHS
!> \param[out]    d             An assisting coeffients for the TDMA scheme.
!-------------------------------------------------------------------------------
  subroutine Buildup_TDMA_LHS_array(n, is_periodic, coeff, a, b, c, d)
    use tridiagonal_matrix_algorithm
    use parameters_constant_mod
    implicit none

    integer, intent(in) :: n
    logical,  intent(in)   :: is_periodic
    real(WP), intent(in)   :: coeff(5, 3, 6)
    real(WP), intent(out)  :: a(n, 6, 6), &
                              b(n, 6, 6), &
                              c(n, 6, 6), &
                              d(n, 6, 6)

    integer :: i, j

    a(:, :, :) =  ZERO
    b(:, :, :) =  ZERO
    c(:, :, :) =  ZERO
    d(:, :, :) =  ZERO

    do j = IBC_PERIODIC, IBC_INTRPL
      do i = IBC_PERIODIC, IBC_INTRPL

        if (j == IBC_PERIODIC .and. i /= IBC_PERIODIC) cycle
        if (j /= IBC_PERIODIC .and. i == IBC_PERIODIC) cycle

        a(1,         i, j) = coeff( 1, 1, i )
        a(2,         i, j) = coeff( 2, 1, i )
        a(3 : n - 2, i, j) = coeff( 3, 1, IBC_PERIODIC )
        a(n - 1,     i, j) = coeff( 4, 1, j )
        a(n,         i, j) = coeff( 5, 1, j )

        b(1,         i, j) = coeff( 1, 2, i )
        b(2,         i, j) = coeff( 2, 2, i )
        b(3 : n - 2, i, j) = coeff( 3, 2, IBC_PERIODIC )
        b(n - 1,     i, j) = coeff( 4, 2, j )
        b(n,         i, j) = coeff( 5, 2, j )

        c(1,         i, j) = coeff( 1, 3, i )
        c(2,         i, j) = coeff( 2, 3, i )
        c(3 : n - 2, i, j) = coeff( 3, 3, IBC_PERIODIC )
        c(n - 1,     i, j) = coeff( 4, 3, j )
        c(n,         i, j) = coeff( 5, 3, j )

        if (is_periodic) then
          call Preprocess_TDMA_coeffs(a(1:n-1, i, j), b(1:n-1, i, j), c(1:n-1, i, j), d(1:n-1, i, j), n-1)
        else 
          call Preprocess_TDMA_coeffs(a(:, i, j), b(:, i, j), c(:, i, j), d(:, i, j), n)
        end if 
      end do
    end do

    return
  end subroutine Buildup_TDMA_LHS_array
!===============================================================================
!> \brief Preparing coefficients for TDMA calculation.
!-------------------------------------------------------------------------------
!> Scope:  mpi    called-freq    xdomain
!>         all    once           all
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           !
!-------------------------------------------------------------------------------
!===============================================================================
  subroutine Prepare_LHS_coeffs_for_operations
    use vars_df_mod, only : domain
    use mpi_mod
    use parameters_constant_mod
    implicit none
    integer :: i, nsz

!===============================================================================
!   building up the basic lhs coeffients for compact schemes, based on the given
!   accuracy
!===============================================================================
    call Prepare_compact_coefficients (domain(1)%iAccuracy)
!-------------------------------------------------------------------------------
!   building up the full size lhs coeffients for compact schemes
!-------------------------------------------------------------------------------
!===============================================================================
! y-direction, with nc unknows
!===============================================================================
    i = 2
    nsz = domain(1)%nc(i)
!-------------------------------------------------------------------------------
!   1st derivative in y direction with nc unknows
!-------------------------------------------------------------------------------
    allocate (ad1y_C2C ( nsz, 6, 6 ) ); ad1y_C2C(:, :, :) = ZERO
    allocate (bd1y_C2C ( nsz, 6, 6 ) ); bd1y_C2C(:, :, :) = ZERO
    allocate (cd1y_C2C ( nsz, 6, 6 ) ); cd1y_C2C(:, :, :) = ZERO
    allocate (dd1y_C2C ( nsz, 6, 6 ) ); dd1y_C2C(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fC2C, &
          ad1y_C2C, bd1y_C2C, cd1y_C2C, dd1y_C2C)

    allocate (ad1y_P2C ( nsz, 6, 6 ) ); ad1y_P2C(:, :, :) = ZERO
    allocate (bd1y_P2C ( nsz, 6, 6 ) ); bd1y_P2C(:, :, :) = ZERO
    allocate (cd1y_P2C ( nsz, 6, 6 ) ); cd1y_P2C(:, :, :) = ZERO
    allocate (dd1y_P2C ( nsz, 6, 6 ) ); dd1y_P2C(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fP2C, &
          ad1y_P2C, bd1y_P2C, cd1y_P2C, dd1y_P2C)
!-------------------------------------------------------------------------------
!   mid-point interpolation in y direction with nc unknows
!-------------------------------------------------------------------------------
    allocate (am1y_P2C ( nsz, 6, 6 ) ); am1y_P2C(:, :, :) = ZERO
    allocate (bm1y_P2C ( nsz, 6, 6 ) ); bm1y_P2C(:, :, :) = ZERO
    allocate (cm1y_P2C ( nsz, 6, 6 ) ); cm1y_P2C(:, :, :) = ZERO
    allocate (dm1y_P2C ( nsz, 6, 6 ) ); dm1y_P2C(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), m1fP2C, &
          am1y_P2C, bm1y_P2C, cm1y_P2C, dm1y_P2C)
!-------------------------------------------------------------------------------
!   2nd order deriviative in y direction with nc unknows
!-------------------------------------------------------------------------------
    allocate (ad2y_C2C ( nsz, 6, 6 ) ); ad2y_C2C(:, :, :) = ZERO
    allocate (bd2y_C2C ( nsz, 6, 6 ) ); bd2y_C2C(:, :, :) = ZERO
    allocate (cd2y_C2C ( nsz, 6, 6 ) ); cd2y_C2C(:, :, :) = ZERO
    allocate (dd2y_C2C ( nsz, 6, 6 ) ); dd2y_C2C(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array( nsz, domain(1)%is_periodic(i), d2fC2C, &
          ad2y_C2C, bd2y_C2C, cd2y_C2C, dd2y_C2C)
!===============================================================================
! y-direction, with np unknows
!===============================================================================
    nsz = domain(1)%np(i)
!-------------------------------------------------------------------------------
!   1st derivative in y direction with np unknows
!-------------------------------------------------------------------------------
    allocate (ad1y_P2P ( nsz, 6, 6 ) ); ad1y_P2P(:, :, :) = ZERO
    allocate (bd1y_P2P ( nsz, 6, 6 ) ); bd1y_P2P(:, :, :) = ZERO
    allocate (cd1y_P2P ( nsz, 6, 6 ) ); cd1y_P2P(:, :, :) = ZERO
    allocate (dd1y_P2P ( nsz, 6, 6 ) ); dd1y_P2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fP2P, &
          ad1y_P2P, bd1y_P2P, cd1y_P2P, dd1y_P2P)

    allocate (ad1y_C2P ( nsz, 6, 6 ) ); ad1y_C2P(:, :, :) = ZERO
    allocate (bd1y_C2P ( nsz, 6, 6 ) ); bd1y_C2P(:, :, :) = ZERO
    allocate (cd1y_C2P ( nsz, 6, 6 ) ); cd1y_C2P(:, :, :) = ZERO
    allocate (dd1y_C2P ( nsz, 6, 6 ) ); dd1y_C2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fC2P, &
          ad1y_C2P, bd1y_C2P, cd1y_C2P, dd1y_C2P) 
!-------------------------------------------------------------------------------
!   mid-point interpolation in y direction with np unknows
!-------------------------------------------------------------------------------
    allocate (am1y_C2P ( nsz, 6, 6 ) ); am1y_C2P(:, :, :) = ZERO
    allocate (bm1y_C2P ( nsz, 6, 6 ) ); bm1y_C2P(:, :, :) = ZERO
    allocate (cm1y_C2P ( nsz, 6, 6 ) ); cm1y_C2P(:, :, :) = ZERO
    allocate (dm1y_C2P ( nsz, 6, 6 ) ); dm1y_C2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), m1fC2P, &
          am1y_C2P, bm1y_C2P, cm1y_C2P, dm1y_C2P)
!-------------------------------------------------------------------------------
! 2nd order deriviative in y direction with np unknows
!-------------------------------------------------------------------------------
    allocate (ad2y_P2P ( nsz, 6, 6 ) ); ad2y_P2P(:, :, :) = ZERO
    allocate (bd2y_P2P ( nsz, 6, 6 ) ); bd2y_P2P(:, :, :) = ZERO
    allocate (cd2y_P2P ( nsz, 6, 6 ) ); cd2y_P2P(:, :, :) = ZERO
    allocate (dd2y_P2P ( nsz, 6, 6 ) ); dd2y_P2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array( nsz, domain(1)%is_periodic(i), d2fP2P, &
          ad2y_P2P, bd2y_P2P, cd2y_P2P, dd2y_P2P)
!===============================================================================
! z-direction, with nc unknows
!===============================================================================
    i = 3
    nsz = domain(1)%nc(i)
!-------------------------------------------------------------------------------
!   1st derivative in z direction with nc unknows
!-------------------------------------------------------------------------------
    allocate (ad1z_C2C ( nsz, 6, 6 ) ); ad1z_C2C(:, : ,:) = ZERO
    allocate (bd1z_C2C ( nsz, 6, 6 ) ); bd1z_C2C(:, : ,:) = ZERO
    allocate (cd1z_C2C ( nsz, 6, 6 ) ); cd1z_C2C(:, : ,:) = ZERO
    allocate (dd1z_C2C ( nsz, 6, 6 ) ); dd1z_C2C(:, : ,:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fC2C, &
          ad1z_C2C, bd1z_C2C, cd1z_C2C, dd1z_C2C)

    allocate (ad1z_P2C ( nsz, 6, 6 ) ); ad1z_P2C(:, :, :) = ZERO
    allocate (bd1z_P2C ( nsz, 6, 6 ) ); bd1z_P2C(:, :, :) = ZERO
    allocate (cd1z_P2C ( nsz, 6, 6 ) ); cd1z_P2C(:, :, :) = ZERO
    allocate (dd1z_P2C ( nsz, 6, 6 ) ); dd1z_P2C(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fP2C, &
          ad1z_P2C, bd1z_P2C, cd1z_P2C, dd1z_P2C)
!-------------------------------------------------------------------------------
!   mid-point interpolation in z direction with nc unknows
!-------------------------------------------------------------------------------
    allocate (am1z_P2C ( nsz, 6, 6 ) ); am1z_P2C(:, :, :) = ZERO
    allocate (bm1z_P2C ( nsz, 6, 6 ) ); bm1z_P2C(:, :, :) = ZERO
    allocate (cm1z_P2C ( nsz, 6, 6 ) ); cm1z_P2C(:, :, :) = ZERO
    allocate (dm1z_P2C ( nsz, 6, 6 ) ); dm1z_P2C(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array( nsz, domain(1)%is_periodic(i), m1fP2C, &
          am1z_P2C, bm1z_P2C, cm1z_P2C, dm1z_P2C)
!-------------------------------------------------------------------------------
!   2nd order deriviative in z direction with nc unknows
!-------------------------------------------------------------------------------
    allocate (ad2z_C2C ( nsz, 6, 6 ) ); ad2z_C2C(:, :, :) = ZERO
    allocate (bd2z_C2C ( nsz, 6, 6 ) ); bd2z_C2C(:, :, :) = ZERO
    allocate (cd2z_C2C ( nsz, 6, 6 ) ); cd2z_C2C(:, :, :) = ZERO
    allocate (dd2z_C2C ( nsz, 6, 6 ) ); dd2z_C2C(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array( nsz, domain(1)%is_periodic(i), d2fC2C, &
          ad2z_C2C, bd2z_C2C, cd2z_C2C, dd2z_C2C)
!===============================================================================
! z-direction, with np unknows
!===============================================================================
    nsz = domain(1)%np(i)
!-------------------------------------------------------------------------------
! 1st derivative in z direction with np unknows
!-------------------------------------------------------------------------------
    allocate (ad1z_P2P ( nsz, 6, 6 ) ); ad1z_P2P(:, :, :) = ZERO
    allocate (bd1z_P2P ( nsz, 6, 6 ) ); bd1z_P2P(:, :, :) = ZERO
    allocate (cd1z_P2P ( nsz, 6, 6 ) ); cd1z_P2P(:, :, :) = ZERO
    allocate (dd1z_P2P ( nsz, 6, 6 ) ); dd1z_P2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fP2P, &
          ad1z_P2P, bd1z_P2P, cd1z_P2P, dd1z_P2P)

    allocate (ad1z_C2P ( nsz, 6, 6 ) ); ad1z_C2P(:, :, :) = ZERO
    allocate (bd1z_C2P ( nsz, 6, 6 ) ); bd1z_C2P(:, :, :) = ZERO
    allocate (cd1z_C2P ( nsz, 6, 6 ) ); cd1z_C2P(:, :, :) = ZERO
    allocate (dd1z_C2P ( nsz, 6, 6 ) ); dd1z_C2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array(nsz, domain(1)%is_periodic(i), d1fC2P, &
          ad1z_C2P, bd1z_C2P, cd1z_C2P, dd1z_C2P)
!-------------------------------------------------------------------------------
! mid-point interpolation in z direction with np unknows
!-------------------------------------------------------------------------------
    allocate (am1z_C2P ( nsz, 6, 6 ) ); am1z_C2P(:, :, :) = ZERO
    allocate (bm1z_C2P ( nsz, 6, 6 ) ); bm1z_C2P(:, :, :) = ZERO
    allocate (cm1z_C2P ( nsz, 6, 6 ) ); cm1z_C2P(:, :, :) = ZERO
    allocate (dm1z_C2P ( nsz, 6, 6 ) ); dm1z_C2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array( nsz, domain(1)%is_periodic(i), m1fC2P, &
          am1z_C2P, bm1z_C2P, cm1z_C2P, dm1z_C2P)
!-------------------------------------------------------------------------------
! 2nd order deriviative in z direction with np unknows
!-------------------------------------------------------------------------------
    allocate (ad2z_P2P ( nsz, 6, 6 ) ); ad2z_P2P(:, :, :) = ZERO
    allocate (bd2z_P2P ( nsz, 6, 6 ) ); bd2z_P2P(:, :, :) = ZERO
    allocate (cd2z_P2P ( nsz, 6, 6 ) ); cd2z_P2P(:, :, :) = ZERO
    allocate (dd2z_P2P ( nsz, 6, 6 ) ); dd2z_P2P(:, :, :) = ZERO
    call Buildup_TDMA_LHS_array( nsz, domain(1)%is_periodic(i), d2fP2P, &
        ad2z_P2P, bd2z_P2P, cd2z_P2P, dd2z_P2P)
!===============================================================================
! x-direction
!===============================================================================
    allocate ( xtdma_lhs (nxdomain) )
    do i = 1, nxdomain
!===============================================================================
! x-direction, with nc unknows
!===============================================================================
      nsz = domain(i)%nc(1)
!-------------------------------------------------------------------------------
! 1st derivative in x direction with nc unknows
!-------------------------------------------------------------------------------
      allocate (xtdma_lhs(i)%ad1x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%ad1x_C2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bd1x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%bd1x_C2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cd1x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%cd1x_C2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dd1x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%dd1x_C2C(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array(nsz, domain(i)%is_periodic(1), d1fC2C, &
            xtdma_lhs(i)%ad1x_C2C, &
            xtdma_lhs(i)%bd1x_C2C, &
            xtdma_lhs(i)%cd1x_C2C, &
            xtdma_lhs(i)%dd1x_C2C)
  
      allocate (xtdma_lhs(i)%ad1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%ad1x_P2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bd1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%bd1x_P2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cd1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%cd1x_P2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dd1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%dd1x_P2C(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array(nsz, domain(i)%is_periodic(1), d1fP2C, &
            xtdma_lhs(i)%ad1x_P2C, &
            xtdma_lhs(i)%bd1x_P2C, &
            xtdma_lhs(i)%cd1x_P2C, &
            xtdma_lhs(i)%dd1x_P2C)
!-------------------------------------------------------------------------------
! 2nd order deriviative in x direction with nc unknows
!-------------------------------------------------------------------------------
      allocate (xtdma_lhs(i)%ad2x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%ad2x_C2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bd2x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%bd2x_C2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cd2x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%cd2x_C2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dd2x_C2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%dd2x_C2C(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array( nsz, domain(i)%is_periodic(1), d2fC2C, &
          xtdma_lhs(i)%ad2x_C2C, &
          xtdma_lhs(i)%bd2x_C2C, &
          xtdma_lhs(i)%cd2x_C2C, &
          xtdma_lhs(i)%dd2x_C2C)
!-------------------------------------------------------------------------------
! mid-point interpolation in x direction with nc unknows
!-------------------------------------------------------------------------------
      allocate (xtdma_lhs(i)%am1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%am1x_P2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bm1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%bm1x_P2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cm1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%cm1x_P2C(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dm1x_P2C ( nsz, 6, 6 ) ); xtdma_lhs(i)%dm1x_P2C(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array(nsz, domain(i)%is_periodic(1), m1fP2C, &
          xtdma_lhs(i)%am1x_P2C, &
          xtdma_lhs(i)%bm1x_P2C, &
          xtdma_lhs(i)%cm1x_P2C, &
          xtdma_lhs(i)%dm1x_P2C)      
!===============================================================================
! x-direction, with np unknows
!===============================================================================
      nsz = domain(i)%np(1)
!-------------------------------------------------------------------------------
! 1st derivative in x direction with np unknows
!-------------------------------------------------------------------------------
      allocate (xtdma_lhs(i)%ad1x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%ad1x_P2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bd1x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%bd1x_P2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cd1x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%cd1x_P2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dd1x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%dd1x_P2P(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array(nsz, domain(i)%is_periodic(1), d1fP2P, &
            xtdma_lhs(i)%ad1x_P2P, &
            xtdma_lhs(i)%bd1x_P2P, &
            xtdma_lhs(i)%cd1x_P2P, &
            xtdma_lhs(i)%dd1x_P2P)
  
      allocate (xtdma_lhs(i)%ad1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%ad1x_C2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bd1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%bd1x_C2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cd1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%cd1x_C2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dd1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%dd1x_C2P(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array(nsz, domain(i)%is_periodic(1), d1fC2P, &
            xtdma_lhs(i)%ad1x_C2P, &
            xtdma_lhs(i)%bd1x_C2P, &
            xtdma_lhs(i)%cd1x_C2P, &
            xtdma_lhs(i)%dd1x_C2P)
!-------------------------------------------------------------------------------
! 2nd order deriviative in x direction with np unknows
!-------------------------------------------------------------------------------
      allocate (xtdma_lhs(i)%ad2x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%ad2x_P2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bd2x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%bd2x_P2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cd2x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%cd2x_P2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dd2x_P2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%dd2x_P2P(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array( nsz, domain(i)%is_periodic(1), d2fP2P, &
            xtdma_lhs(i)%ad2x_P2P, &
            xtdma_lhs(i)%bd2x_P2P, &
            xtdma_lhs(i)%cd2x_P2P, &
            xtdma_lhs(i)%dd2x_P2P)
!-------------------------------------------------------------------------------
! mid-point interpolation in x direction with np unknows
!-------------------------------------------------------------------------------
      allocate (xtdma_lhs(i)%am1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%am1x_C2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%bm1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%bm1x_C2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%cm1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%cm1x_C2P(:, :, :) = ZERO
      allocate (xtdma_lhs(i)%dm1x_C2P ( nsz, 6, 6 ) ); xtdma_lhs(i)%dm1x_C2P(:, :, :) = ZERO
      call Buildup_TDMA_LHS_array(nsz, domain(i)%is_periodic(1), m1fC2P, &
          xtdma_lhs(i)%am1x_C2P, &
          xtdma_lhs(i)%bm1x_C2P, &
          xtdma_lhs(i)%cm1x_C2P, &
          xtdma_lhs(i)%dm1x_C2P)      

    end do

    return
  end subroutine Prepare_LHS_coeffs_for_operations
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for interpolation.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the interpolation.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       given pencil     needed         specified   private
!-------------------------------------------------------------------------------
!>  index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! interpolation. P2C: periodic
! [ 1     alpha                       alpha] [f_1] = [a/2 * ( f_{1'}   + f_{2'}   ) + b/2 * ( f_{3'}   + f_{n'}   ) ]
! [       alpha  1      alpha              ] [f_2]   [a/2 * ( f_{2'}   + f_{3'}   ) + b/2 * ( f_{4'}   + f_{1'}   ) ]
! [              alpha  1      alpha       ] [f_i]   [a/2 * ( f_{i'}   + f_{i'+1} ) + b/2 * ( f_{i'+2} + f_{i'-1} ) ]
! [                     alpha  1      alpha] [f_4]   [a/2 * ( f_{n'-1} + f_{n'}   ) + b/2 * ( f_{1'}   + f_{n'-2} ) ]
! [alpha                       alpha  1    ] [f_5]   [a/2 * ( f_{n'}   + f_{1'}   ) + b/2 * ( f_{2'}   + f_{n'-1} ) ]
!-------------------------------------------------------------------------------
! interpolation. P2C: Dirichlet
! [ 1     alpha1                              ] [f_1] = [a1 *  f_{1'}   + b1 * f_{2'} + c1 * f_{3'}                       ]
! [alpha2 1      alpha2                       ] [f_2]   [a2/2 * ( f_{2'}   + f_{3'}   )                                 ] 
! [              alpha  1       alpha         ] [f_i]   [a/2 * ( f_{i'}   + f_{i'+1} ) + b/2 * ( f_{i'+2} + f_{i'-1} ) ]
! [                     alpha2  1       alpha2] [f_4]   [a2/2 * ( f_{n'-1} + f_{n'}   )                                 ]
! [                             alpha1  1     ] [f_5]   [a1 *  f_{n'+1}   + b1 * f_{n'}  + c1 * f_{n'-1}                  ]
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is nc
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_interp_P2C_RHS_array(fi, fo, nc, coeff, ibc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: nc ! unknow numbers, nc
    real(WP), intent(out) :: fo(nc)
    real(WP), intent(in ) :: coeff(5, 3, 6)
    integer,  intent(in ) :: ibc(2)

    integer :: i, m, l

    fo(:) = ZERO
!-------------------------------------------------------------------------------
!   i = bulk
!-------------------------------------------------------------------------------
    do i = 3, nc - 2
      fo(i) = coeff( 3, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( 3, 2, IBC_PERIODIC ) * ( fi(i - 1) + fi(i + 2) )
    end do
!-------------------------------------------------------------------------------
!   for either perirodic or non-periodic of nc
!   ---(-1')-(-1)-(0')-(0)-(|1')-(1)-(2')-(2)-(3')---(i')---(nc'-1)-(nc-1)-(nc')-(nc)-(nc'+1)|-(nc+1)-(nc'+2)---
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   i = 1
!-------------------------------------------------------------------------------
    i = 1
    m = 1
    l = 1
    if ( ibc(m) == IBC_PERIODIC) then
      ! 0' = nc = np
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(nc   ) + fi(i + 2) )
    else if (ibc(m) == IBC_SYMMETRIC ) then
      ! 0' = 2'
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(2    ) + fi(i + 2) )
    else if (ibc(m) == IBC_ASYMMETRIC) then
      ! 0' = 2'
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( l, 2, ibc(m) ) * (-fi(2    ) + fi(i + 2) )
    else
      fo(i) = coeff( l, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( l, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( l, 3, IBC_INTRPL) * fi(i + 2) 
    end if
!-------------------------------------------------------------------------------
!   i = 2
!-------------------------------------------------------------------------------
    i = 2
    l = 2
    fo(i) = coeff( l, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i + 1) ) + &
            coeff( l, 2, IBC_PERIODIC ) * ( fi(i - 1) + fi(i + 2) )
!-------------------------------------------------------------------------------
!   i = nc - 1
!-------------------------------------------------------------------------------
    i = nc - 1
    m = 2
    l = 4
    if ( ibc(m) == IBC_PERIODIC) then
      ! i + 2 = nc' + 1 = 1
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i - 1) + fi(1    ) )
    else
      fo(i) = coeff( l, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( l, 2, IBC_PERIODIC ) * ( fi(i - 1) + fi(i + 2) )
    end if
!-------------------------------------------------------------------------------
!   i = nc
!-------------------------------------------------------------------------------
    i = nc
    m = 2
    l = 5
    if ( ibc(m) == IBC_PERIODIC) then
      ! i + 1 = nc' + 1 = 1
      ! i + 2 = nc' + 2 = 2
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i - 1) + fi(2) )
    else if (ibc(m) == IBC_SYMMETRIC ) then
      ! i + 2 = nc' + 2 = nc'
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i - 1) + fi(nc   ) )
    else if (ibc(m) == IBC_ASYMMETRIC) then
      ! i + 2 = nc' + 2 = nc'
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i - 1) - fi(nc   ) )
    else
      fo(i) = coeff( l, 1, IBC_INTRPL) * fi(i + 1) + &
              coeff( l, 2, IBC_INTRPL) * fi(i    ) + &
              coeff( l, 3, IBC_INTRPL) * fi(i - 1) 
    end if
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    ! nothing.
    return
  end subroutine Prepare_TDMA_interp_P2C_RHS_array

!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for interpolation.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the interpolation.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   private
!-------------------------------------------------------------------------------
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! interpolation. C2P: periodic
! [ 1    alpha                   alpha][f_1']=[a/2 * (f_{1}   + f_{n})   + b/2 * (f_{2}   + f_{n-1})]
! [      alpha 1     alpha            ][f_2'] [a/2 * (f_{2}   + f_{1})   + b/2 * (f_{3}   + f_{n})  ]
! [            alpha 1     alpha      ][f_i'] [a/2 * (f_{i}   + f_{i-1}) + b/2 * (f_{i+1} + f_{i-2})]
! [                  alpha 1     alpha][f_4'] [a/2 * (f_{n-1} + f_{n-2}) + b/2 * (f_{n}   + f_{n-3})]
! [alpha                   alpha 1    ][f_5'] [a/2 * (f_{1}   + f_{n-2}) + b/2 * (f_{n}   + f_{n-1})]
!-------------------------------------------------------------------------------
! interpolation. C2P: Dirichlet
! [ 1    alpha1                          ][f_1']=[a1 * f_{1} + b1 * f_{2} + c1 * f_{3}  ]
! [      alpha2 1     alpha2             ][f_2'] [a2/2 * (f_{2}   + f_{1})]
! [             alpha 1      alpha       ][f_i'] [a/2 * (f_{i}   + f_{i-1}) + b/2 * (f_{i+1} + f_{i-2})]
! [                   alpha2 1     alpha2][f_4'] [a2/2 * (f_{n-1}   + f_{n-2})]
! [                          alpha1 1    ][f_5'] [a1 * f_{n-1} + b1 * f_{n-2} + c1 * f_{n-3}]
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is np
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_interp_C2P_RHS_array(fi, fo, np, coeff, ibc, fbc)
    use parameters_constant_mod
    implicit none

    real(WP),           intent(in ) :: fi(:)
    integer,            intent(in ) :: np ! unknow numbers, np
    real(WP),           intent(out) :: fo(np)
    real(WP),           intent(in ) :: coeff(5, 3, 6)
    integer,            intent(in ) :: ibc(2)
    real(WP), optional, intent(in)  :: fbc(2) ! used for Dirichlet B.C.

    integer :: i, m, l

    fo(:) = ZERO
!-------------------------------------------------------------------------------
!   i = bulk
!-------------------------------------------------------------------------------
    do i = 3, np - 2
      fo(i) = coeff( 3, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( 3, 2, IBC_PERIODIC ) * ( fi(i + 1) + fi(i - 2) )
    end do
!-------------------------------------------------------------------------------
!   for non-periodic:
!   ---(-1')-(-1)-(0')-(0)-(|1')-(1)-(2')-(2)-(3')---(i')---(np'-1)-(np-1)-(np'|)-(np)-(np'+1)-(np+1)-(np'+2)---
!   for periodic:
!   ---(-1')-(-1)-(0')-(0)-(|1')-(1)-(2')-(2)-(3')---(i')---(np'-1)-(np-1)-(np')-(np)-(np'+1|)-(np+1)-(np'+2)---
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   i = 1
!-------------------------------------------------------------------------------
    i = 1
    m = 1
    l = 1
    if ( ibc(m) == IBC_PERIODIC) then
      ! 0 = np
      !-1 = np - 1
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(np    ) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i + 1) + fi(np - 1) )
    else if (ibc(m) == IBC_SYMMETRIC ) then
      ! 0 = 1
      !-1 = 2
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i + 1) + fi(2) )
    else if (ibc(m) == IBC_ASYMMETRIC) then
      fo(i) = ZERO
    else if (ibc(m) == IBC_DIRICHLET) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_DIRICHLET')
      fo(i) = fbc(m)
    else
      fo(i) = coeff( l, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( l, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( l, 3, IBC_INTRPL) * fi(i + 2) 
    end if
!-------------------------------------------------------------------------------
!   i = 2
!-------------------------------------------------------------------------------
    i = 2
    m = 1
    l = 2
    if ( ibc(m) == IBC_PERIODIC) then
      ! 0 = np
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i + 1) + fi(np   ) )
    else if (ibc(m) == IBC_SYMMETRIC ) then
      ! 0 = 1
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i + 1) + fi(1    ) )
    else if (ibc(m) == IBC_ASYMMETRIC) then
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i + 1) - fi(1    ) )
    else if (ibc(m) == IBC_DIRICHLET) then
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i - 1) )
    else
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i - 1) )
    end if
!-------------------------------------------------------------------------------
!   i = np - 1
!-------------------------------------------------------------------------------
    i = np - 1
    m = 2
    l = 4
    if ( ibc(m) == IBC_PERIODIC) then
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i     ) + fi(i - 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(i + 1 ) + fi(i - 2 ) )
    else if (ibc(m) == IBC_SYMMETRIC ) then
      ! np     = np - 1
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i     ) + fi(i - 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(np - 1) + fi(i - 2) )
    else if (ibc(m) == IBC_ASYMMETRIC) then
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i     ) + fi(i - 1) ) + &
              coeff( l, 2, ibc(m) ) * (-fi(np - 1) + fi(i - 2) )
    else if (ibc(m) == IBC_DIRICHLET) then
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i     ) + fi(i - 1) )
    else
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i     ) + fi(i - 1) )
    end if
!-------------------------------------------------------------------------------
!   i = np
!-------------------------------------------------------------------------------
    i = np
    m = 2
    l = 5
    if ( ibc(m) == IBC_PERIODIC) then
      ! np + 1 = 1
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i    ) + fi(i - 1 ) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(1    ) + fi(i - 2) )
    else if (ibc(m) == IBC_SYMMETRIC ) then
      ! np     = np - 1
      ! np + 1 = np - 2
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(np - 1) + fi(i - 1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(np - 2) + fi(i - 2) )
    else if (ibc(m) == IBC_ASYMMETRIC) then
      fo(i) = ZERO
    else if (ibc(m) == IBC_DIRICHLET) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_DIRICHLET')
      fo(i) = fbc(m)
    else
      fo(i) = coeff( l, 1, IBC_INTRPL) * fi(i - 1) + &
              coeff( l, 2, IBC_INTRPL) * fi(i - 2) + &
              coeff( l, 3, IBC_INTRPL) * fi(i - 3) 
    end if
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    ! nothing
    return
  end subroutine Prepare_TDMA_interp_C2P_RHS_array
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   private
!-------------------------------------------------------------------------------
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! 1st derivative on collocated grids, C2C/P2P coefficients : Periodic or Symmetric B.C.
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!-------------------------------------------------------------------------------
! 1st-derivative. C2C or P2P: periodic
! [ 1    alpha                   alpha][f'_1]=[a/2 * (f_{2}   - f_{n})/h   + b/4 * (f_{3}   - f_{n-1})/h]
! [      alpha 1     alpha            ][f'_2] [a/2 * (f_{3}   - f_{1})/h   + b/4 * (f_{4}   - f_{n})/h  ]
! [            alpha 1     alpha      ][f'_i] [a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                  alpha 1     alpha][f'_4] [a/2 * (f_{n}   - f_{n-2})/h + b/4 * (f_{1}   - f_{n-3})/h]
! [alpha                   alpha 1    ][f'_5] [a/2 * (f_{1}   - f_{n-1})/h + b/4 * (f_{2}   - f_{n-2})/h]
!-------------------------------------------------------------------------------
! 1st-derivative : C2C or P2P : Dirichlet B.C.
! [ 1     alpha1                            ][f'_1]=[a1 * f_{1}/h  + b1 * f_{2}/h + c1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2] [a2/2 * (f_{3} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i] [ a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4] [a2/2 * (f_{n} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5] [-a1 * f_{n}/h  - b1 * f_{n-1}/h - c1 * f_{n-2}/h]
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is np
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     dd            1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_1deri_C2C_RHS_array(fi, fo, n, inbr, coeff, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: n ! unknow numbers
    real(WP), intent(out) :: fo(n)
    integer,  intent(in ) :: inbr(4, 4)
    real(WP), intent(in ) :: coeff(5, 4, 6)
    real(WP), intent(in ) :: dd
    integer,  intent(in ) :: ibc(2)
    real(WP), optional, intent(in)  :: fbc(2) ! used for Dirichlet B.C.

    integer :: i, j, m, l
    integer :: im2, im1, ip1, ip2

    fo(:) = ZERO

!-------------------------------------------------------------------------------
!   boundaries
!-------------------------------------------------------------------------------
    do j = 1, 4
      ! preparation
      if(j < 3) then
        m = 1
        l = j
        i = j
      else if (j == 3) then
        m = 2
        l = 4
        i = n - 1
      else if (j == 4) then
        m = 2
        l = 5
        i = n
      else
      end if
      im2 = inbr(j, 1)
      im1 = inbr(j, 2)
      ip1 = inbr(j, 3)
      ip2 = inbr(j, 4)

      ! assign to different bc
      if ( ibc(m) == IBC_PERIODIC  .or. &
           ibc(m) == IBC_SYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( fi(ip2) - fi(im2) )
      else if ( ibc(m) == IBC_ASYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( SIGN(ONE, real(n - i - 1, wp) ) * fi(ip1) - &
                                          SIGN(ONE, real(i - 1 - 1, wp) ) * fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( SIGN(ONE, real(n - i - 2, wp) ) * fi(ip2) - &
                                          SIGN(ONE, real(i - 2 - 1, wp) ) * fi(im2) )
      else if (ibc(m) == IBC_DIRICHLET) then
        if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_DIRICHLET')
        if(i == 1 .or. i == n) then
          fo(i) = coeff( l, 1, ibc(m) ) * fbc(m)                + &
                  coeff( l, 2, ibc(m) ) * fi(i)                 + &
                  coeff( l, 3, ibc(m) ) * fi(i + SIGN(1, 1-m))  + &
                  coeff( l, 4, ibc(m) ) * fi(i + SIGN(2, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - fi(im1) )
        else
        end if
      else! all other b.c. = interpolation = no specified bc = neumann
        if(i == 1 .or. i == n) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * fi(i)                 + &
                  coeff( l, 2, IBC_INTRPL ) * fi(i + SIGN(1, 1-m))  + &
                  coeff( l, 3, IBC_INTRPL ) * fi(i + SIGN(2, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * ( fi(ip1) - fi(im1) )
        else
        end if
      end if
      
    end do
!-------------------------------------------------------------------------------
!   bulk area
!-------------------------------------------------------------------------------
    l = 3
    do i = 3, n - 2
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2
      fo(i) = coeff( l, 1, IBC_PERIODIC ) * ( fi(ip1) - fi(im1) ) + &
              coeff( l, 2, IBC_PERIODIC ) * ( fi(ip2) - fi(im2) )
    end do
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_C2C_RHS_array
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   private
!-------------------------------------------------------------------------------
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! 1st derivative on collocated grids, C2C/P2P coefficients : Periodic or Symmetric B.C.
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!-------------------------------------------------------------------------------
! 1st-derivative. C2C or P2P: periodic
! [ 1    alpha                   alpha][f'_1]=[a/2 * (f_{2}   - f_{n})/h   + b/4 * (f_{3}   - f_{n-1})/h]
! [      alpha 1     alpha            ][f'_2] [a/2 * (f_{3}   - f_{1})/h   + b/4 * (f_{4}   - f_{n})/h  ]
! [            alpha 1     alpha      ][f'_i] [a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                  alpha 1     alpha][f'_4] [a/2 * (f_{n}   - f_{n-2})/h + b/4 * (f_{1}   - f_{n-3})/h]
! [alpha                   alpha 1    ][f'_5] [a/2 * (f_{1}   - f_{n-1})/h + b/4 * (f_{2}   - f_{n-2})/h]
!-------------------------------------------------------------------------------
! 1st-derivative : C2C or P2P : Dirichlet B.C.
! [ 1     alpha1                            ][f'_1]=[a1 * f_{1}/h  + b1 * f_{2}/h + c1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2] [a2/2 * (f_{3} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i] [ a/2 * (f_{i+1} - f_{i-1})/h + b/4 * (f_{i+2} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4] [a2/2 * (f_{n} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5] [-a1 * f_{n}/h  - b1 * f_{n-1}/h - c1 * f_{n-2}/h]
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is np
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     dd            1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_1deri_P2P_RHS_array(fi, fo, n, inbr, coeff, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP),           intent(in ) :: fi(:)
    integer,            intent(in ) :: n ! unknow numbers
    real(WP),           intent(out) :: fo(n)
    integer,            intent(in ) :: inbr(4, 4)
    real(WP),           intent(in ) :: coeff(5, 3, 6)
    real(WP),           intent(in ) :: dd
    integer,            intent(in ) :: ibc(2)
    real(WP), optional, intent(in ) :: fbc(2) ! used for IBC_NEUMANN

    integer :: m, l, i, j
    integer :: im1, im2, ip1, ip2

    fo(:) = ZERO
!-------------------------------------------------------------------------------
!   boundaries
!-------------------------------------------------------------------------------
    do j = 1, 4
      ! preparation
      if(j < 3) then
        m = 1
        l = j
        i = j
      else if (j == 3) then
        m = 2
        l = 4
        i = n - 1
      else if (j == 4) then
        m = 2
        l = 5
        i = n
      else
      end if
      im2 = inbr(j, 1)
      im1 = inbr(j, 2)
      ip1 = inbr(j, 3)
      ip2 = inbr(j, 4)

      ! assign to different bc
      if ( ibc(m) == IBC_PERIODIC  .or. &
           ibc(m) == IBC_SYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( fi(ip2) - fi(im2) )
      else if ( ibc(m) == IBC_ASYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( SIGN(ONE, real(n - i - 1, wp) ) * fi(ip1) - &
                                          SIGN(ONE, real(i - 1 - 1, wp) ) * fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( SIGN(ONE, real(n - i - 2, wp) ) * fi(ip2) - &
                                          SIGN(ONE, real(i - 2 - 1, wp) ) * fi(im2) )
      else if (ibc(m) == IBC_NEUMANN) then
        if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_NEUMANN')
        if(i == 1 .or. i == n) then
          fo(i) = fbc(m)
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - fi(im1) )
        else
        end if
      else! all other b.c. = interpolation = no specified bc = dirichlet
        if(i == 1 .or. i == n) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * fi(i)                 + &
                  coeff( l, 2, IBC_INTRPL ) * fi(i + SIGN(1, 1-m))  + &
                  coeff( l, 3, IBC_INTRPL ) * fi(i + SIGN(2, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * ( fi(ip1) - fi(im1) )
        else
        end if
      end if
      
    end do
!-------------------------------------------------------------------------------
!   bulk area
!-------------------------------------------------------------------------------
    l = 3
    do i = 3, n - 2
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2
      fo(i) = coeff( l, 1, IBC_PERIODIC ) * ( fi(ip1) - fi(im1) ) + &
              coeff( l, 2, IBC_PERIODIC ) * ( fi(ip2) - fi(im2) )
    end do
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_P2P_RHS_array
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   private
!-------------------------------------------------------------------------------
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! 1st derivative on staggered grids C2P
! C2P ==>
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(h ) * ( f_{i}   - f_{i-1} ) + &
!                                                 b/(3h) * ( f_{i+1} - f_{i-2} )
!-------------------------------------------------------------------------------
! 1st-derivative. C2P: periodic
! [ 1    alpha                   alpha][f_1']=[a/2 * (f_{1}   + f_{n})   + b/2 * (f_{2}   + f_{n-1})]
! [      alpha 1     alpha            ][f_2'] [a/2 * (f_{2}   + f_{1})   + b/2 * (f_{3}   + f_{n})  ]
! [            alpha 1     alpha      ][f_i'] [a/2 * (f_{i}   + f_{i-1}) + b/2 * (f_{i+1} + f_{i-2})]
! [                  alpha 1     alpha][f_4'] [a/2 * (f_{n-1} + f_{n-2}) + b/2 * (f_{n}   + f_{n-3})]
! [alpha                   alpha 1    ][f_5'] [a/2 * (f_{1}   + f_{n-2}) + b/2 * (f_{n}   + f_{n-1})]
!-------------------------------------------------------------------------------
! 1st-derivative : C2P : Dirichlet B.C.
! [ 1     alpha1                            ][f'_1']=[a1 * f_{1}/h  + b1 * f_{2}/h + c1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2'] [a2 * (f_{2} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i'] [a *  (f_{i} - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4'] [a2 * (f_{n-1} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5'] [-a1 * f_{n-1}/h  - b1 * f_{n-2}/h - c1 * f_{n-3}/h]
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is np
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     dd            1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_1deri_C2P_RHS_array(fi, fo, n, inbr, coeff, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP),           intent(in ) :: fi(:)
    integer,            intent(in ) :: n ! unknow numbers, np
    real(WP),           intent(out) :: fo(n)
    integer,            intent(in ) :: inbr(4, 4)
    real(WP),           intent(in ) :: coeff(5, 4, 6)
    real(WP),           intent(in ) :: dd
    integer,            intent(in ) :: ibc(2)
    real(WP), optional, intent(in ) :: fbc(2) ! used for IBC_NEUMANN

    integer :: i, j, k, m, l, ncmax
    integer :: im1, im2, ip1

    fo(:) = ZERO

!-------------------------------------------------------------------------------
!   boundaries
!-------------------------------------------------------------------------------
    if(ibc(1) == IBC_PERIODIC .or. ibc(2) == IBC_PERIODIC) then
      ncmax = n
    else
      ncmax = n - 1
    end if

    do j = 1, 4
      ! preparation
      if(j < 3) then
        m = 1
        l = j
        i = j
      else if (j == 3) then
        m = 2
        l = 4
        i = n - 1
      else if (j == 4) then
        m = 2
        l = 5
        i = n
      else
      end if
      im2 = inbr(j, 1)
      im1 = inbr(j, 2)
      ip1 = inbr(j, 3)
      !ip2 = inbr(j, 4)

      ! assign to different bc
      if ( ibc(m) == IBC_PERIODIC  .or. &
           ibc(m) == IBC_SYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i  ) - fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( fi(ip1) - fi(im2) )
      else if ( ibc(m) == IBC_ASYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * (                                       fi(i  ) - &
                                          SIGN(ONE, real(i     - 1 - 1, wp) ) * fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( SIGN(ONE, real(ncmax - i - 1, wp) ) * fi(ip1) - &
                                          SIGN(ONE, real(i     - 2 - 1, wp) ) * fi(im2) )
      else if (ibc(m) == IBC_NEUMANN) then
        if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_NEUMANN')
        if(i == 1 .or. i == n) then
          fo(i) = fbc(m)
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i) - fi(im1) )
        else
        end if
      else if (ibc(m) == IBC_DIRICHLET) then
        if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_DIRICHLET')
        if (i == 1) k = i
        if (i == n) k = ncmax
        if(i == 1 .or. i == n) then
          fo(i) = coeff( l, 1, ibc(m) ) * fbc(m)                + &
                  coeff( l, 2, ibc(m) ) * fi(k)                 + &
                  coeff( l, 3, ibc(m) ) * fi(k + SIGN(1, 1-m))  + &
                  coeff( l, 4, ibc(m) ) * fi(k + SIGN(2, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, ibc(m) ) * ( fi(i) - fi(im1) )
        else
        end if
      else! all other b.c. = interpolation = no specified bc = dirichlet
        if(i == 1 .or. i == n) then
          if (i == 1) k = i
          if (i == n) k = ncmax
          fo(i) = coeff( l, 1, IBC_INTRPL ) * fi(k)                 + &
                  coeff( l, 2, IBC_INTRPL ) * fi(k + SIGN(1, 1-m))  + &
                  coeff( l, 3, IBC_INTRPL ) * fi(k + SIGN(2, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * ( fi(i) - fi(im1) )
        else
        end if
      end if
      
    end do
!-------------------------------------------------------------------------------
!   bulk area
!-------------------------------------------------------------------------------
    l = 3
    do i = 3, n - 2
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      !ip2 = i + 2
      fo(i) = coeff( l, 1, ibc(1) ) * ( fi(i)   - fi(im1) ) + &
              coeff( l, 2, ibc(1) ) * ( fi(ip1) - fi(im2) )
    end do
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    fo(:) = fo(:) * dd
    return
  end subroutine Prepare_TDMA_1deri_C2P_RHS_array

!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   private
!-------------------------------------------------------------------------------
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
! 1st derivative on staggered grids P2C and C2P : Periodic or Symmetric B.C.
! P2C ==>
! alpha * f'_{i-1} +  f'_i +  alpha * f'_{i+1}  = a/(h ) * ( f_{i'+1} - f_{i'} ) + &
!                                                 b/(3h) * ( f_{i'+2} - f_{i'-1} )
!-------------------------------------------------------------------------------
! 1st-derivative. P2C: Dirichlet
! [ 1     alpha1                            ][f'_1]=[a1 * f_{1'}/h  + b1 * f_{2'}/h + c1 * f_{3'}/h  ]
! [alpha2 1      alpha2                     ][f'_2] [a2 * (f_{3'} - f_{2'})/h  ]
! [       alpha  1      alpha               ][f'_i] [a *  (f_{i'+1} - f_{i'})/h + b/3 * (f_{i'+2} - f_{i'-1})/h]
! [                     alpha2 1      alpha2][f'_4] [a2 * (f_{n'} - f_{n'-1})/h]
! [                            alpha1 1     ][f'_5] [-a1 * f_{n'+1}/h  - b1 * f_{n'}/h - c1 * f_{n'-1}/h]
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is np
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     dd            1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_1deri_P2C_RHS_array(fi, fo, n, inbr, coeff, dd, ibc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: n ! unknow numbers, nc
    real(WP), intent(out) :: fo(n)
    integer,  intent(in ) :: inbr(4, 4)
    real(WP), intent(in ) :: coeff(5, 3, 6)
    real(WP), intent(in ) :: dd
    integer,  intent(in ) :: ibc(2)

    integer :: i, j, m, l, k, npmax
    integer :: im1, ip1, ip2

    fo(:) = ZERO

!-------------------------------------------------------------------------------
!   boundaries
!-------------------------------------------------------------------------------
    if(ibc(1) == IBC_PERIODIC .or. ibc(2) == IBC_PERIODIC) then
      npmax = n
    else
      npmax = n + 1
    end if

    do j = 1, 4
      ! preparation
      if(j < 3) then
        m = 1
        l = j
        i = j
      else if (j == 3) then
        m = 2
        l = 4
        i = n - 1
      else if (j == 4) then
        m = 2
        l = 5
        i = n
      else
      end if
      !im2 = inbr(j, 1)
      im1 = inbr(j, 2)
      ip1 = inbr(j, 3)
      ip2 = inbr(j, 4)

      ! assign to different bc
      if ( ibc(m) == IBC_PERIODIC  .or. &
           ibc(m) == IBC_SYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - fi(i  ) ) + &
                coeff( l, 2, ibc(m) ) * ( fi(ip2) - fi(im1) )
      else if ( ibc(m) == IBC_ASYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( SIGN(ONE, real(npmax - i - 1, wp) ) * fi(ip1) - &
                                                                                fi(i  ) ) + &
                coeff( l, 2, ibc(m) ) * ( SIGN(ONE, real(npmax - i - 2, wp) ) * fi(ip2) - &
                                          SIGN(ONE, real(i     - 1 - 1, wp) ) * fi(im1) )
      else! all other b.c. = interpolation = no specified bc = dirichlet
        if(i == 1 .or. i == n) then
          if (i == 1) k = i
          if (i == n) k = npmax
          fo(i) = coeff( l, 1, IBC_INTRPL ) * fi(k)                 + &
                  coeff( l, 2, IBC_INTRPL ) * fi(k + SIGN(1, 1-m))  + &
                  coeff( l, 3, IBC_INTRPL ) * fi(k + SIGN(2, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * ( fi(ip1) - fi(i  ) ) + &
                  coeff( l, 2, IBC_INTRPL ) * ( fi(ip2) - fi(im1) )
        else
        end if
      end if
      
    end do
!-------------------------------------------------------------------------------
!   bulk area
!-------------------------------------------------------------------------------
    l = 3
    do i = 3, n - 2
      !im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2
      fo(i) = coeff( l, 1, IBC_PERIODIC ) * ( fi(ip1) - fi(i  ) ) + &
              coeff( l, 2, IBC_PERIODIC ) * ( fi(ip2) - fi(im1) )
    end do
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_P2C_RHS_array
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   private
!-------------------------------------------------------------------------------
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! 2nd derivative on collocated grids, C2C/P2P coefficients : Periodic or Symmetric B.C.
! alpha * f"_{i-1} + f"_i + alpha * f"_{i+1} = a/(2h) * ( f_{i+1} - 2f_{i} + f(i-1) ) + &
!                                              b/(4h) * ( f_{i+2} - 2f_{i} + f_{i-2} )
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is nc
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     dd            1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_2deri_C2C_RHS_array(fi, fo, n, inbr, coeff, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: n ! unknow numbers
    real(WP), intent(out) :: fo(n)
    integer,  intent(in ) :: inbr(4, 4)
    real(WP), intent(in ) :: coeff(5, 4, 6)
    real(WP), intent(in ) :: dd
    integer,  intent(in ) :: ibc(2)
    real(WP), intent(in ), optional :: fbc(2)

    integer :: i, l, m, j
    integer :: im2, im1, ip1, ip2

    fo(:) = ZERO
!-------------------------------------------------------------------------------
!   boundaries
!-------------------------------------------------------------------------------
    do j = 1, 4
      ! preparation
      if(j < 3) then
        m = 1
        l = j
        i = j
      else if (j == 3) then
        m = 2
        l = 4
        i = n - 1
      else if (j == 4) then
        m = 2
        l = 5
        i = n
      else
      end if
      im2 = inbr(j, 1)
      im1 = inbr(j, 2)
      ip1 = inbr(j, 3)
      ip2 = inbr(j, 4)

      ! assign to different bc
      if ( ibc(m) == IBC_PERIODIC  .or. &
           ibc(m) == IBC_SYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - TWO * fi(i) + fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( fi(ip2) - TWO * fi(i) + fi(im2) )
      else if ( ibc(m) == IBC_ASYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( SIGN(ONE, real(n - i - 1, wp) ) * fi(ip1) - &
                                                                      TWO * fi(i  ) + &
                                          SIGN(ONE, real(i - 1 - 1, wp) ) * fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( SIGN(ONE, real(n - i - 2, wp) ) * fi(ip2) - &
                                                                      TWO * fi(i  ) + &
                                          SIGN(ONE, real(i - 2 - 1, wp) ) * fi(im2) )
      else if (ibc(m) == IBC_DIRICHLET) then
        if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_DIRICHLET')
        if(i == 1 .or. i == n) then
          fo(i) = coeff( l, 1, ibc(m) ) * fbc(m)                + &
                  coeff( l, 2, ibc(m) ) * fi(i)                 + &
                  coeff( l, 3, ibc(m) ) * fi(i + SIGN(1, 1-m))  + &
                  coeff( l, 4, ibc(m) ) * fi(i + SIGN(2, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - TWO * fi(i) + fi(im1) )
        else
        end if
      else! all other b.c. = interpolation = no specified bc = neumann
        if(i == 1 .or. i == n) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * fi(i)                 + &
                  coeff( l, 2, IBC_INTRPL ) * fi(i + SIGN(1, 1-m))  + &
                  coeff( l, 3, IBC_INTRPL ) * fi(i + SIGN(2, 1-m))  + &
                  coeff( l, 4, IBC_INTRPL ) * fi(i + SIGN(3, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * ( fi(ip1) - TWO * fi(i) + fi(im1) )
        else
        end if
      end if
      
    end do
!-------------------------------------------------------------------------------
!   bulk area
!-------------------------------------------------------------------------------
    l = 3
    do i = 3, n - 2
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - TWO * fi(i) + fi(im1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(ip2) - TWO * fi(i) + fi(im2) )
    end do
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    fo(:) = fo(:) * dd ! dd = (1/dx)^2

    return
  end subroutine Prepare_TDMA_2deri_C2C_RHS_array
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   private
!-------------------------------------------------------------------------------
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! 2nd derivative on collocated grids, C2C/P2P coefficients : Periodic or Symmetric B.C.
! alpha * f"_{i-1} + f"_i + alpha * f"_{i+1} = a/(2h) * ( f_{i+1} - 2f_{i} + f(i-1) ) + &
!                                              b/(4h) * ( f_{i+2} - 2f_{i} + f_{i-2} )
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknowns, here is nc
!> \param[in]    ibc            the b.c. type at two ends of the unknown array
!> \param[in]    fbc            the b.c. values for the given ibc
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     dd            1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!===============================================================================
  subroutine Prepare_TDMA_2deri_P2P_RHS_array(fi, fo, n, inbr, coeff, dd, ibc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: n ! unknow numbers
    real(WP), intent(out) :: fo(n)
    integer,  intent(in ) :: inbr(4, 4)
    real(WP), intent(in ) :: coeff(5, 4, 6)
    real(WP), intent(in ) :: dd
    integer,  intent(in ) :: ibc(2)

    integer  :: im2, im1, ip1, ip2
    integer  :: i, j, l, m

    fo(:) = ZERO

!-------------------------------------------------------------------------------
!   boundaries
!-------------------------------------------------------------------------------
    do j = 1, 4
      ! preparation
      if(j < 3) then
        m = 1
        l = j
        i = j
      else if (j == 3) then
        m = 2
        l = 4
        i = n - 1
      else if (j == 4) then
        m = 2
        l = 5
        i = n
      else
      end if
      im2 = inbr(j, 1)
      im1 = inbr(j, 2)
      ip1 = inbr(j, 3)
      ip2 = inbr(j, 4)

      ! assign to different bc
      if ( ibc(m) == IBC_PERIODIC  .or. &
           ibc(m) == IBC_SYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - TWO * fi(i) + fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( fi(ip2) - TWO * fi(i) + fi(im2) )
      else if ( ibc(m) == IBC_ASYMMETRIC ) then
        fo(i) = coeff( l, 1, ibc(m) ) * ( SIGN(ONE, real(n - i - 1, wp) ) * fi(ip1) - &
                                                                      TWO * fi(i  ) + &
                                          SIGN(ONE, real(i - 1 - 1, wp) ) * fi(im1) ) + &
                coeff( l, 2, ibc(m) ) * ( SIGN(ONE, real(n - i - 2, wp) ) * fi(ip2) - &
                                                                      TWO * fi(i  ) + &
                                          SIGN(ONE, real(i - 2 - 1, wp) ) * fi(im2) )
      else! all other b.c. = interpolation = no specified bc = neumann
        if(i == 1 .or. i == n) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * fi(i)                 + &
                  coeff( l, 2, IBC_INTRPL ) * fi(i + SIGN(1, 1-m))  + &
                  coeff( l, 3, IBC_INTRPL ) * fi(i + SIGN(2, 1-m))  + &
                  coeff( l, 4, IBC_INTRPL ) * fi(i + SIGN(3, 1-m)) 
        else if (i == 2 .or. i == n - 1) then
          fo(i) = coeff( l, 1, IBC_INTRPL ) * ( fi(ip1) - TWO * fi(i) + fi(im1) )
        else
        end if
      end if
      
    end do
!-------------------------------------------------------------------------------
!   bulk area
!-------------------------------------------------------------------------------
    l = 3
    do i = 3, n - 2
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2
      fo(i) = coeff( l, 1, ibc(m) ) * ( fi(ip1) - TWO * fi(i) + fi(im1) ) + &
              coeff( l, 2, ibc(m) ) * ( fi(ip2) - TWO * fi(i) + fi(im2) )
    end do
!-------------------------------------------------------------------------------
!   mesh-based scaling
!-------------------------------------------------------------------------------
    fo(:) = fo(:) * dd ! dd = (1/dx)^2

    return
  end subroutine Prepare_TDMA_2deri_P2P_RHS_array
!===============================================================================
!> \brief To caculate the mid-point interpolation in 1D.
!> This subroutine is called as required to get the mid-point interpolation.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!_______________________________________________________________________________
  subroutine Get_x_midp_C2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)

    integer :: ixsub, nsz

    ixsub = dm%idom
    
    nsz = size(fo)
    fo = ZERO
    
    call Prepare_TDMA_interp_C2P_RHS_array(fi(:), fo(:), nsz, dm%icnbr_c2p(:,:), m1rC2P(:, :, :), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%am1x_C2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bm1x_C2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cm1x_C2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dm1x_C2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_midp_C2P_1D
!===============================================================================
  subroutine Get_x_midp_P2C_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: ixsub, nsz

    nsz = size(fo)
    fo = ZERO
    ixsub = dm%idom

    call Prepare_TDMA_interp_P2C_RHS_array(fi(:), fo(:), nsz, dm%ipnbr(:,:), m1rP2C(:, :, :), ibc(:))

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%am1x_P2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bm1x_P2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cm1x_P2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dm1x_P2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_midp_P2C_1D
!===============================================================================
  subroutine Get_y_midp_C2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use tridiagonal_matrix_algorithm
    use udf_type_mod
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_interp_C2P_RHS_array(fi(:), fo(:), nsz, dm%jcnbr_c2p(:,:), m1rC2P(:, :, :), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          am1y_C2P(:, ibc(1), ibc(2)), &
          bm1y_C2P(:, ibc(1), ibc(2)), &
          cm1y_C2P(:, ibc(1), ibc(2)), &
          dm1y_C2P(:, ibc(1), ibc(2)), &
          nsz)
! stretching? No stretching conversion
    return
  end subroutine Get_y_midp_C2P_1D
!===============================================================================
  subroutine Get_y_midp_P2C_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: nsz
    
    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_interp_P2C_RHS_array(fi(:), fo(:), nsz, dm%jpnbr(:, :), m1rP2C(:, :, :), ibc(:) )

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          am1y_P2C(:, ibc(1), ibc(2)), &
          bm1y_P2C(:, ibc(1), ibc(2)), &
          cm1y_P2C(:, ibc(1), ibc(2)), &
          dm1y_P2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_y_midp_P2C_1D
!===============================================================================
  subroutine Get_z_midp_C2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO
    
    call Prepare_TDMA_interp_C2P_RHS_array(fi(:), fo(:), nsz, dm%kcnbr_c2p(:,:), m1rC2P(:, :, :), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          am1z_C2P(:, ibc(1), ibc(2)), &
          bm1z_C2P(:, ibc(1), ibc(2)), &
          cm1z_C2P(:, ibc(1), ibc(2)), &
          dm1z_C2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_midp_C2P_1D
!===============================================================================
  subroutine Get_z_midp_P2C_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_interp_P2C_RHS_array(fi(:), fo(:), nsz, dm%kpnbr(:, :), m1rP2C(:, :, :), ibc(:))

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          am1z_P2C(:, ibc(1), ibc(2)), &
          bm1z_P2C(:, ibc(1), ibc(2)), &
          cm1z_P2C(:, ibc(1), ibc(2)), &
          dm1z_P2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_midp_P2C_1D
!===============================================================================
!> \brief To caculate the 1st derivative in 1D.
!> This subroutine is called as required to get the 1st derivative
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!===============================================================================
  subroutine Get_x_1st_derivative_C2C_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP),           intent(in   ) :: fbc(2)
    integer :: ixsub, nsz

    ixsub = dm%idom
    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_C2C_RHS_array(fi(:), fo(:), nsz, dm%icnbr_c2c(:, :), d1rC2C(:, :, :), dm%h1r(1), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%ad1x_C2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bd1x_C2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cd1x_C2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dd1x_C2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_1st_derivative_C2C_1D
!===============================================================================
  subroutine Get_x_1st_derivative_P2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)

    integer :: ixsub, nsz

    nsz = size(fo)
    fo = ZERO
    ixsub = dm%idom

    call Prepare_TDMA_1deri_P2P_RHS_array(fi(:), fo(:), nsz, dm%ipnbr_p2p(:, :), d1rP2P(:, :, :), dm%h1r(1), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%ad1x_P2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bd1x_P2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cd1x_P2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dd1x_P2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_1st_derivative_P2P_1D
!===============================================================================
  subroutine Get_x_1st_derivative_C2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)

    integer :: ixsub, nsz

    nsz = size(fo)
    fo = ZERO

    ixsub = dm%idom

    call Prepare_TDMA_1deri_C2P_RHS_array(fi(:), fo(:), nsz, dm%icnbr_c2p(:, :), d1rC2P(:, :, :), dm%h1r(1), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%ad1x_C2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bd1x_C2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cd1x_C2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dd1x_C2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_1st_derivative_C2P_1D
!===============================================================================
  subroutine Get_x_1st_derivative_P2C_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: nsz, ixsub

    nsz = size(fo)
    fo = ZERO

    ixsub = dm%idom

    call Prepare_TDMA_1deri_P2C_RHS_array(fi(:), fo(:), nsz, dm%ipnbr_p2c(:, :), d1rP2C(:, :, :), dm%h1r(1), ibc(:) )

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%ad1x_P2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bd1x_P2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cd1x_P2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dd1x_P2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_1st_derivative_P2C_1D
!===============================================================================
! y - Get_1st_derivative_1D
!===============================================================================
  subroutine Get_y_1st_derivative_C2C_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP),           intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_C2C_RHS_array(fi(:), fo(:), nsz, dm%jcnbr_c2c(:, :), d1rC2C(:, :, :), dm%h1r(2), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          ad1y_C2C(:, ibc(1), ibc(2)), &
          bd1y_C2C(:, ibc(1), ibc(2)), &
          cd1y_C2C(:, ibc(1), ibc(2)), &
          dd1y_C2C(:, ibc(1), ibc(2)), &
          nsz)

    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingcc(:, 1)

    return
  end subroutine Get_y_1st_derivative_C2C_1D
!===============================================================================
  subroutine Get_y_1st_derivative_P2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_P2P_RHS_array(fi(:), fo(:), nsz, dm%jpnbr_p2p(:, :), d1rP2P(:, :, :), dm%h1r(2), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          ad1y_P2P(:, ibc(1), ibc(2)), &
          bd1y_P2P(:, ibc(1), ibc(2)), &
          cd1y_P2P(:, ibc(1), ibc(2)), &
          dd1y_P2P(:, ibc(1), ibc(2)), &
          nsz)
    
    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingpt(:, 1)

    return
  end subroutine Get_y_1st_derivative_P2P_1D
!===============================================================================
  subroutine Get_y_1st_derivative_C2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_C2P_RHS_array(fi(:), fo(:), nsz, dm%jcnbr_c2p(:, :), d1rC2P(:, :, :), dm%h1r(2), ibc(:), fbc(:) )

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          ad1y_C2P(:, ibc(1), ibc(2)), &
          bd1y_C2P(:, ibc(1), ibc(2)), &
          cd1y_C2P(:, ibc(1), ibc(2)), &
          dd1y_C2P(:, ibc(1), ibc(2)), &
          nsz)
  
    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingpt(:, 1)

    return
  end subroutine Get_y_1st_derivative_C2P_1D
!===============================================================================
  subroutine Get_y_1st_derivative_P2C_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_P2C_RHS_array(fi(:), fo(:), nsz, dm%jpnbr_p2c(:, :), d1rP2C(:, :, :), dm%h1r(2), ibc(:))

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          ad1y_P2C(:, ibc(1), ibc(2)), &
          bd1y_P2C(:, ibc(1), ibc(2)), &
          cd1y_P2C(:, ibc(1), ibc(2)), &
          dd1y_P2C(:, ibc(1), ibc(2)), &
          nsz)

    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingcc(:, 1)

    return
  end subroutine Get_y_1st_derivative_P2C_1D
!===============================================================================
! z - Get_1st_derivative_1D
!===============================================================================
  subroutine Get_z_1st_derivative_C2C_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP),           intent(in   ) :: fbc(2)

    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_C2C_RHS_array(fi(:), fo(:), nsz, dm%kcnbr_c2c(:,:), d1rC2C(:, :, :), dm%h1r(3), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          ad1z_C2C(:, ibc(1), ibc(2)), &
          bd1z_C2C(:, ibc(1), ibc(2)), &
          cd1z_C2C(:, ibc(1), ibc(2)), &
          dd1z_C2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_1st_derivative_C2C_1D
!===============================================================================
  subroutine Get_z_1st_derivative_P2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_P2P_RHS_array(fi(:), fo(:), nsz, dm%kpnbr_p2p(:,:), d1rP2P(:, :, :), dm%h1r(3), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          ad1z_P2P(:, ibc(1), ibc(2)), &
          bd1z_P2P(:, ibc(1), ibc(2)), &
          cd1z_P2P(:, ibc(1), ibc(2)), &
          dd1z_P2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_1st_derivative_P2P_1D
!===============================================================================
  subroutine Get_z_1st_derivative_C2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_C2P_RHS_array(fi(:), fo(:), nsz, dm%kcnbr_c2p(:,:), d1rC2P(:, :, :), dm%h1r(3), ibc(:), fbc(:) )

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          ad1z_C2P(:, ibc(1), ibc(2)), &
          bd1z_C2P(:, ibc(1), ibc(2)), &
          cd1z_C2P(:, ibc(1), ibc(2)), &
          dd1z_C2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_1st_derivative_C2P_1D
!===============================================================================
  subroutine Get_z_1st_derivative_P2C_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_1deri_P2C_RHS_array(fi(:), fo(:), nsz, dm%kpnbr_p2c(:,:), d1rP2C(:, :, :), dm%h1r(3), ibc(:))

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          ad1z_P2C(:, ibc(1), ibc(2)), &
          bd1z_P2C(:, ibc(1), ibc(2)), &
          cd1z_P2C(:, ibc(1), ibc(2)), &
          dd1z_P2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_1st_derivative_P2C_1D
!===============================================================================
!> \brief To caculate the 2nd derivative in 1D.
!> This subroutine is called as required to get the 2nd derivative
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!===============================================================================
  subroutine Get_x_2nd_derivative_C2C_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP),           intent(in   ) :: fbc(2)
    integer :: nsz, ixsub

    nsz = size(fo)
    fo = ZERO
    ixsub = dm%idom

    call Prepare_TDMA_2deri_C2C_RHS_array(fi(:), fo(:), nsz, dm%icnbr_c2c(:,:), d2rC2C(:, :, :), dm%h2r(1), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%ad2x_C2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bd2x_C2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cd2x_C2C(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dd2x_C2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_2nd_derivative_C2C_1D
!===============================================================================
  subroutine Get_x_2nd_derivative_P2P_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: nsz, ixsub

    nsz = size(fo)
    fo = ZERO
    ixsub = dm%idom

    call Prepare_TDMA_2deri_P2P_RHS_array(fi(:), fo(:), nsz, dm%ipnbr_p2p(:,:), d2rP2P(:, :, :), dm%h2r(1), ibc(:))

    call Solve_TDMA(dm%is_periodic(1), fo(:), &
          xtdma_lhs(ixsub)%ad2x_P2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%bd2x_P2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%cd2x_P2P(:, ibc(1), ibc(2)), &
          xtdma_lhs(ixsub)%dd2x_P2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_x_2nd_derivative_P2P_1D
!===============================================================================
! y - Get_2nd_derivative_1D
!===============================================================================
  subroutine Get_y_2nd_derivative_C2C_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP),           intent(in   ) :: fbc(2)
    integer :: nsz
    real(WP), allocatable :: fo1(:)

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_2deri_C2C_RHS_array(fi(:), fo(:), nsz, dm%jcnbr_c2c(:,:), d2rC2C(:, :, :), dm%h2r(2), ibc(:), fbc(:) )

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          ad2y_C2C(:, ibc(1), ibc(2)), &
          bd2y_C2C(:, ibc(1), ibc(2)), &
          cd2y_C2C(:, ibc(1), ibc(2)), &
          dd2y_C2C(:, ibc(1), ibc(2)), &
          nsz)

    if(dm%is_stretching(2)) then 
      allocate ( fo1(nsz) ); fo1(:) = ZERO
      call Prepare_TDMA_1deri_C2C_RHS_array(fi(:), fo1(:), nsz, dm%jcnbr_c2c(:,:), d1rC2C(:, :, :), dm%h1r(2), ibc(:), fbc(:))
      call Solve_TDMA(dm%is_periodic(2), fo1(:), &
           ad1y_C2C(:, ibc(1), ibc(2)), &
           bd1y_C2C(:, ibc(1), ibc(2)), &
           cd1y_C2C(:, ibc(1), ibc(2)), &
           dd1y_C2C(:, ibc(1), ibc(2)), &
           nsz)
      fo(:) = fo(:) * dm%yMappingcc(:, 2) + fo1(:) * dm%yMappingcc(:, 3)
      deallocate (fo1)
    end if

    return
  end subroutine Get_y_2nd_derivative_C2C_1D
!===============================================================================
  subroutine Get_y_2nd_derivative_P2P_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP), optional, intent(in   ) :: fbc(2)

    integer :: nsz
    real(WP), allocatable :: fo1(:)

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_2deri_P2P_RHS_array(fi(:), fo(:), nsz, dm%jpnbr_p2p(:,:), d2rP2P(:, :, :), dm%h2r(2), ibc(:))

    call Solve_TDMA(dm%is_periodic(2), fo(:), &
          ad2y_P2P(:, ibc(1), ibc(2)), &
          bd2y_P2P(:, ibc(1), ibc(2)), &
          cd2y_P2P(:, ibc(1), ibc(2)), &
          dd2y_P2P(:, ibc(1), ibc(2)), &
          nsz)

    if(dm%is_stretching(2)) then 
      allocate ( fo1(nsz) ); fo1(:) = ZERO

      call Prepare_TDMA_1deri_P2P_RHS_array(fi(:), fo1(:), nsz, dm%jpnbr_p2p(:,:), d1rP2P(:, :, :), dm%h1r(2), ibc(:), fbc(:))

      call Solve_TDMA(dm%is_periodic(2), fo1(:), &
           ad1y_P2P(:, ibc(1), ibc(2)), &
           bd1y_P2P(:, ibc(1), ibc(2)), &
           cd1y_P2P(:, ibc(1), ibc(2)), &
           dd1y_P2P(:, ibc(1), ibc(2)), &
           nsz)
      fo(:) = fo(:) * dm%yMappingpt(:, 2) + fo1(:) * dm%yMappingpt(:, 3)
      deallocate (fo1)
    end if

    return
  end subroutine Get_y_2nd_derivative_P2P_1D
!===============================================================================
! z - Get_2nd_derivative_1D
!===============================================================================
  subroutine Get_z_2nd_derivative_C2C_1D (fi, fo, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    real(WP),           intent(in   ) :: fbc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_2deri_C2C_RHS_array(fi(:), fo(:), nsz, dm%kcnbr_c2c(:,:), d2rC2C(:, :, :), dm%h2r(3), ibc(:), fbc(:))

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          ad2z_C2C(:, ibc(1), ibc(2)), &
          bd2z_C2C(:, ibc(1), ibc(2)), &
          cd2z_C2C(:, ibc(1), ibc(2)), &
          dd2z_C2C(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_2nd_derivative_C2C_1D
!===============================================================================
  subroutine Get_z_2nd_derivative_P2P_1D (fi, fo, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in   ) :: fi(:)
    real(WP),           intent(inout) :: fo(:)
    type(t_domain),     intent(in   ) :: dm
    integer,            intent(in   ) :: ibc(2)
    integer :: nsz

    nsz = size(fo)
    fo = ZERO

    call Prepare_TDMA_2deri_P2P_RHS_array(fi(:), fo(:), nsz, dm%kpnbr_p2p(:,:), d2rP2P(:, :, :), dm%h2r(3), ibc(:))

    call Solve_TDMA(dm%is_periodic(3), fo(:), &
          ad2z_P2P(:, ibc(1), ibc(2)), &
          bd2z_P2P(:, ibc(1), ibc(2)), &
          cd2z_P2P(:, ibc(1), ibc(2)), &
          dd2z_P2P(:, ibc(1), ibc(2)), &
          nsz)

    return
  end subroutine Get_z_2nd_derivative_P2P_1D
!===============================================================================
!> \brief To caculate the mid-point interpolation in 3D.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!===============================================================================
  subroutine Get_x_midp_C2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use tridiagonal_matrix_algorithm
    use udf_type_mod
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  default : x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(ibc(1) == IBC_DIRICHLET .or. ibc(2) == IBC_DIRICHLET) then
          call Get_x_midp_C2P_1D (fi, fo, dm, ibc, fbc)
        else
          call Get_x_midp_C2P_1D (fi, fo, dm, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_midp_C2P_3D
!===============================================================================
  subroutine Get_x_midp_P2C_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  default : x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Get_x_midp_P2C_1D (fi, fo, dm, ibc)
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_midp_P2C_3D
!===============================================================================
  subroutine Get_y_midp_C2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
!-------------------------------------------------------------------------------
!  default : y-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(ibc(1) == IBC_DIRICHLET .or. ibc(2) == IBC_DIRICHLET) then
          call Get_y_midp_C2P_1D (fi, fo, dm, ibc, fbc)
        else
          call Get_y_midp_C2P_1D (fi, fo, dm, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_midp_C2P_3D
!===============================================================================
  subroutine Get_y_midp_P2C_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
!-------------------------------------------------------------------------------
!  y-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Get_y_midp_P2C_1D (fi, fo, dm, ibc)
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_midp_P2C_3D
  !===============================================================================
  subroutine Get_z_midp_C2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  default : z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(ibc(1) == IBC_DIRICHLET .or. ibc(2) == IBC_DIRICHLET) then
          call Get_z_midp_C2P_1D (fi, fo, dm, ibc, fbc)
        else
          call Get_z_midp_C2P_1D (fi, fo, dm, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_midp_C2P_3D
!===============================================================================
  subroutine Get_z_midp_P2C_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  default : z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Get_z_midp_P2C_1D (fi, fo, dm, ibc)
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_midp_P2C_3D
!===============================================================================
!> \brief To caculate the 1st-deriviate in 3D.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!===============================================================================
  subroutine Get_x_1st_derivative_C2C_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP),           intent(in) :: fbc(2)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Get_x_1st_derivative_C2C_1D(fi, fo, dm, ibc, fbc)
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_C2C_3D
!===============================================================================
  subroutine Get_x_1st_derivative_P2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(ibc(1) == IBC_NEUMANN .or. ibc(2) == IBC_NEUMANN) then
          call Get_x_1st_derivative_P2P_1D(fi, fo, dm, ibc, fbc)
        else 
          call Get_x_1st_derivative_P2P_1D(fi, fo, dm, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_P2P_3D
!===============================================================================
  subroutine Get_x_1st_derivative_C2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(ibc(1) == IBC_NEUMANN .or. ibc(2) == IBC_NEUMANN) then
          call Get_x_1st_derivative_C2P_1D(fi, fo, dm, ibc, fbc)
        else 
          call Get_x_1st_derivative_C2P_1D(fi, fo, dm, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_C2P_3D
!===============================================================================
  subroutine Get_x_1st_derivative_P2C_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Get_x_1st_derivative_P2C_1D(fi, fo, dm, ibc)
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_P2C_3D
!===============================================================================
  subroutine Get_y_1st_derivative_C2C_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP),           intent(in) :: fbc(2)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i

    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Get_y_1st_derivative_C2C_1D(fi, fo, dm, ibc, fbc)
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_C2C_3D
!===============================================================================
  subroutine Get_y_1st_derivative_P2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
!!-------------------------------------------------------------------------------
!  y-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Get_y_1st_derivative_P2P_1D(fi, fo, dm, ibc, fbc)
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_P2P_3D
!===============================================================================
  subroutine Get_y_1st_derivative_C2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
!-------------------------------------------------------------------------------
!  y-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(ibc(1) == IBC_NEUMANN .or. ibc(2) == IBC_NEUMANN) then
          call Get_y_1st_derivative_C2P_1D(fi, fo, dm, ibc, fbc)
        else 
          call Get_y_1st_derivative_C2P_1D(fi, fo, dm, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_C2P_3D

!===============================================================================
  subroutine Get_y_1st_derivative_P2C_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
!-------------------------------------------------------------------------------
!  y-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Get_y_1st_derivative_P2C_1D(fi, fo, dm, ibc)
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_P2C_3D
!===============================================================================
  subroutine Get_z_1st_derivative_C2C_3D (fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP),           intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Get_z_1st_derivative_C2C_1D(fi, fo, dm, ibc, fbc)
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_C2C_3D

!===============================================================================
  subroutine Get_z_1st_derivative_P2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(ibc(1) == IBC_NEUMANN .or. ibc(2) == IBC_NEUMANN) then
          call Get_z_1st_derivative_P2P_1D(fi, fo, dm, ibc, fbc)
        else 
          call Get_z_1st_derivative_P2P_1D(fi, fo, dm, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_P2P_3D
!===============================================================================
  subroutine Get_z_1st_derivative_C2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(ibc(1) == IBC_NEUMANN .or. ibc(2) == IBC_NEUMANN) then
          call Get_z_1st_derivative_C2P_1D(fi, fo, dm, ibc, fbc)
        else 
          call Get_z_1st_derivative_C2P_1D(fi, fo, dm, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_C2P_3D
  !===============================================================================
  subroutine Get_z_1st_derivative_P2C_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Get_z_1st_derivative_P2C_1D(fi, fo, dm, ibc)
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_P2C_3D
!===============================================================================
!> \brief To caculate the 2nd-deriviate in 3D.
!------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!===============================================================================
  subroutine Get_x_2nd_derivative_C2C_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP),           intent(in) :: fbc(2)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Get_x_2nd_derivative_C2C_1D(fi, fo, dm, ibc, fbc)
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_x_2nd_derivative_C2C_3D
  !===============================================================================
  subroutine Get_x_2nd_derivative_P2P_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
!-------------------------------------------------------------------------------
!  x-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Get_x_2nd_derivative_P2P_1D(fi, fo, dm, ibc)
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_x_2nd_derivative_P2P_3D
!===============================================================================
  subroutine Get_y_2nd_derivative_C2C_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP),           intent(in) :: fbc(2)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
!-------------------------------------------------------------------------------
!  y-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Get_y_2nd_derivative_C2C_1D(fi, fo, dm, ibc, fbc)
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_y_2nd_derivative_C2C_3D
!===============================================================================
  subroutine Get_y_2nd_derivative_P2P_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc(2)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
!-------------------------------------------------------------------------------
!  y-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(ibc(1) == IBC_NEUMANN .or. ibc(2) == IBC_NEUMANN) then
          call Get_y_2nd_derivative_P2P_1D(fi, fo, dm, ibc, fbc)
        else 
          call Get_y_2nd_derivative_P2P_1D(fi, fo, dm, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_y_2nd_derivative_P2P_3D
!===============================================================================
  subroutine Get_z_2nd_derivative_C2C_3D(fi3d, fo3d, dm, ibc, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)
    real(WP),           intent(in) :: fbc(2)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Get_z_2nd_derivative_C2C_1D(fi, fo, dm, ibc, fbc)
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return
  end subroutine Get_z_2nd_derivative_C2C_3D
!===============================================================================
  subroutine Get_z_2nd_derivative_P2P_3D(fi3d, fo3d, dm, ibc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: ibc(2)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
!-------------------------------------------------------------------------------
!  z-pencil calculation
!-------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Get_z_2nd_derivative_P2P_1D(fi, fo, dm, ibc)
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return
  end subroutine Get_z_2nd_derivative_P2P_3D

end module