module parsolvers
implicit none
!
    real*8, parameter :: tol_abs = 1.0d-14   ! Tolerance of residual.
    real*8, parameter :: tol_rel = 1.0d-15   !Tolerance realitive initial rho
! bug in the maximum iteration, fix it
    integer*4, parameter :: maxits = 400     ! maximum iteration permited
    integer*4, parameter :: mr = 100         ! Inner iterations on MGMRES
    integer*4, parameter :: lfil = 4         ! fill-in elem, good:(5<lfil<10)
    real*8, parameter :: tol = 0.005d0       ! tolerance of fill-in, good: 0.005
    integer*4, parameter :: iout = 1         ! if iout NOT zero print residual file
!
!
end module parsolvers
!******************************************************************************
module parmesh
implicit none
!
! Number of outer iterations on main program
! Steady state set this = 0
    integer*4, parameter :: itr_max = 2000
! Norm of temperature between two time steps
    real*8, parameter :: epsolon = 10.0d-6
!
! msh are options: mesh100, mesh500, mesh1500, mesh2500, 
! mesh5000, mesh10000, mesh15000, mesh40000
  character(*), parameter :: extension = '.dat'
  character(*), parameter :: mshname = 'mesh500'
! concatenating the mesh name and the extension
  character(*), parameter :: meshfile = mshname//extension
! Chosing the iterative method and the preconditioner
! 1 == Gauss Siedel
! 2 == GMRES w/out preconditioner
! 3 == GMRES w/ ILU0
! 4 == GMRES w/ ILUT
! 5 == GMRES w/ MILU0
! 6 == GMRES w/ ILU(K)
    integer(kind=4), parameter :: choice = 6
!
! Chosing the type of analysis
! 1 == unsteady state
! 2 == steady state second order discretization
! 3 == steady state fourth order discretization
!
! REMEMBER IF IS IN STATE STATE SET ITR_MAX = 0
    integer(kind=4), parameter :: choice2 = 1
!
!
! This parameter is for subroutine log and sol_exact
! if the problem has a analytical solution this parameter
! should be set to 1, if don't set to 0
    integer(kind=4), parameter :: choice3 = 1
!
end module parmesh
!******************************************************************************
module chars
implicit none
!
!
! Defining names and directories of the solution files
    character(*), parameter :: filename = 'sol_aproximate.plt'
    character(*), parameter :: filename_exact = 'sol_exact.plt'
    character(*), parameter :: savename = 'save_sol.dat'
    character(*), parameter :: mshdir = '/home/vedat/HTS2D/mesh/'
!
end module chars
!******************************************************************************
module physics
!
! This module contains all physics constants of the problem
!
   real*8 :: h = 1.0d0                      ! W/cm²K
   real*8 :: e = 1.0d0                      ! W/cm³
   real*8 :: alpha = 1.0d0                  ! cm²/s
   real*8 :: temp_inf = 273.15d0            ! K
   real*8 :: k = 1.0d0                      ! W/cmK
   real*8 :: deltat = 0.00518d0             ! [SECS]
!
end module physics
