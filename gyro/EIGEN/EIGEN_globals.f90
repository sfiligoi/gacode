!---------------------------------------
! EIGEN_globals.f90
!
! PURPOSE:
!  Module containing most shared variables
!  for EIGEN routines.
!---------------------------------------

module EIGEN_globals

  !------------------------------------
  ! Character variables
  !
  ! in EIGEN_matrix_read:
  !
  character(len=12), parameter :: &
       file_eigen_restart = &
       "EIGEN_matrix"
  !
  ! in EIGEN_do:
  !
  character(len=14), parameter :: &
       file_eigen_freq = &
       "EIGEN_freq.out"
  !----------------------------------

  !------------------------------------------------------------
  ! Matrix and vector storage variables
  !
  complex, dimension(:,:,:,:,:), allocatable :: eigensolve_vec
  !
  complex, dimension(:), allocatable :: eigensolve_value
  complex, dimension(:,:), allocatable :: eigensolve_matrix
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Iteration variables
  !
  ! The eigensolver requires two indices for each phase space
  ! degree of freedom.
  !
  ! Radial, orbit, species indices
  integer :: i1
  integer :: i2
  integer :: m1
  integer :: m2
  integer :: is1
  integer :: is2
  !
  ! Eigenvalues
  integer :: iev
  !
  ! Total state vector index
  integer :: istate
  integer :: istate_1
  integer :: jstate
  !
  ! n_nek = n_n*n_energy*n_lambda
  !
  integer :: p_nek1
  integer :: p_nek2
  !
  integer :: p_nek_loc1
  integer :: p_nek_loc2
  !
  !
  ! Processors
  integer :: i_proc_e
  !---------------------------------------------------------

end module EIGEN_globals
