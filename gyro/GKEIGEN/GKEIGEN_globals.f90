!---------------------------------------
! GKEIGEN_globals.f90
!
! PURPOSE:
!  Module containing most shared variables
!  for GKEIGEN routines.
!---------------------------------------

module GKEIGEN_globals

  !------------------------------------
  ! Character variables
  !
  ! in GKEIGEN_matrix_read:
  !
  character(len=14), parameter :: &
       file_eigen_restart = &
       "GKEIGEN_matrix"
  !
  ! in GKEIGEN_do:
  !
  character(len=16), parameter :: &
       file_eigen_freq = &
       "GKEIGEN_freq.out"
  !----------------------------------

  !------------------------------------------------------------
  ! Matrix and vector storage variables
  !
  complex, dimension(:,:,:,:,:), allocatable :: eigensolve_vec
  !
  complex, dimension(:), allocatable :: eigensolve_value
  complex, dimension(:,:), allocatable :: eigensolve_matrix
!  complex, dimension(:,:), allocatable :: M_loc
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
  real :: ev_omega
  real :: ev_gamma
  !
  ! Total state vector index
  integer :: istate          !-------------- 
  integer :: istate_1        !
  integer :: istate_start    ! column index 
  integer :: istate_end      !
  integer :: ielem           !--------------
  integer :: jstate          !--------------
  integer :: jstate_start    !
  integer :: jstate_end      ! row indices
  integer :: jset_jstart     !
  integer :: jset_end        !
  integer :: jelem           !
  integer :: jelem0          !
  integer :: i_index         !
  integer :: jelem_start     !
  integer :: jelem_end       !--------------
  integer :: iseq         ! Merged row and
  !                    ---! column index   
  !
  ! Processors in GYRO_COMM_WORLD
  integer :: i_proc_e
  !
  ! Parallel instances of GYRO_COMM_WORLD
  integer :: gkeigen_j_set_e0
  integer :: gkeigen_j_set_e1
  integer :: gkeigen_j_set_e2
  integer :: j_sr
  integer :: i_sr
  !
  ! Rank of processor over MPI_COMM_WORLD
  integer :: j_proc_tot_e1
  integer :: j_proc_tot_0
  integer :: j_send_proc
  integer :: j_recv_proc
  !
  ! n_nek = n_n*n_energy*n_lambda
  !
  integer :: p_nek1
  integer :: p_nek2
  !
  ! Locally owned n_nek points
  integer :: p_nek_loc1
  integer :: p_nek_loc2
  !
  !---------------------------------------------------------

  ! Logicals
  logical :: l_in_i_block
  logical :: l_in_j_block
  logical :: l_print
  logical :: l_set       ! Controls MPI of matrix
  logical :: l_send      ! blocks of PETSc object
  logical :: l_recv      ! time_derivative_matrix.
  !
  ! Restart control
!  integer :: irestart
!  complex, dimension(500,h_length_loc) :: restart_cols
  !

  integer, dimension(8) :: time_array

end module GKEIGEN_globals
