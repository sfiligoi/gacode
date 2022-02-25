!\
! module for MPI version of gftm
!/

      MODULE gftm_mpi

        use mpi

        IMPLICIT NONE

        ! local communicator
        integer :: iCommgftm    = -1
        integer :: iProcgftm    = -1
        integer :: nProcgftm    = -1
        integer :: iProc0gftm   = 0
        integer :: iGroupIDgftm = 1

      END MODULE gftm_mpi
!
!______________________________________________________
!
!--------------------------------------------------------------
! gftm_init.f90
!
! PURPOSE:
!  Initialize external gftm interface.
!---------------------------------------------------------------

subroutine gftm_init(path_in,mpi_comm_in)

  use gftm_interface
  use gftm_mpi

  implicit none

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer,          intent(in) :: mpi_comm_in
  integer :: ierr

  ! Local variable
  logical :: inputdat_flag

  gftm_path_in = path_in
  iCommgftm    = mpi_comm_in

  call MPI_COMM_RANK(iCommgftm,iProcgftm,ierr)
  call MPI_COMM_SIZE(iCommgftm,nProcgftm,ierr)

  inquire(file=trim(path_in)//'input.gftm.gen',exist=inputdat_flag)

  if (inputdat_flag .eqv. .true.) then
     if (gftm_quiet_flag_in .eqv. .false.) then
        print '(a,a,a)', 'INFO: (gftm) gftm_init reading ',trim(path_in),'input.gftm.gen'
     endif
     call gftm_read_input
  endif

end subroutine gftm_init

