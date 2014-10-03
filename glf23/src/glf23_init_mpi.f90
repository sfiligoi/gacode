!\
! module for MPI version of TGLF
!/

      MODULE glf23_mpi

        use mpi

        IMPLICIT NONE

        ! local communicator
        integer :: iCommGlf    = -1
        integer :: iProcGlf    = -1
        integer :: nProcGlf    = -1
        integer :: iProc0Glf   = 0
        integer :: iGroupIDGlf = 1

      END MODULE glf23_mpi
!
!______________________________________________________
!
!--------------------------------------------------------------
! glf23_init_mpi.f90
!
! PURPOSE:
!  Initialize external TGLF interface.
!---------------------------------------------------------------

subroutine glf23_init(path_in,mpi_comm_in)

  use glf23_interface
  use glf23_mpi

  implicit none

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer,          intent(in) :: mpi_comm_in
  integer :: ierr

  ! Local variable
  logical :: inputdat_flag

  glf23_path_in = path_in
  iCommGlf    = mpi_comm_in

  call MPI_COMM_RANK(iCommGlf,iProcGlf,ierr)
  call MPI_COMM_SIZE(iCommGlf,nProcGlf,ierr)

  inquire(file=trim(path_in)//'input.glf23.gen',exist=inputdat_flag)

  if (inputdat_flag .eqv. .true.) then
     if (glf23_quiet_flag_in .eqv. .false.) then
        print '(a,a,a)', 'INFO: (GLF23) glf23_init reading ',trim(path_in),'input.glf23.gen'
     endif
     call glf23_read_input
  endif

end subroutine glf23_init

