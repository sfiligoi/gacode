!--------------------------------------------------------------
! tglf_init.f90
!
! PURPOSE:
!  Initialize external TGLF interface.
!---------------------------------------------------------------

subroutine tglf_init(path_in,mpi_comm_in)

  use tglf_interface

  implicit none

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer,          intent(in) :: mpi_comm_in

  ! Local variable
  logical :: inputdat_flag

  tglf_path_in = path_in

  inquire(file=trim(path_in)//'input.tglf.gen',exist=inputdat_flag)

  if (inputdat_flag .eqv. .true.) then
     print '(a,a,a)', '[tglf_init reading ',trim(path_in),'input.tglf.gen]'
     call tglf_read_input
  else
     !print '(a,a,a)', '[tglf_init NOT reading ',trim(path_in),'input.tglf.gen]'
  endif

end subroutine tglf_init

