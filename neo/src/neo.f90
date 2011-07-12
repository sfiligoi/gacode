program neo

  use neo_globals, only : path
  implicit none

  integer :: ierr

  include 'mpif.h'

  call MPI_INIT(ierr)
  
  ! Path is cwd:
  path= './'


  call neo_read_input
  call neo_do

  call MPI_FINALIZE(ierr)

end program neo
