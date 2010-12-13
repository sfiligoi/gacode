program neo

  implicit none

  integer :: ierr

  include 'mpif.h'

  call MPI_INIT(ierr)

  call neo_read_input
  call neo_do

  call MPI_FINALIZE(ierr)

end program neo
