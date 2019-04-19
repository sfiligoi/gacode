subroutine vpro_icomm(c,c0,p)

  use mpi

  implicit none

  integer, intent(inout) :: p
  integer :: ierr,iproc
  logical :: flag
  character*22 :: c,c0

  if (c /= c0) return
  
  call MPI_INITIALIZED(flag,ierr)

  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,*) p
     call MPI_BCAST(p,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  else
     read(1,*) p
  endif

end subroutine vpro_icomm

