subroutine EXPRO_palloc(comm,path,flag)

  use EXPRO_interface

  implicit none

  integer, intent(in) :: comm
  character(len=*), intent(in) :: path
  integer, intent(in) :: flag
  integer :: i_proc
  integer :: ierr

  include 'mpif.h'

  call MPI_COMM_RANK(comm,i_proc,ierr)
  
  if (i_proc == 0) then
     ! Read dimensions and do allocations on 0
     call EXPRO_alloc_control(i_proc,path,flag)
  endif

  if (flag == 1) then

     ! Broadcast dimensions

     call MPI_BCAST(EXPRO_n_exp,&
          1,&
          MPI_INTEGER,&
          0,&
          comm,&
          ierr)

     call MPI_BCAST(EXPRO_nfourier,&
          1,&
          MPI_INTEGER,&
          0,&
          comm,&
          ierr)

  endif

  if (i_proc /= 0) then
     ! Do allocations on remaining processes
     call EXPRO_alloc_control(i_proc,path,flag)
  endif

end subroutine EXPRO_palloc
