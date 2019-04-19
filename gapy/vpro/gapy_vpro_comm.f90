subroutine vpro_icomm(p)

  use mpi

  implicit none

  integer, intent(inout) :: p
  integer :: ierr,iproc
  logical :: flag
  
  call MPI_INITIALIZED(flag,ierr)

  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,*) p
     call MPI_BCAST(p,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  else
     read(1,*) p
  endif

end subroutine vpro_icomm

subroutine vpro_rcomm(x)

  use mpi

  implicit none

  double precision, intent(inout) :: x
  integer :: ierr,iproc
  logical :: flag
  
  call MPI_INITIALIZED(flag,ierr)
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,10) x
     call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     read(1,10) x
  endif

10 format(1pe14.7)

end subroutine vpro_rcomm

subroutine vpro_scomm(x,n)

  use mpi

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: ierr,iproc
  logical :: flag
  
  call MPI_INITIALIZED(flag,ierr)
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,10) x
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     read(1,10) x
  endif

10 format(10(1pe14.7))

end subroutine vpro_scomm

subroutine vpro_vcomm(x,n)

  use mpi

  implicit none

  integer :: idum,i
  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: ierr,iproc
  logical :: flag

  call MPI_INITIALIZED(flag,ierr)

  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) then
        do i=1,n
           read(1,10) idum,x(i)
        enddo
     endif
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     do i=1,n
        read(1,10) idum,x(i)
     enddo
  endif

10 format(i3,1x,1pe14.7)

end subroutine vpro_vcomm

subroutine vpro_acomm(x,m,n)

  use mpi

  implicit none

  integer :: idum,i
  integer, intent(in) :: m,n
  double precision, intent(inout), dimension(m,n) :: x
  integer :: ierr,iproc
  logical :: flag

  call MPI_INITIALIZED(flag,ierr)

  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) then
        do i=1,n
           read(1,10) idum,x(:,i)
        enddo
     endif
     call MPI_BCAST(x,m*n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     do i=1,n
        read(1,10) idum,x(:,i)
     enddo
  endif

10 format(i3,1x,10(1pe14.7,1x))

end subroutine vpro_acomm

