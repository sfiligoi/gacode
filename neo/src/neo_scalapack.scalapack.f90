subroutine neo_scalapack(info)

  use neo_globals
  use mpi

  implicit none

  integer, intent(inout) :: info

  real, dimension(:), allocatable :: work
  integer :: lwork

  integer :: desca(7)
  integer :: descb(7)

  call SL_INIT(ictxt,1,n_proc)  

  lwork = 2*(npb+bw)*bw+36*bw**2+npb+6*bw
  allocate(work(lwork))

  ! Matrix (block-column distribution)
  desca(1) = 501
  desca(2) = ictxt
  desca(3) = n_row
  desca(4) = npb
  desca(5) = 0
  desca(6) = ldab
  desca(7) = 0

  ! Righthand-side vectors (block-row distribution)
  descb(1) = 502
  descb(2) = ictxt
  descb(3) = n_row
  descb(4) = npb
  descb(5) = 0
  descb(6) = npb
  descb(7) = 0

  call PDGBSV(n_row,bw,bw,1,ab,1,desca,ipiv,g_loc,1,descb,work,lwork,info)

  ! Construct global solution vector
  call MPI_ALLGATHER(g_loc,&
       npb,&
       MPI_DOUBLE_PRECISION,&
       g,&
       npb,&
       MPI_DOUBLE_PRECISION,&
       MPI_COMM_WORLD)

end subroutine neo_scalapack
