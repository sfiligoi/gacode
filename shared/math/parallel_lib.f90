
!  parallel_lib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
!  parallel_lib_r -> g(nj_loc,ni) -> f(ni_loc,nj)

module parallel_lib

  implicit none

  integer, private :: nproc,iproc
  integer, private :: ni,nj
  integer, private :: ni_loc,nj_loc
  integer, private :: transp_comm
  integer, private :: nsend

  complex, dimension(:,:,:), allocatable, private :: fsendf
  complex, dimension(:,:,:), allocatable, private :: fsendr

contains

  !=========================================================

  subroutine parallel_lib_init(ni_in,nj_in,ni_loc_out,nj_loc_out,comm)

    use mpi

    implicit none

    integer, intent(in) :: ni_in,nj_in
    integer, intent(in) :: comm
    integer, intent(inout) :: ni_loc_out,nj_loc_out
    integer, external :: parallel_dim
    integer :: ierr

    transp_comm = comm

    call MPI_COMM_RANK(TRANSP_COMM,iproc,ierr)
    call MPI_COMM_SIZE(TRANSP_COMM,nproc,ierr)

    ni = ni_in
    nj = nj_in

    ni_loc = parallel_dim(ni,nproc)
    nj_loc = parallel_dim(nj,nproc)

    ni_loc_out = ni_loc
    nj_loc_out = nj_loc

    nsend = ni*nj/nproc**2

    allocate(fsendf(nj_loc,ni_loc,nproc))
    allocate(fsendr(ni_loc,nj_loc,nproc))

  end subroutine parallel_lib_init

  !=========================================================

  subroutine parallel_lib_f(f,ft)

    use mpi

    implicit none

    complex, intent(in), dimension(ni_loc,nj) :: f
    complex, intent(inout), dimension(nj_loc,ni) :: ft
    integer :: ierr,i_loc,i,j,k


    do k=1,nproc
       i_loc = 0
       do i=1+iproc*ni_loc,(1+iproc)*ni_loc
          i_loc = i_loc+1
          do j=1,nj_loc
             fsendf(j,i_loc,k) = f(i_loc,j+(k-1)*nj_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         transp_comm, &
         ierr)

  end subroutine parallel_lib_f

  !=========================================================

  subroutine parallel_lib_r(ft,f)

    use mpi

    implicit none

    complex, intent(in), dimension(nj_loc,ni) :: ft
    complex, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k

    do k=1,nproc
       j_loc = 0
       do j=1+iproc*nj_loc,(1+iproc)*nj_loc
          j_loc = j_loc+1
          do i=1,ni_loc
             fsendr(i,j_loc,k) = ft(j_loc,i+(k-1)*ni_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         transp_comm, &
         ierr)

  end subroutine parallel_lib_r

  !=========================================================

end module parallel_lib
