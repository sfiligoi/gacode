module parallel_lib

  implicit none

  ! lib
  logical :: use_alltoall_lib_f = .true.
  logical :: use_alltoall_lib_r = .true.
  logical :: use_alltoall_slib_f = .true.
  logical :: use_alltoall_slib_r = .true.

  public :: use_alltoall_lib_f 
  public :: use_alltoall_lib_r 
  public :: use_alltoall_slib_f 
  public :: use_alltoall_slib_r 

  integer, parameter :: idebug = 0
  integer, parameter :: TAG_UB = 32767
  complex, parameter :: czero = (0.0,0.0)

  integer, private :: nproc,iproc
  integer, private :: ni,nj
  integer, private :: ni_loc
  integer, private :: nj_loc
  integer, private :: lib_comm
  integer, private :: nsend

  complex, dimension(:,:,:), allocatable, private :: fsendf
  complex, dimension(:,:,:), allocatable, private :: fsendr
  real, dimension(:,:,:), allocatable, private :: fsendr_real

  ! slib

  integer, private :: nsproc,isproc
  integer, private :: nn
  integer, private :: nexch
  integer, private :: nkeep
  integer, private :: nsplit
  integer, private :: slib_comm

  integer, private, parameter :: default_size = 1024*1024*32
  private :: assert

contains

  subroutine assert(lcond,msg,ival)
    implicit none
    logical,intent(in) :: lcond
    character*(*), intent(in) :: msg
    integer, intent(in) :: ival

    if (.not.lcond) then
       write(*,*) msg,ival
       stop 'assertion failed'
    endif

  end subroutine assert


  !=========================================================

  !  parallel_lib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_lib_r -> g(nj_loc,ni) -> f(ni_loc,nj)

  subroutine parallel_lib_init(ni_in,nj_in,ni_loc_out,nj_loc_out,comm)

    use mpi

    implicit none

    integer, intent(in) :: ni_in,nj_in
    integer, intent(in) :: comm
    integer, intent(inout) :: ni_loc_out,nj_loc_out
    integer, external :: parallel_dim
    integer :: ierr

    lib_comm = comm

    call MPI_COMM_RANK(lib_comm,iproc,ierr)
    call MPI_COMM_SIZE(lib_comm,nproc,ierr)

    ni = ni_in
    nj = nj_in

    ni_loc = parallel_dim(ni,nproc)
    nj_loc = parallel_dim(nj,nproc)

    ni_loc_out = ni_loc
    nj_loc_out = nj_loc

    nsend = ni*nj/nproc**2

    allocate(fsendf(nj_loc,ni_loc,nproc))
    allocate(fsendr(ni_loc,nj_loc,nproc))
    allocate(fsendr_real(ni_loc,nj_loc,nproc))

  end subroutine parallel_lib_init

  !=========================================================

  subroutine parallel_lib_f(f,ft)

    use mpi

    implicit none

    complex, intent(in), dimension(ni_loc,nj) :: f
    complex, intent(inout), dimension(nj_loc,ni) :: ft
    integer :: ierr,i_loc,i,j,k, i1,i2

    integer :: irequest,source, dest, tag
    integer :: request_list(nproc+nproc)

    if (idebug >= 1) then
       print 9010, iproc,nsend,ni_loc,nj,nj_loc,ni
9010   format('parallel_lib_f: iproc=',i5,' nsend=',i5,                   &
            ' ni_loc=',i5,' nj=',i5,' nj_loc=',i5,' ni=',i5)
    endif

    if (use_alltoall_lib_f) then

       do k=1,nproc
          i1 = 1+iproc*ni_loc
          i2 = (1+iproc)*ni_loc
          do i=i1,i2
             i_loc = i-i1+1 
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
            lib_comm, &
            ierr)
    else

       !     ----------------
       !     post mpi_irecv()
       !     ----------------

       do k=1,nproc
          source = (k-1)
          dest = iproc
          tag = mod(source + dest*nproc,TAG_UB)
          call MPI_Irecv( ft(1,k), nsend, MPI_DOUBLE_COMPLEX,         &
               source, tag, lib_comm, irequest,ierr)
          call assert(ierr.eq.MPI_SUCCESS,                            &
               'parallel_lib_f: MPI_Irecv return ierr',ierr) 
          request_list(k) = irequest
       enddo

       !     ----------------
       !     post mpi_isend()
       !     ----------------

       do k=1,nproc
          i1 = 1+iproc*ni_loc
          i2 = (1+iproc)*ni_loc
          do i=i1,i2
             i_loc = i-i1+1 
             do j=1,nj_loc
                fsendf(j,i_loc,k) = f(i_loc,j+(k-1)*nj_loc) 
             enddo
          enddo

          source = iproc
          dest = (k-1)
          tag = mod(source + dest*nproc,TAG_UB)
          call MPI_Isend(fsendf(1,1,k),nsend,MPI_DOUBLE_COMPLEX,            &
               dest, tag, lib_comm, irequest, ierr)
          call assert(ierr.eq.MPI_SUCCESS,                                  &
               'parallel_lib_f:MPI_Isend return ierr',ierr)
          request_list(nproc+k) = irequest

       enddo

       call MPI_Waitall(2*nproc,request_list, MPI_STATUSES_IGNORE,ierr)
       call assert(ierr.eq.MPI_SUCCESS,                                   &
            'parallel_lib_f:MPI_Waitall return ierr',ierr)

    endif

  end subroutine parallel_lib_f

  !=========================================================

  subroutine parallel_lib_r(ft,f)

    use mpi

    implicit none

    complex, intent(in), dimension(nj_loc,ni) :: ft
    complex, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k, j1,j2

    if (idebug >= 1) then
       print 9010, iproc,nsend,ni_loc,nj,nj_loc,ni
9010   format('parallel_lib_f: iproc=',i5,' nsend=',i5,                   &
            ' ni_loc=',i5,' nj=',i5,' nj_loc=',i5,' ni=',i5)
    endif

!!$omp  parallel do if (size(fsendr) >= default_size) default(none) &
!!$omp& shared(nproc,iproc,nj_loc,ni_loc) &
!!$omp& private(k,j,j_loc,i,j1,j2) &
!!$omp& shared(ft) &
!!$omp& shared(fsendr)
    do k=1,nproc
       j1 = 1+iproc*nj_loc
       j2 = (1+iproc)*nj_loc
       do j=j1,j2
          j_loc = j-j1+1 
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
         lib_comm, &
         ierr)

  end subroutine parallel_lib_r

  subroutine parallel_lib_rtrans(fin,f)
    ! -----------------------------------------
    ! transpose version of parallel_lib_r(fin,f)
    ! -----------------------------------------
    use mpi

    implicit none

    complex, intent(in), dimension(:,:) :: fin
    complex, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k, j1,j2

!!$omp  parallel do if (size(fsendr) >= default_size) default(none) &
!!$omp& shared(nproc,iproc,nj_loc,ni_loc) &
!!$omp& private(k,j,j_loc,i,j1,j2) &
!!$omp& shared(fin) &
!!$omp& shared(fsendr)
    do k=1,nproc
       j1 = 1+iproc*nj_loc
       j2 = (1+iproc)*nj_loc
       do j=j1,j2
          j_loc = j-j1+1 
          do i=1,ni_loc
             fsendr(i,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_rtrans

  !=========================================================

  subroutine parallel_lib_r_real(ft,f)

    use mpi

    implicit none

    real, intent(in), dimension(nj_loc,ni) :: ft
    real, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2

!!$omp  parallel do if (size(fsendr_real) >= default_size) default(none) &
!!$omp& shared(nproc,iproc,nj_loc,ni_loc) &
!!$omp& private(k,j,j_loc,i,j1,j2) &
!!$omp& shared(ft) &
!!$omp& shared(fsendr_real)
    do k=1,nproc
       j1 = 1+iproc*nj_loc
       j2 = (1+iproc)*nj_loc
       do j=j1,j2
          j_loc = j-j1+1
          do i=1,ni_loc
             fsendr_real(i,j_loc,k) = ft(j_loc,i+(k-1)*ni_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendr_real, &
         nsend, &
         MPI_DOUBLE_PRECISION,&
         f, &
         nsend, &
         MPI_DOUBLE_PRECISION, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_r_real

  subroutine parallel_lib_rtrans_real(fin,f)

    use mpi

    implicit none

    real, intent(in), dimension(:,:) :: fin
    real, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2

!!$omp  parallel do if (size(fsendr_real) >= default_size) default(none) &
!!$omp& shared(nproc,iproc,nj_loc,ni_loc) &
!!$omp& private(k,j,j_loc,i,j1,j2) &
!!$omp& shared(fin) &
!!$omp& shared(fsendr_real)
    do k=1,nproc
       j1 = 1+iproc*nj_loc
       j2 = (1+iproc)*nj_loc
       do j=j1,j2
          j_loc = j-j1+1
          do i=1,ni_loc
             fsendr_real(i,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendr_real, &
         nsend, &
         MPI_DOUBLE_PRECISION,&
         f, &
         nsend, &
         MPI_DOUBLE_PRECISION, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_rtrans_real

  !=========================================================

  !  parallel_slib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_slib_r -> g(nj_loc,ni) -> f(ni_loc,nj)

  subroutine parallel_slib_init(nn_in,nexch_in,nkeep_in,nsplit_out,comm)

    use mpi

    !-------------------------------------------
    implicit none
    !
    integer, intent(in) :: nn_in
    integer, intent(in) :: nexch_in
    integer, intent(in) :: nkeep_in
    integer, intent(inout) :: nsplit_out
    integer, intent(in) :: comm
    integer :: ierr
    !-------------------------------------------

    slib_comm = comm

    !-------------------------------------------------
    ! Get rank and number of processors from slib_comm 
    ! (which is assumed to already exist).
    !
    call MPI_COMM_RANK(slib_comm,isproc,ierr)
    call MPI_COMM_SIZE(slib_comm,nsproc,ierr)
    !-----------------------------------------------

    nn  = nn_in
    nexch = nexch_in
    nkeep = nkeep_in

    nsplit = 1+(nexch-1)/nn

    nsplit_out = nsplit

  end subroutine parallel_slib_init

  !=========================================================

  subroutine parallel_slib_f(x_in,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nexch) :: x_in
    complex, dimension(nkeep,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep*nsplit*nn) :: xt
    !
    integer :: j
    integer :: ierr
    integer :: icount, ip, irequest, k, source, dest, tag, offset
    integer :: request_list(nsproc+nsproc)
    logical :: isok
    !-------------------------------------------------------
    if (idebug >= 1) then
       print 9010, isproc,nkeep,nexch,nsplit,nn
9010   format('parallel_slib_f: isproc=',i5,' nkeep=',i5,                  &
            ' nexch=',i5,' nsplit=',i5,' nn=',i5)
    endif

    isok = (nsproc <= nn)
    if (.not.isok) then
       print 9020,isproc,nsproc,nsplit,nn
9020   format('parallel_slib_f:invalid nn, isproc=',i5,                   &
            ' nsproc=',i5,' nsplit=',i5,' nn=',i5)
    endif
    call assert(isok,'parallel_slib_f: invalid nn=',nn)

    if (use_alltoall_slib_f) then
       do j=1,nexch
          x(:,j) = x_in(:,j)
       enddo

       do j=nexch+1,nsplit*nn
          x(:,j) = czero
       enddo

       call MPI_ALLTOALL(x, &
            nkeep*nsplit, &
            MPI_DOUBLE_COMPLEX, &
            xt, &
            nkeep*nsplit, &
            MPI_DOUBLE_COMPLEX, &
            slib_comm, &
            ierr)
    else


       !   --------------
       !   post MPI_Irecv
       !   --------------
       do k=1,nsproc
          source = (k-1)
          dest = isproc
          tag = mod( source + dest*nsproc,TAG_UB)
          offset = 1 + (k-1)*(nkeep*nsplit)
          icount = nkeep*nsplit
          call MPI_Irecv( xt(offset),icount,MPI_DOUBLE_COMPLEX,              &
               source, tag, slib_comm, irequest,ierr)
          call assert(ierr.eq.MPI_SUCCESS,                                   &
               'parallel_slib_f: MPI_Irecv return ierr',ierr)
          request_list(k) = irequest
       enddo

       !    ---------------------
       !    copy data into buffer
       !    ---------------------
       do j=1,nexch
          x(:,j) = x_in(:,j)
       enddo

       do j=nexch+1,nsplit*nn
          x(:,j) = czero
       enddo

       !   -----------------
       !   perform MPI_Isend
       !   -----------------

       do k=1,nsproc
          source = isproc
          dest = (k-1)
          tag = mod(source+dest*nsproc,TAG_UB)
          ip = 1 + (k-1)*nsplit
          icount = nkeep*nsplit
          call MPI_Isend( x(1,ip),icount,MPI_DOUBLE_COMPLEX,                 &
               dest, tag, slib_comm, irequest, ierr)
          call assert(ierr.eq.MPI_SUCCESS, &
               'parallel_slib_f: MPI_Isend return ierr',ierr)
          request_list(nsproc+k) = irequest
       enddo

       call MPI_Waitall(2*nsproc,request_list,MPI_STATUSES_IGNORE,ierr)
       call assert(ierr.eq.MPI_SUCCESS,                                     &
            'parallel_slib_f:MPI_Waitall return ierr',ierr)


    endif

  end subroutine parallel_slib_f

  !=========================================================

  subroutine parallel_slib_r(xt,x_out)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    !complex, intent(in), dimension(nkeep,nsplit,nn) :: xt
    complex, intent(in), dimension(nkeep*nsplit*nn) :: xt
    complex, intent(inout), dimension(nkeep,nexch) :: x_out
    complex, dimension(nkeep,nsplit*nn) :: x
    !
    integer :: j
    integer :: ierr
    integer :: source, dest, tag, icount, ip, irequest
    integer :: request_list(nsproc+nsproc)
    logical :: isok
    !-------------------------------------------------------

    isok = (nsproc <= nn)
    if (.not.isok) then
       print 9020,isproc,nkeep,nsplit,nexch,nn
9020   format('parallel_slib_r: invalid nn, isproc=',i5,                  &
            ' nkeep=',i5,' nsplit=',i5,' nexch=',i5,' nn=',i5)
    endif
    call assert(isok,'parallel_slib_r: invalid nn=',nn)

    if (use_alltoall_slib_r) then

       call MPI_ALLTOALL(xt, &
            nkeep*nsplit, &
            MPI_DOUBLE_COMPLEX, &
            x, &
            nkeep*nsplit, &
            MPI_DOUBLE_COMPLEX, &
            slib_comm, &
            ierr)

       do j=1,nexch
          x_out(:,j) = x(:,j)
       enddo

    else

       !   --------------
       !   post MPI_Irecv
       !   --------------

       do j=1,nsproc
          source = (j-1)
          dest = isproc
          tag = mod(source + dest*nsproc,TAG_UB)
          ip = 1 + (j-1)*nsplit
          icount = nkeep*nsplit
          call MPI_Irecv( x(1,ip),icount, MPI_DOUBLE_COMPLEX,                &
               source, tag, slib_comm, irequest, ierr )
          call assert(ierr.eq.MPI_SUCCESS,                                   &
               'parallel_slib_r: MPI_Irecv return ierr',ierr)
          request_list(j) = irequest
       enddo

       !   --------------
       !   post MPI_Isend
       !   --------------

       do j=1,nsproc
          source = isproc
          dest = (j-1)
          tag = mod(source + dest*nsproc,TAG_UB)
          icount = nkeep*nsplit
          ip = 1 + (j-1)*(nkeep*nsplit)
          call MPI_Isend( xt(ip),icount,MPI_DOUBLE_COMPLEX,                 &
               dest, tag, slib_comm, irequest, ierr)
          call assert(ierr.eq.MPI_SUCCESS,                                  &
               'parallel_slib_r: MPI_Isend return ierr',ierr)

          request_list(nsproc+j) = irequest
       enddo

       call MPI_Waitall(2*nsproc,request_list,MPI_STATUSES_IGNORE,ierr)
       call assert(ierr.eq.MPI_SUCCESS,                                      &
            'parallel_slib_r: MPI_Waitall return ierr',ierr)
       do j=1,nexch
          x_out(:,j) = x(:,j)
       enddo

    endif
  end subroutine parallel_slib_r

end module parallel_lib
