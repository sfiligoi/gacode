!-----------------------------------------------------------------
! cgyro_parallel_lib.F90
!
! PURPOSE:
!  Routines to manage the two types (lib,slib) of ALLTOALL
!  communication required by CGYRO.  This includes both GPU
!  and non-GPU options, as well as GPUDirect MPI options.
!-----------------------------------------------------------------

module parallel_lib

  implicit none

  ! lib

  integer :: nproc,iproc
  integer, private :: ni,nj
  integer, private :: ni_loc
  integer, private :: nj_loc
  integer, private :: lib_comm
  integer, private :: nsend
  real, dimension(:,:,:), allocatable, private :: fsendr_real

  ! (expose these)
  complex, dimension(:,:,:), allocatable :: fsendf
  complex, dimension(:,:,:), allocatable :: fsendr

  ! slib

  integer, private :: nsproc,isproc
  integer, private :: nn
  integer, private :: nexch
  integer, private :: nkeep
  integer, private :: nsplit
  integer, private :: slib_comm

  integer, private, parameter :: default_size = 1

contains

  !=========================================================
  !  parallel_lib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_lib_r -> g(nj_loc,ni) -> f(ni_loc,nj)
  !=========================================================

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
    if (.not. allocated(fsendr_real)) allocate(fsendr_real(ni_loc,nj_loc,nproc))

!$acc enter data create(fsendf,fsendr)

  end subroutine parallel_lib_init

  !=========================================================

  subroutine parallel_lib_f_i_do(ft)

    use mpi

    implicit none

    complex, intent(inout), dimension(nj_loc,ni) :: ft
    integer :: ierr

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_f_i_do

  !=========================================================

  subroutine parallel_lib_f_i_do_gpu(ft)

    use mpi

    implicit none

    complex, intent(inout), dimension(nj_loc,ni) :: ft
    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(fsendf)
#else
!$acc host_data use_device(fsendf,ft)
#endif

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(ft)
#else
!$acc end host_data
#endif

  end subroutine parallel_lib_f_i_do_gpu

  !=========================================================

  subroutine parallel_lib_r_do(f)

    use mpi

    implicit none

    complex, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_r_do

  !=========================================================

  subroutine parallel_lib_r_do_gpu(f)

    use mpi

    implicit none

    complex, intent(inout), dimension(ni_loc,nj) :: f

    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(fsendr)
#else
!$acc host_data use_device(fsendr,f)
#endif

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(f)
#else
!$acc end host_data
#endif

  end subroutine parallel_lib_r_do_gpu

  !=========================================================

  subroutine parallel_lib_r_pack(ft)

    use mpi

    implicit none

    complex, intent(in), dimension(nj_loc,ni) :: ft
    integer :: j_loc,i,j,k,j1,j2

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do if (size(fsendr) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc) &
!$omp& private(j,j_loc,i) &
!$omp& shared(ft,fsendr)
    do k=1,nproc
       do j=j1,j2
          j_loc = j-j1+1 
          do i=1,ni_loc
             fsendr(i,j_loc,k) = ft(j_loc,i+(k-1)*ni_loc) 
          enddo
       enddo
    enddo

  end subroutine parallel_lib_r_pack

  !=========================================================

  subroutine parallel_lib_rtrans_pack(fin)
! -----------------------------------------
! transpose version of parallel_lib_r_pack(fin)
! -----------------------------------------
    use mpi

    implicit none

    complex, intent(in), dimension(:,:) :: fin
    integer :: j_loc,i,j,k,j1,j2

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do if (size(fsendr) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc) &
!$omp& private(j,j_loc,i) &
!$omp& shared(fin,fsendr)
    do k=1,nproc
       do j=j1,j2
          j_loc = j-j1+1 
          do i=1,ni_loc
             fsendr(i,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc) 
          enddo
       enddo
    enddo

  end subroutine parallel_lib_rtrans_pack

  !=========================================================

  subroutine parallel_lib_rtrans_pack_gpu(fin)
! -----------------------------------------
! transpose version of parallel_lib_r_pack(fin)
! -----------------------------------------

    use mpi

    implicit none

    complex, intent(in), dimension(:,:) :: fin
    integer :: j_loc,i,j,k,j1,j2

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc
!$acc parallel loop collapse(3) independent private(j_loc) &
!$acc&         present(fsendr,fin) default(none)
    do k=1,nproc
       do j=j1,j2
          do i=1,ni_loc
             j_loc = j-j1+1
             fsendr(i,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc)
          enddo
       enddo
    enddo

  end subroutine parallel_lib_rtrans_pack_gpu

  !=========================================================

  subroutine parallel_lib_rtrans(fin,f)
! -----------------------------------------
! transpose version of parallel_lib_r(fin,f)
! -----------------------------------------
    use mpi

    implicit none

    complex, intent(in), dimension(:,:) :: fin
    complex, intent(inout), dimension(ni_loc,nj) :: f

    call parallel_lib_rtrans_pack(fin)
    call parallel_lib_r_do(f)

  end subroutine parallel_lib_rtrans

  !=========================================================

  subroutine parallel_lib_rtrans_real(fin,f)

    use mpi

    implicit none

    real, intent(in), dimension(:,:) :: fin
    real, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do if (size(fsendr_real) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc) &
!$omp& private(j,j_loc,i) &
!$omp& shared(fin,fsendr_real)
    do k=1,nproc
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

  ! exchange indexes
  ! This is always done using CPU memory
  subroutine parallel_slib_f_idxs(nels,x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels
    integer, intent(in), dimension(nels*nn) :: x
    integer, intent(inout), dimension(nels,nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

    call MPI_ALLTOALL(x, &
         nels, &
         MPI_INTEGER, &
         xt, &
         nels, &
         MPI_INTEGER, &
         slib_comm, &
         ierr)

  end subroutine parallel_slib_f_idxs

  ! exchange indexes
  ! This is always done using CPU memory
  subroutine parallel_slib_r_idxs(nels,xt,x)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels
    integer, intent(in), dimension(nels,nn) :: xt
    integer, intent(inout), dimension(nels,nn) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

    call MPI_ALLTOALL(xt, &
         nels, &
         MPI_INTEGER, &
         x, &
         nels, &
         MPI_INTEGER, &
         slib_comm, &
         ierr)

  end subroutine parallel_slib_r_idxs

  ! This is always done using CPU memory
  subroutine parallel_slib_cpu_maxval_int(val)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(inout) :: val
    !
    integer :: ierr
    integer, dimension(1) :: valin,valout
    !-------------------------------------------------------

    valin(1) = val

    call MPI_ALLREDUCE(valin, &
         valout, &
         1, &
         MPI_INTEGER, &
         MPI_MAX, &
         slib_comm, &
         ierr)

    val = valout(1)

  end subroutine parallel_slib_cpu_maxval_int

  ! test an async req, to progress async operations
  subroutine parallel_slib_test(req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(inout) :: req
    !
    logical :: iflag
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

    call MPI_REQUEST_GET_STATUS(req, iflag, istat, ierr)
    ! we discard all the outputs... it is just a way to progress async mpi

  end subroutine parallel_slib_test

!=========================================================

  subroutine parallel_slib_f_nc(x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nsplit,nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

    call MPI_ALLTOALL(x, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

  end subroutine parallel_slib_f_nc

  subroutine parallel_slib_f_nc_async(x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nsplit,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------
#ifdef _OPENACC
!$acc data present(x,xt)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(x)
#else
!$acc host_data use_device(x,xt)
#endif
#endif

   call MPI_IALLTOALL(x, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifdef _OPENACC
#ifdef DISABLE_GPUDIRECT_MPI
   !do nothing yet, async
#else
!$acc end host_data
#endif

!$acc end data
#endif

  end subroutine parallel_slib_f_nc_async

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_nc_wait(x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nsplit,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifdef _OPENACC
!$acc data present(xt)
#endif

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef _OPENACC
#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(xt)
#endif

!$acc end data
#endif

  end subroutine parallel_slib_f_nc_wait

 !=========================================================

  subroutine parallel_slib_r_nc (xt,x)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nsplit,nn) :: xt
    complex, intent(inout), dimension(nkeep,nsplit*nn) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef _OPENACC
!$acc data present(xt,x)
#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(xt)
#else
!$acc host_data use_device(xt,x)
#endif
#endif

    call MPI_ALLTOALL(xt, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         x, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef _OPENACC
#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(x)
#else
!$acc end host_data
#endif
!$acc end data
#endif

  end subroutine parallel_slib_r_nc

!=========================================================

  ! Automaticallt use CPU or GPU version, based on presence of ACC
  subroutine parallel_slib_f_fd(nels1,nels2,nels3,x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(in), dimension(nels1,nels2,nels3,nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef _OPENACC
!$acc data present(xt,x)
!$acc host_data use_device(xt,x)
#endif

    call MPI_ALLTOALL(x, &
         nels1*nels2*nels3, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nels1*nels2*nels3, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef _OPENACC
!$acc end host_data
!$acc end data
#endif

  end subroutine parallel_slib_f_fd

  subroutine parallel_slib_f_fd_async(nels1,nels2,nels3,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(in), dimension(nels1,nels2,nels3,nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef _OPENACC
!$acc data present(xt,x)
!$acc host_data use_device(xt,x)
#endif

   call MPI_IALLTOALL(x, &
         nels1*nels2*nels3, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nels1*nels2*nels3, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifdef _OPENACC
!$acc end host_data
!$acc end data
#endif

  end subroutine parallel_slib_f_fd_async

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_fd_wait(nels1,nels2,nels3,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(in), dimension(nels1,nels2,nels3,nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifdef _OPENACC
!$acc data present(x,xt)
#endif

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef _OPENACC
!$acc end data
#endif

  end subroutine parallel_slib_f_fd_wait

!=========================================================

  ! x is logically (nels,nn)
  subroutine parallel_slib_distribute_real(nels,x)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels
    real, intent(inout), dimension(*) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef _OPENACC
!$acc host_data use_device(x)
#endif

    call MPI_ALLTOALL(MPI_IN_PLACE, &
         nels, &
         MPI_DOUBLE, &
         x, &
         nels, &
         MPI_DOUBLE, &
         slib_comm, &
         ierr)

#ifdef _OPENACC
!$acc end host_data
#endif

  end subroutine parallel_slib_distribute_real

!=========================================================

  subroutine parallel_lib_nj_loc(nj_loc_in)

    integer, intent(inout) :: nj_loc_in

    nj_loc_in = nj_loc
    
  end subroutine parallel_lib_nj_loc

  subroutine parallel_lib_clean
    implicit none
    
    if(allocated(fsendf)) then
!$acc exit data delete(fsendf)
       deallocate(fsendf)
    endif

    if(allocated(fsendr)) then
!$acc exit data delete(fsendr)
       deallocate(fsendr)
    endif

    if(allocated(fsendr_real))     deallocate(fsendr_real)
    
  end subroutine parallel_lib_clean
  
end module parallel_lib

