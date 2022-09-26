program cgyro_test_alltotall
  use mpi
  implicit none

  integer :: supported, i_err
  integer :: i_proc,n_proc
  integer :: i,n,j

  integer :: splitkey, i_group
  integer :: NEW_COMM_1

  integer :: total_bufsize
  integer :: bufsize
  real, dimension(:,:), allocatable :: t1,t2

  real:: m1,m2,m3

  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,i_err)

  if (i_err /= 0) then
     write(*,*) 'ERROR: MPI_INIT_THREAD failed'
     call exit(1)
  endif

  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)

  if (n_proc .LT. 2) then
     write(*,*) i_proc, 'ERROR: n_proc < 2', n_proc
     call exit(1)
  endif

#ifdef DISABLE_MPI_SPLIT
  NEW_COMM_1 = MPI_COMM_WORLD
  if (i_proc==0) then
    write(*,*) "INFO: Using MPI_COMM_WORLD"
  endif
#else
  splitkey =  i_proc
  i_group = 1
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       i_group,&
       splitkey,&
       NEW_COMM_1, &
       i_err)
  if (i_err /= 0) then
     write(*,*) i_proc, 'ERROR: MPI_COMM_SPLIT failed'
     call exit(1)
  endif

  call MPI_COMM_RANK(NEW_COMM_1,i_proc,i_err)
  call MPI_COMM_SIZE(NEW_COMM_1,n_proc,i_err)

  if (n_proc .LT. 2) then
     write(*,*) i_proc, 'ERROR: split n_proc < 2', n_proc
     call exit(1)
  endif
#endif


  if (modulo(n_proc, 2) /= 0) then
     write(*,*) i_proc, 'ERROR: split n_proc not even', n_proc
     call exit(1)
  endif

  !-------------------
  ! Initialize buffers
  !-------------------

  total_bufsize = 8*6*9*64*1024
  bufsize = total_bufsize/n_proc

  if (modulo(total_bufsize, n_proc*2) /=  0) then
    write(*,*) i_proc, "ERROR: total_bufsize not dividable by n_proc*2", total_bufsize, n_proc*2
    call exit(1)
  endif

  if (i_proc==0) then
    write(*,*) "n_proc=",n_proc,"bufsize=",bufsize
  endif

  allocate(t1(bufsize,n_proc))
  allocate(t2(bufsize,n_proc))

!$omp parallel do collapse(2)
  do n=1, n_proc
    do i=1,bufsize/2
      t1(i,n) = 100.0*i_proc + i + 0.01*n
      t2(i,n) = 0.0
    enddo
  enddo

!$acc enter data copyin(t2)
!$acc enter data copyin(t1)

#ifdef _OPENACC
!$acc parallel loop collapse(2) gang vector independent present(t1,t2)
#else
!$omp parallel do collapse(2)
#endif
do n=1, n_proc
    do i=(1+bufsize/2),bufsize
      t1(i,n) = 100.0*i_proc + i + 0.01*n
      t2(i,n) = 0.0
    enddo
  enddo

!$acc update host(t1)
  if (t1(1+bufsize/2,n_proc/2) /= (100.0*i_proc + 1+bufsize/2 + 0.01*n_proc/2)) then
    write(*,*) i_proc, "ERROR: Initial t1(bufsize,n_proc) validation failed!", &
               t1(1+bufsize/2,n_proc/2), (100.0*i_proc + 1+bufsize/2+ 0.01*n_proc/2)
    call exit(1)
  endif

  m1 = MPI_Wtime()

!$acc host_data use_device(t1,t2)
    call MPI_ALLTOALL(t1, &
         bufsize, &
         MPI_DOUBLE, &
         t2, &
         bufsize, &
         MPI_DOUBLE, &
         NEW_COMM_1, &
         i_err)
!$acc end host_data

!$acc update host(t2)
  if (t2(1+bufsize/2,n_proc/2) /= (100.0*(n_proc/2-1) + 1+bufsize/2 + 0.01*(i_proc+1))) then
    write(*,*) i_proc, "ERROR: Initial t2(bufsize,n_proc) validation failed!", &
               t2(1+bufsize/2,n_proc/2), (100.0*(n_proc/2-1) + 1+bufsize/2+ 0.01*(i_proc+1))
    call exit(1)
  endif

  ! -------------------
  ! Now run many times for benchmark purposes
  m2 = MPI_Wtime()

  do j=1,10000
    if ((i_proc==0) .AND. (modulo(j,1000) == 0)) then
      m1 = MPI_Wtime()
      write(*,*) "j=", j, "dt=", m1-m2
      call flush(6)
    endif
#ifdef _OPENACC
!$acc parallel loop gang vector independent present(t2) private(n)
#else
!$omp parallel do private(n)
#endif
    do i=1,bufsize
!$acc loop seq
      do n=2, n_proc
        t2(i,n) = t2(i,n) +  t2(i,1)
      enddo
    enddo

!$acc host_data use_device(t1,t2)
      call MPI_ALLTOALL(t2, &
           bufsize, &
           MPI_DOUBLE, &
           t1, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_1, &
           i_err)
!$acc end host_data

#ifdef _OPENACC
!$acc parallel loop gang vector independent present(t1) private(n)
#else
!$omp parallel do private(n)
#endif
    do i=1,bufsize
!$acc loop seq
    do n=1, n_proc-1
        t1(i,n) = t1(i,n) + t1(i,n_proc)
      enddo
    enddo

!$acc host_data use_device(t1,t2)
      call MPI_ALLTOALL(t1, &
           bufsize, &
           MPI_DOUBLE, &
           t2, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_1, &
           i_err)
!$acc end host_data

  enddo
  m3 = MPI_Wtime()

!$acc exit data delete(t1)

!$acc exit data copyout(t2)
  write(*,*) i_proc, "Success", t2(1,i_proc+1), "dt=", m3-m2

end program

