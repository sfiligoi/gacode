!------------------------------------------------------------
! make_MPI_grid.f90
!
! PURPOSE:
!  This routine manages most of the MPI initialization chore.
!------------------------------------------------------------

subroutine gyro_mpi_grid

  use mpi
  use gyro_globals
  use gyro_pointers

  !----------------------------------------
  implicit none
  !
  integer, external :: parallel_dim
  integer :: i_group_mumps
  integer :: splitkey
  !----------------------------------------

  !------------------------------------------------
  ! Check for grid validity
  !
  if (modulo(n_proc,n_n) /= 0) then
     call catch_error('ERROR: bad processor count.')
  endif
  !
  if (linsolve_method == 3 .and. modulo(n_energy,n_proc) /= 0) then
     call catch_error('ERROR: bad processor count for LINSOLVE_METHOD=3.')
  endif
  !------------------------------------------------

  !-------------------------------
  ! Assign subgroup dimensions:
  !
  n_proc_2 = n_n
  n_proc_1 = n_proc/n_n
  n_n_1    = 1
  !-------------------------------

  !------------------------------------------------
  ! Local group indices:
  !
  i_group_1 = i_proc/n_proc_1
  i_group_2 = modulo(i_proc,n_proc_1)
  !
  !------------------------------------------------

  !-----------------------------------------------------------
  ! Split up GYRO_COMM_WORLD into groups and adjoint:
  !
  !             NEW_COMM_1  and  NEW_COMM_2
  !
  splitkey = i_proc

  call MPI_COMM_SPLIT(GYRO_COMM_WORLD,&
       i_group_1,& 
       splitkey,&
       NEW_COMM_1, &
       i_err)
  if (i_err /= 0) then
     print *,'NEW_COMM_1 creation status',i_err
     call MPI_FINALIZE(i_err)
     stop
  endif

  ! Local adjoint Group number

  call MPI_COMM_SPLIT(GYRO_COMM_WORLD,&
       i_group_2,& 
       splitkey,&
       NEW_COMM_2, &
       i_err)
  if (i_err /= 0) then
     print *,'NEW_COMM_2 creation status',i_err
     call MPI_FINALIZE(i_err)
     stop
  endif
  !
  call MPI_COMM_RANK(NEW_COMM_1,i_proc_1,i_err)
  call MPI_COMM_RANK(NEW_COMM_2,i_proc_2,i_err)
  !
  !------------------------------------------------------------

  ! Dimensions used for SSUB library:
  ! (note that n_proc_2 == n_n).

  nv1_SSUB    = n_x
  nv2_SSUB    = parallel_dim(n_energy*n_lambda,n_proc_1)
  msplit_SSUB = 1+(n_stack*nv2_SSUB-1)/n_n

  call SSUB_init(n_n,n_stack,nv1_SSUB,nv2_SSUB,NEW_COMM_2)  

  if (debug_flag == 1) then
     print *,'----------------------------------------'
     print '(t2,a,2(i3,1x))','i_proc,i_group_1 : ',i_proc,i_group_1
     print *,'----------------------------------------'
  endif

  !--------------------------------------------------------------
  ! Create MUMPS_COMM if using MUMPS:
  !
  if (sparse_method == 2) then
     i_group_mumps = i_group_2/(min(n_proc_1,n_mumps_max))
     call MPI_COMM_SPLIT(NEW_COMM_1,&
          i_group_mumps,& 
          i_proc_1,&
          MUMPS_COMM, &
          i_err)
     if (i_err /= 0) then
        print *,'MUMPS_COMM creation status',i_err
        call MPI_FINALIZE(i_err)
        stop
     endif
  endif
  !--------------------------------------------------------------  

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_mpi_grid done]'
  endif

end subroutine gyro_mpi_grid
