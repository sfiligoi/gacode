!------------------------------------------
! GKEIGEN_split_com.f90
!
! PURPOSE:
!   Split the MPI communicator into
!   gkeigen_proc_mult separate instances
!   of GYRO_COMM_WORLD.
!------------------------------------------

subroutine GKEIGEN_split_comm

  use gyro_globals
  use mpi

  integer :: ierr

  !----------------------------------------------------------------
  ! Split the communicator into gkeigen_proc_mult parallel instances
  ! of gyro each calculating a fraction of the columns of the
  ! gkeigen matrix.
  !
  ! j_proc_tot : processor rank in MPI_COMM_WORLD
  ! i_proc     : processor rank in GYRO_COMM_WORLD
  ! n_proc     : (total process)/gkeigen_proc_mult
  !               number of processes in each GYRO_COMM_WORLD
  !
  ! If gkeigen_proc_mult > 1, the user must specify a processor 
  ! number that is an integer N such that:
  !
  !        N = (n_proc)*(gkeigen_proc_mult)
  ! and
  !
  !        n_nek = (n_proc)*M, where M is any integer.
  !
  if ((linsolve_method/=2) .AND. (gkeigen_proc_mult/=1)) then
    gkeigen_proc_mult = 1
    call send_line('gkeigen_proc_mult reset to 1 (non-GKEIGEN run)')
  endif

  gkeigen_j_set = Mod(i_proc,gkeigen_proc_mult)

  call MPI_COMM_SPLIT(GYRO_COMM_WORLD,&
       gkeigen_j_set,&
       i_proc,&
       GKEIGEN_J_SUBSET, &
       i_err)
  !
  GYRO_COMM_WORLD = GKEIGEN_J_SUBSET

  !---------------------------------------------------------------
  ! Reset process ranks if MPI_COMM_WORLD was split.
  call MPI_COMM_RANK(GYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(GYRO_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------

  call MPI_COMM_RANK(MPI_COMM_WORLD,j_proc_tot,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc_tot,i_err)

  ! Generate a separate communicator for all process sharing
  ! the same value of i_proc for use in MPI_ALLTOALL.
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
                      i_proc,&
                      j_proc_tot,&
                      GYRO_COMM_UNIPROC,&
                      ierr)

end subroutine GKEIGEN_split_comm
