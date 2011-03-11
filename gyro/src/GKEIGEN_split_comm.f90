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

  integer :: ierr

  include 'mpif.h'

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
  if (linsolve_method /= 2) then
    gkeigen_proc_mult = 1
    call send_line('gkeigen_proc_mult reset to 1 (non-GKEIGEN run)')
  endif

  if (n_proc < gkeigen_proc_mult**2) then
    gkeigen_transpose_flag = 1
    gkeigen_j_set = Int(i_proc*gkeigen_proc_mult/n_proc)
  else
    gkeigen_transpose_flag = 0
    gkeigen_j_set = Mod(i_proc,gkeigen_proc_mult)
  endIf

  call MPI_COMM_SPLIT(GYRO_COMM_WORLD,&
       gkeigen_j_set,&
       i_proc,&
       GKEIGEN_J_SUBSET, &
       i_err)
  !
  GYRO_COMM_WORLD = GKEIGEN_J_SUBSET

  call MPI_COMM_RANK(MPI_COMM_WORLD,j_proc_tot,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc_tot,i_err)

end subroutine GKEIGEN_split_comm
