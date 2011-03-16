!------------------------------------------------
! catch_blowup.f90
!
! PURPOSE:
!  This routine traps a blow-up of the numerical
!  solution by monitoring the potential amplitude. 
!------------------------------------------------

subroutine catch_blowup

  use mpi
  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  real, parameter :: phi_limit = 0.5
  real :: phi_local
  real :: phi_max
  !---------------------------------------------------

  !------------------------------------------------
  ! Halt on large field
  !
  phi_local = maxval(abs(field_blend))
  !
  call MPI_ALLREDUCE(phi_local, &
       phi_max, &
       1, &
       MPI_DOUBLE_PRECISION, &
       MPI_MAX, &
       NEW_COMM_2, &
       i_err)
  
  if (phi_max > phi_limit .and. nonlinear_flag == 1) then
     call set_exit_status('field blowup',1)
  endif
  !
  !------------------------------------------------

end subroutine catch_blowup
