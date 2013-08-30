!-------------------------------------------------------------
! gyro_initialize_timestep.f90
!
! PURPOSE:
!  Set number of steps based on maximum simulation time
!  (time_max).  For single-mode runs, the run will be
!  terminated early if freq_err drops below a given 
!  criteria (freq_tol).
!
!  Also, set integrator parameters and initialize timers.
!-------------------------------------------------------------

subroutine gyro_initialize_timestep

  use gyro_globals

  !------------------------
  implicit none
  !
  real :: new_max
  !------------------------

  !----------------------------------------------------
  ! Mikkelsen's time-step tolerance function
  !
  ! Note: freq_err is updated, and freq_err test 
  ! executed, ONLY for single-mode runs.
  !
  if (time_max < 0.0) then

     ! Initialize status to unconverged.

     call gyro_set_exit_status('unconverged',1)

     if (abs(time_max) > 1.0) then

        ! Example:
        !
        ! time_max=-1000.001 then error tolerance is 
        ! 0.001 and time_max -> 1000.0

        new_max  = int(abs(time_max))
        freq_tol = abs(time_max) - new_max
        time_max = real(new_max)
        freq_err = 1.0


     else

        ! Example:
        !
        ! time_max=-0.001 then error tolerance is 
        ! 0.001 and time_max -> 50.0

        freq_tol = abs(time_max)
        freq_err = 1.0
        time_max = 50.0

     endif

  else

     ! Initialize status:

     call gyro_set_exit_status('clean exit',0)

     freq_tol = 1e-10
     freq_err = 1.0

  endif

  nstep = nint(time_max/dt)
  !----------------------------------------------------

  !----------------------------------------------------
  ! Compute integrator substep parameter:
  !
  if (integrator_method > 1 .and. electron_method == 2) then

     n_substep = integrator_method-1

  else 

     ! We will always have n_substep=0 for explicit 
     ! integration.

     n_substep = 0

  endif

  a_SDIRK = 0.5
  !----------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_initialize_timestep done]'
  endif

end subroutine gyro_initialize_timestep
