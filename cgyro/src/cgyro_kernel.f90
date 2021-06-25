!-----------------------------------------------------------------
! cgyro_kernel.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!
! NOTES:
!  This can be called directly using the driver routine cgyro 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using cgyro_sub.
!-----------------------------------------------------------------

subroutine cgyro_kernel
  
  use timer_lib
  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none

  if (test_flag == 1) return

  !---------------------------------------------------------------------------
  !
  ! Time-stepping count
  n_time = nint(max_time/delta_t)

  ! Time-averaged fluxes for this time-stepping segment
  tave_step  = 0
  tave_min   = t_current
  tave_max   = t_current
  cflux_tave = 0.0
  gflux_tave = (0.0,0.0)
  
  do i_time=1,n_time

     call timer_lib_in('TOTAL')

     !------------------------------------------------------------
     ! Time advance
     !
     t_current = t_current+delta_t

     ! Collisionless step: returns new h_x, cap_h_x, fields
     
     select case(delta_t_method)
     case(1)
        call cgyro_step_gk_ck
     case(2)
        call cgyro_step_gk_bs5
     case(3)
        call cgyro_step_gk_v76 
     case default
        ! Normal timestep
        call cgyro_step_gk
     end select
     
     ! Collision step: returns new h_x, cap_h_x, fields
     if (collision_model == 5) then
        call cgyro_step_collision_simple
     else
        call cgyro_step_collision
     endif

     if (shear_method == 1) then
        ! Discrete shift (Hammett) 
        call cgyro_shear_hammett
     endif

     ! field will not be modified in GPU memory for the rest of the loop
!$acc update host(field) async(3)

     if (mod(i_time,print_step) == 0) then
       ! cap_h_c will not be modified in GPU memory for the rest of the loop
!$acc update host(cap_h_c) async(4)
     endif

     call cgyro_source
    !------------------------------------------------------------

     !------------------------------------------------------------
     ! Diagnostics
     !
     ! NOTE: Fluxes are calculated in cgyro_write_timedata

  call timer_lib_in('coll_mem')
  ! wait for fields to be synched into system memory, used by cgyro_error_estimate
!$acc wait(3)
  call timer_lib_out('coll_mem')

     ! Error estimate
     call cgyro_error_estimate
     ! Exit if error too large
     if (error_status > 0) exit
     !------------------------------------------------------------

     !---------------------------------------
     ! IO
     !
     if (mod(i_time,print_step) == 0) then
       call timer_lib_in('coll_mem')
       ! wait for cap_h_c to be synched into system memory, used by cgyro_write_timedata
!$acc wait(4)
       call timer_lib_out('coll_mem')

       call timer_lib_in('io')

       ! Write simulation data
       call cgyro_write_timedata

       call timer_lib_out('io')
     endif

     call timer_lib_in('io')

     ! Write restart data
     call cgyro_write_restart

     call timer_lib_out('io')
     !---------------------------------------

     call timer_lib_out('TOTAL')

     ! Don't wrap timer output in a timer
     if (mod(i_time,print_step) == 0) call write_timers(trim(path)//runfile_timers)

     ! Exit if convergenced
     if (signal == 1) exit

  enddo
  !---------------------------------------------------------------------------

end subroutine cgyro_kernel
