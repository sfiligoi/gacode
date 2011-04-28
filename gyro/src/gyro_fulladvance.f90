!--------------------------------------------------------
! gyro_fulladvance.f90
!
! PURPOSE:
!  Control routine for a single timestep of length dt.
!--------------------------------------------------------

subroutine gyro_fulladvance

  use gyro_globals
  use gyro_pointers

  !----------------------------------------
  implicit none
  !----------------------------------------

  call proc_time(CPU_C_in)

  h_old(:,:,:,:) = h(:,:,:,:)

  h_cap_old2(:,:,:,:) = h_cap_old(:,:,:,:)
  h_cap_old(:,:,:,:)  = h_cap(:,:,:,:)

  field_blend_old2(:,:,:) = field_blend_old(:,:,:)
  field_blend_old(:,:,:)  = field_blend(:,:,:)

  gyro_uv_old2(:,:,:,:,:) = gyro_uv_old(:,:,:,:,:)
  gyro_uv_old(:,:,:,:,:)  = gyro_uv(:,:,:,:,:)

  !----------------------
  ! Start time step cycle
  !----------------------

  !------------------------------------------------------
  ! MANAGE pitch-angle scattering:
  !
  if (collision_flag == 1) then

     select case (collision_method)

     case (1) 

        ! Attempts to update fields; should not be used

        call gyro_collision

     case (2) 

        ! Default

        call gyro_collision_jc

     case (3,4)

        ! For testing only, not for production use

        call gyro_collision_ebelli

     end select

  endif
  !------------------------------------------------------

  call proc_time(CPU_C_out)   
  CPU_C = CPU_C + (CPU_C_out - CPU_C_in)

  !------------------------------------------------------
  ! MANAGE time-integrator
  !
  select case (electron_method) 

  case (1,3,4)

     ! All species explicit

     call timestep_explicit

  case (2)

     ! Drift-kinetic electrons

     call timestep_SSP_322

  end select
  !------------------------------------------------------

  call proc_time(CPU_diag_in)
  CPU_ts = CPU_ts + (CPU_diag_in - CPU_C_out)

  call gyro_timestep_error

  !------------------------------------------------------
  ! MANAGE diagnostics
  !
  ! * Note that (2,3,4) are called once after making 
  !   initial conditions (make_initial_h).
  !
  ! 1. Compute flux-surface average of (phi,a)
  ! 
  call get_field_fluxave
  !
  ! 2. Compute (phi,a) for plotting; we want to call 
  !    get_field_plot every timestep because it is 
  !    averaged (with a weight) over every timestep.
  !
  !    To compute E_parallel
  !
  call gyro_field_time_derivative
  !
  !    then 
  !
  call get_field_plot
  !
  ! 3. Compute (n,T) moments for plotting if user has 
  !    selected OUTPUT_METHOD > 1.
  !
  call gyro_moments_plot 
  if (io_method > 1 .and. time_skip_wedge > 0) call gyro_moments_plot_wedge
  !
  ! 4. Compute (phi,a) at r=r0 for plotting (if user 
  !    has selected FIELD_RO_FLAG=1).
  !
  if (field_r0_flag == 1) call get_field_r0_plot 
  !------------------------------------------------------

  !------------------------------------------------------
  ! Increment micro-step counter
  !
  alltime_index = alltime_index+1
  !------------------------------------------------------

  call proc_time(CPU_diag_mid)
  CPU_diag_a = CPU_diag_a + (CPU_diag_mid - CPU_diag_in)

  !------------------------------------------------------
  ! MANAGE time:
  !
  ! Increment simulation time; this is the actual time 
  ! in the gyrokinetic equation:
  !  
  t_current = t_current+dt
  !------------------------------------------------------

  !-------------------------------------------------------------------
  ! MANAGE data output: 
  !
  if(time_skip_wedge > 0) then
    if (modulo(step,time_skip_wedge) == 0 .and. io_method > 1) then
     call write_hdf5_wedge_timedata(2)
    endif
  endif
  if (modulo(step,time_skip) == 0) then

     ! Counter for number of data output events.

     data_step = data_step+1

     !------------------------------------------------
     ! Compute nonlinear transfer and turbulent energy 
     !------------------------------------------------

     if (n_substep == 0 .and. nonlinear_transfer_flag == 1) then

        rhs(:,:,:,:) = (0.0,0.0)

        call get_nonlinear_advance
        call gyro_nonlinear_transfer

     endif

     ! Reset micro-step counter.

     alltime_index = 0

     ! Main data I/O handler

     call gyro_write_master(2)
     if (io_method > 1) call write_hdf5_timedata(2)

     !--------------------------------------------------
     ! Update diffusivity and flux time-record for TGYRO 
     !--------------------------------------------------

     if (transport_method == 2) then
        diff_vec(:,:,:,data_step+1) = diff(:,:,:)
        gbflux_vec(:,:,:,data_step+1) = gbflux(:,:,:)
     endif

     ! Restart test:

     if (modulo(data_step,restart_data_skip) == 0 &
          .and. restart_method >= 0) then

        call write_restart

     endif

  endif ! modulo(step,time_skip) test
  !-------------------------------------------------------------------

  !----------------------------------------------------
  ! Convergence check for single-n simulation:
  ! freq_err calculated in write_freq
  !  
  if (freq_err < freq_tol) call set_exit_status('converged',2)
  !----------------------------------------------------

  !--------------------------------------
  ! ** proc_time call for CPU_diag_outp in 
  ! previous call to gyro_write_master:
  !--------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print '(t2,2(a,i5,3x))','-> step =',step,'data_step =',data_step
     print *,'**[gyro_fulladvance done]'
  endif
  call proc_time(CPU_diag_out)
  CPU_diag_b = CPU_diag_b + (CPU_diag_out - CPU_diag_mid)

end subroutine gyro_fulladvance
