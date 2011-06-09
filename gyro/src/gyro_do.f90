!-----------------------------------------------------------------
! gyro_do.f90
!
! PURPOSE:
!  Subroutinized main gyro program. 
!-----------------------------------------------------------------

subroutine gyro_do

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants
  use GEO_interface

  !--------------------------------------
  implicit none
  !
  logical :: rfe
  !--------------------------------------

  ! Begin with clean exit status
  !
  ! gyro_exit_status=0 : unset
  ! gyro_exit_status=1 : error
  ! gyro_exit_status=2 : clean
  !
  gyro_exit_status  = 0
  gyro_exit_message = 'unset'

  ! Prepend path:
  runfile  = trim(path)//trim(baserunfile)
  precfile = trim(path)//trim(baseprecfile)

  if (baserunfile == 'out.gyro.run')  then
     if (i_proc==0 .AND. output_flag==1) THEN
        inquire(file=trim(runfile),exist=rfe)
        if (.not.rfe) then
           open(unit=99,file=trim(runfile),status='unknown')
           close(99)
        endif
     ENDIF
  endif

  if (i_proc==0 .and. gkeigen_j_set==0) print *,runfile

  CPU_0 = 0.0
  CPU_1 = 0.0
  CPU_2 = 0.0
  CPU_3 = 0.0
  CPU_4 = 0.0
  CPU_5 = 0.0
  CPU_6 = 0.0
  CPU_7 = 0.0

  call proc_time(CPU_0)

  !-------------------------
  ! Early initializations:
  !
  total_memory  = 0.0
  alltime_index = 0
  !-------------------------

  !----------------------------------------------------------------
  ! If running in eigensolve mode.
  !
  if (linsolve_method == 2) then
     if (nonlinear_flag == 1) then
        if (i_proc==0 .and. gkeigen_j_set==0) then
           print *, "Eigensolver unavailable in nonlinear mode."
        endif
        stop 
     else
        if (i_proc==0 .and. gkeigen_j_set==0) then
           print *, "GYRO is running in eigensolve mode."
        endif
        eigensolve_restart_flag = restart_method
        restart_method = 0
        if (electron_method /= 1) then
           electron_method = 4
        endif
     endif
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  if (debug_flag == 1) then
     !  Dump the global variables that can be read
     call gyro_dump_input
     !  Dump the interface variables for comparison
     call gyro_dump_interface
  endif
  !----------------------------------------------------------------

  !------------------------------------------------------------
  ! The order of these routines is critical:
  !
  ! Sort through and check all combinations of operational modes
  !
  call gyro_select_methods
  !
  ! Set parameters connected with timestepping.
  !
  call gyro_initialize_timestep
  !
  ! Generate theta grid dimensions (no operators yet).
  !
  call make_theta_grid
  !
  ! Parallel setup 
  !
  if (gyrotest_flag == 0) then
     call gyro_mpi_grid
  else
     n_n_1     = 1
     i_group_1 = 0
  endif
  !
  ! Read, generate or otherwise construct equilibrium profiles.  If experimental 
  ! profiles are used, GEO will be allocate/deallocated with all settings 
  ! determined in EXPRO.
  !
  call gyro_alloc_profile_sim(1)
  call gyro_profile_init
  !
  ! Construct pitch-angle (lambda) and energy grids
  ! and integration weights using Gauss-Legendre rules:
  !
  call gyro_alloc_velocity(1)
  call gyro_alloc_orbit(1)
  !
  ! Set geometry (GEO) library control variables
  !
  GEO_nfourier_in = n_fourier_geo 
  GEO_model_in    = geometry_method
  GEO_signb_in    = -btccw
  call GEO_alloc(1)
  !
  ! Lambda (pitch-angle) weights (GEO needed again, so just reallocate)
  call gyro_lambda_grid
  !
  ! Energy weights
  !
  call energy_integral(n_energy,energy_max,n_kinetic,energy,w_energy)
  !
  ! Compute myriad arrays (for example, map between poloidal
  ! angle and orbit time) for poloidal discretization.
  !
  call gyro_banana_operators
  !
  ! Total velocity-space-tau weights.
  !
  call gyro_set_phase_space(trim(path)//'out.gyro.phase_space',1)
  !
  do i=1,n_x
     call gyro_to_geo(i)
     if (i_proc == 0 .and. i == ir_norm .and. debug_flag == 1) then
        call GEO_write(trim(path)//'gyro_geo_diagnostic.out',1)
     endif
  enddo
  !
  ! Generate geometry-dependent factors using model or Miller equilibrium:
  !
  call make_geometry_arrays
  !
  if (gyrotest_flag == 0) then

     call gyro_set_pointer_dim
     call gyro_alloc_distrib(1)

     ! Make pointers for use with parallelization scheme.
     call gyro_set_pointers

     ! Compute drift and diamagnetic frequency coefficients
     ! for GKE solution.
     call gyro_omegas

     call proc_time(CPU_1)

     ! Generate discrete Fourier and Fourier-sine transform 
     ! matrices for spectral Poisson and adaptive source methods, 
     ! compute radial derivative and gyroaverage operators,
     ! and then print them.

     call gyro_radial_operators
     call proc_time(CPU_2)
     !
     ! Precomputation of arrays which depend on blending 
     ! coefficients.  These are used in the Maxwell solves.
     call gyro_set_blend_arrays

     call proc_time(CPU_3)
     if (electron_method == 2) then

        ! Make advection operators for electrons
        call make_implicit_advect(0)

     endif

     call proc_time(CPU_4)
     !
     ! Collision setup (if required; note PNL hook)
     !
     if (collision_flag == 1) then
        call gyro_collision_setup
     endif

     !------------------------------------------------------
     ! FIELD MATRICES:
     !
     ! Manage building of field-solve matrices.  Must have 
     ! built the gyro-operators beforehand:
     !
     ! Explicit (Poisson,Ampere)
     !
     call proc_time(CPU_5)
     if (n_field == 3) then
        call make_poissonaperp_matrix
     else
        call make_poisson_matrix
     endif
     !
     if (n_field > 1) then
        call make_electron_current(0)
        call make_ampere_matrix
     endif
     !
     ! Implicit (Maxwell)
     !
     if (electron_method == 2) then
        call make_maxwell_matrix
     endif
     call proc_time(CPU_6)
     !
     !------------------------------------------------------

     !------------------------------------------------------
     ! Make a few trivial quantities for the nonlinear step:
     !
     if (nonlinear_flag == 1) then
        call gyro_alloc_nl(1)
        call make_nl
     endif
     !
     ! Write information about radial stencils:
     !
     call write_radial_operators(trim(path)//'r_operators.out',1)
     !
     ! Allocate most large arrays:
     !
     call gyro_alloc_big(1)
     call gyro_initialize_arrays
     !
     ! Make attempt to estimate storage:
     !
     call gyro_memory_usage(trim(path)//'alloc.out',30)
     !
     !------------------------------------------------------------

     !------------------------------------------------------------
     ! Open files and write values for t=0:
     !
     call gyro_read_restart
     !------------------------------------------------------------

     call proc_time(CPU_7)

     !---------------------------------------
     ! Write precision data:
     !
     call gyro_write_precision(10,sum(abs(gbflux)))
     !---------------------------------------

  endif

  !---------------------------------------
  ! Dump input parameters runfile
  !
  call gyro_write_input
  !---------------------------------------

  !---------------------------------------------------------------
  ! I/O control for time-independent initial data
  !
  if (io_method == 1) then
     call gyro_write_initdata(&
          trim(path)//'profile_vugyro.out',&
          trim(path)//'units.out',&
          trim(path)//'geometry_arrays.out',1)
  else
     call gyro_write_initdata_hdf5(trim(path)//'out.gyro.initdata.h5')
  endif
  !
  ! Close geometry (GEO) library
  call GEO_alloc(0)
  !---------------------------------------------------------------

  !------------------------------------------------------------
  ! By this point, we are finished with running the
  ! 'gyro -t' (test) mode.
  !
  if (gyrotest_flag == 1) then
     call write_efficiency(trim(path)//'efficiency.out',1)
     call set_exit_status('test complete',2)
     return
  endif
  !------------------------------------------------------------

  if (restart_method < 1) then
     ! Open
     io_control = output_flag*1
  else
     ! Rewind
     io_control = output_flag*3
  endif
  if (gkeigen_j_set == 0) then
     if (io_method == 1) then
        call gyro_write_timedata
     else
        call gyro_write_timedata_hdf5
     endif
  endif

  !-------------------------------------------------
  ! NEW SIMULATION ONLY: write *initial conditions*
  !
  if (restart_method /= 1) then
     io_control = output_flag*2
     if (gkeigen_j_set == 0) then
        if (io_method == 1) then
           call gyro_write_timedata
        else
           call gyro_write_timedata_hdf5
           if (time_skip_wedge > 0) call gyro_write_timedata_wedge_hdf5
        endif
     endif
  endif
  !--------------------------------------------

  !--------------------------------------------
  ! Start elapsed time.
  !
  if (step == 0) then
     call system_clock(clock_count,clock_rate,clock_max)
     elapsed_time = clock_count*1.0/clock_rate
  endif
  !--------------------------------------------------

  select case (linsolve_method)

  case (1)

     !------------------------------
     ! INITIAL VALUE SOLVER
     !------------------------------

     !-----------------------------------------------
     ! start time-stepping

     do step=1,nstep

        call gyro_fulladvance

        !-------------------------------------
        ! Check for premature exit conditions
        ! 
        call catch_blowup
        call catch_halt_signal
        if (gyro_exit_status > 0) return
        !-------------------------------------

     enddo
     gyro_exit_status = 2

     ! end time-stepping   
     !-----------------------------------------------

  case (2)

     !------------------------------
     ! GK EIGENVALUE SOLVER
     !------------------------------

     ! Determine true range of p_nek_loc for each processor
     n_nek_loc_true = n_nek_1/n_proc_1
     if (i_proc_1 < mod(n_nek_1,n_proc_1) ) then
        n_nek_loc_true = n_nek_loc_true + 1
     endif
     ! Length of total state vector
     h_length = n_kinetic * n_nek_1 * n_x * n_stack
     h_length_loc = n_kinetic * n_nek_loc_true * n_x * n_stack
     h_width_loc = h_length / gkeigen_proc_mult
     h_length_block = h_length_loc / gkeigen_proc_mult
     h_length_block_t = h_width_loc / n_proc
     seq_length = h_length_block * h_width_loc
     seq_length_t = h_length_block_t * h_length_loc

     call GKEIGEN_do                                       

  case (3)

     !------------------------------
     ! FIELD EIGENVALUE SOLVER
     !------------------------------

     call gyro_fieldeigen

  end select

end subroutine gyro_do
