!-----------------------------------------------------------------
! gyro_do.f90
!
! PURPOSE:
!  Subroutinized main gyro program. 
!
! NOTES:
!  << BigScience >> was the legacy name for main dating back to 
!  1999.  The executable name is also BigScience.
!-----------------------------------------------------------------

subroutine gyro_do(skipinit)

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !--------------------------------------
  implicit none
  !
  integer, optional :: skipinit
  logical :: rfe
  !--------------------------------------

  !-------------------------------------
  ! Handling of optional arguments
  !
  if (present(skipinit)) then
     lskipinit = skipinit
  else
     lskipinit = 0
  endif
  !-------------------------------------

  ! Begin with clean exit status
  !
  ! gyro_exit_status=0 : unset
  ! gyro_exit_status=1 : error
  ! gyro_exit_status=2 : clean
  !
  gyro_exit_status  = 0
  gyro_exit_message = 'unset'

  if (lskipinit == 1) then
     goto 100
  endif

  ! Prepend path:
  runfile  = trim(path)//trim(baserunfile)
  precfile = trim(path)//trim(baseprecfile)

  if (baserunfile == 'out.gyro.run')  then
     IF(i_proc==0 .AND. output_flag==1) THEN
        inquire(file=trim(runfile),exist=rfe)
        if (.not.rfe) then
            open(unit=99,file=trim(runfile),status='unknown')
            close(99)
         endif
     ENDIF
  endif



  if ((i_proc==0).AND.(gkeigen_j_set==0)) print *,runfile


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

  !-----------------------------------------------------------
  ! If running in eigensolve mode.
  !
  if (linsolve_method == 2) then
     if (nonlinear_flag == 1) then
        if ((i_proc==0) .and. (gkeigen_j_set==0)) print *, "Eigensolver unavailable in nonlinear mode."
        stop
     else
        if ((i_proc==0) .and. (gkeigen_j_set==0)) print *, "GYRO is running in eigensolve mode."
        eigensolve_restart_flag = restart_method
        restart_method = 0
        if (electron_method /= 1) then
           electron_method = 4
        endif
     endif
  endif
  if (debug_flag == 1) then
  !----------------------------------------------------------
  !sv dump the global variables that can be read
  call gyro_dump_input
  !sv dump the inteface variables for comparison
  call gyro_dump_interface
  endif


  !------------------------------------------------------------
  ! The order of these routines is critical:
  !
  ! (1) Sort through and check all combinations of 
  !     operational modes
  !
  call gyro_select_methods
  !
  ! (2) Set parameters connected with timestepping.
  !
  call initialize_timestep
  !
  ! (3) Generate theta grid dimensions (no operators yet).
  !
  call make_theta_grid
  !
  ! (4) Parallel setup 
  !
  if (gyrotest_flag == 0) then
     call gyro_mpi_grid
  else
     n_n_1       = 1
     i_group_1   = 0
  endif
  !
  ! (5) Read, generate or otherwise construct equilibrium
  !     profiles
  !
  call gyro_alloc_profile_sim(1)
  call gyro_profile_init
  !
  ! Construct pitch-angle (lambda) and energy grids
  ! and integration weights using Gauss-Legendre rules:
  call gyro_alloc_velocity(1)
  call gyro_alloc_orbit(1)
  !
  ! Lambda (pitch-angle) weights 
  call gyro_lambda_grid
  !
  ! Energy weights
  call energy_integral(n_energy,energy_max,n_kinetic,energy,w_energy)
  !
  ! Compute myriad arrays (for example, map between poloidal
  ! angle and orbit time) for poloidal discretization.
  call gyro_banana_operators
  !
  ! Total velocity-space-tau weights.
  call gyro_set_phase_space(trim(path)//'out.gyro.phase_space',1)
  !
  ! Generate geometry-dependent factors using model or
  ! Miller equilibrium:
  call make_geometry_arrays
  if (io_method > 1) call write_hdf5_data(trim(path)//'gyro_data.h5',1)
  !
  ! Deallocate GEO
  call GEO_alloc(0)
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
     if (lskipinit == 0) then

        call read_restart

     else

        ! We will retain the value of h in this case.

        step = 0
        call get_field_explicit

     endif
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

  !------------------------------------------------
  ! Large data dump for vugyro
  !
  call write_profile_vugyro(trim(path)//'profile_vugyro.out',1)
  !------------------------------------------------

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

  if ((lskipinit == 0) .and. (gkeigen_j_set==0)) call gyro_write_master(1)

  !-------------------------------------------------
  ! NEW SIMULATION ONLY:
  !
  ! Write the initial conditions:
  !
  if (restart_method /= 1) then

     if (lskipinit == 0) then
      if (gkeigen_j_set==0) call gyro_write_master(2)
      if (io_method > 1 ) call write_hdf5_timedata(2)
      if (io_method > 1 .and. time_skip_wedge > 0) call write_hdf5_wedge_timedata(2)
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

100 continue

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
