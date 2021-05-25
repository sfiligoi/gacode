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

  character(len=30) :: final_msg

  character(8)  :: sdate
  character(10) :: stime
  character(len=64)  :: platform
  integer(KIND=8) :: start_time,aftermpi_time,beforetotal_time,exit_time
  integer(KIND=8) :: count_rate,count_max
  real :: mpi_dt,init_dt,exit_dt
  integer :: statusfd

  ! initialize tiny float
  small = tiny(0.0)
  
  ! the time_lib relies on MPI being initalized, so need to use lower level functions for this
  call system_clock(start_time,count_rate,count_max)

  i_time = 0

  ! Need to initialize the info runfile very early
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_info,status='replace')
  endif

  ! 1. MPI setup
  call cgyro_mpi_grid
  if (error_status > 0) goto 100

  call system_clock(aftermpi_time,count_rate,count_max)
  if (aftermpi_time > start_time) then
     mpi_dt = (aftermpi_time-start_time)/real(count_rate)
  else
     mpi_dt = (aftermpi_time-start_time+count_max)/real(count_rate)
  endif

  ! 2. Profile setup
  call cgyro_make_profiles
  if (error_status > 0) goto 100

  ! 3. Parameter consistency checks
  call cgyro_check
  if (error_status > 0) goto 100

  ! 4. Array initialization and construction
  !    NOTE: On exit, field_old = field 

  call cgyro_init_manager
  if (error_status /=0 ) then
     ! something went terribly wrong, hard abort, as things may be
     ! in weird state
     call abort
  endif
  if (test_flag == 1) return

  !---------------------------------------------------------------------------
  !
  ! Time-stepping
  n_time = nint(max_time/delta_t)

  if (restart_flag == 1) then
     io_control = 3*(1-silent_flag)
  else
     io_control = 1*(1-silent_flag)
  endif

  call timer_lib_in('io_init')
  call cgyro_write_timedata
  call timer_lib_out('io_init')
  call write_timers(trim(path)//runfile_timers)

  io_control = 2*(1-silent_flag)

  call system_clock(beforetotal_time,count_rate,count_max)
  if (beforetotal_time > start_time) then
    init_dt = (beforetotal_time-start_time)/real(count_rate)
  else
    init_dt = (beforetotal_time-start_time+count_max)/real(count_rate)
  endif

  if (i_proc == 0) then
    call date_and_time(sdate,stime);
    call get_environment_variable('GACODE_PLATFORM',platform)
    open(NEWUNIT=statusfd,FILE=trim(path)//runfile_startups,action='write',status='unknown',position='append')
    write(statusfd,'(14(a),f7.3,a,f7.3,a)') &
         sdate(1:4),'/',sdate(5:6),'/',sdate(7:8),' ', &
         stime(1:2),':',stime(3:4),':',stime(5:6),' ', &
         trim(platform),' [  INITIALIZATION] Time =',init_dt,' (mpi init: ', mpi_dt, ')'
    close(statusfd)
  endif

  ! Initialize adaptive time-stepping parameter
  delta_t_gk = delta_t
  
  ! GPU versions of step_gk and coll work on the following in the GPU memory
  call timer_lib_in('str_mem')
!$acc update device(field,psi,cap_h_c,chi,h_x,source)
  call timer_lib_out('str_mem')

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

     ! field will not be modified in GPU memory for the rest fo the loop
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

100 continue


  call timer_lib_in('coll_mem')
!$acc update host(field,psi,cap_h_c,chi,h_x,source,rhs(:,:,1))
  call timer_lib_out('coll_mem')

  ! Manage exit message

  if (error_status == 0) then
     if (nonlinear_flag == 1) then
        final_msg = 'Normal'
     else
        if (signal == 1) then
           final_msg = 'Linear converged'
        else
           final_msg = 'Linear terminated at max time'
        endif
     endif
     if (silent_flag == 0 .and. i_proc == 0) then
        open(unit=io,file=trim(path)//runfile_info,status='old',position='append')
        write(io,'(a)') 'EXIT: (CGYRO) '//trim(final_msg)
        close(io)
     endif
  endif

  call cgyro_cleanup

  call system_clock(exit_time,count_rate,count_max)
  if (exit_time.gt.start_time) then
    exit_dt = (exit_time-start_time)/real(count_rate)
  else
    exit_dt = (exit_time-start_time+count_max)/real(count_rate)
  endif

end subroutine cgyro_kernel
