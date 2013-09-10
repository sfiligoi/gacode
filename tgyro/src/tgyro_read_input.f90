!----------------------------------------------------------------
! tgyro_read_input.f90
!
! PURPOSE:
!  Read input from input.dat including paths for each instance
!  of GYRO.
!
! NOTES:
!  - transp_fix_density_flag signals holding density profile 
!    fixed
!----------------------------------------------------------------

subroutine tgyro_read_input

  use mpi
  use tgyro_globals

  !-----------------------------------------------
  implicit none
  !
  integer :: ind
  integer :: ioerr
  integer :: i
  logical :: fileex
  character(len=20) :: filename='input.tgyro.gen'
  character(len=80) :: ipath
  !-----------------------------------------------


  if (i_proc_global == 0) then

     inquire(file=filename,exist=fileex)
     if (fileex) then
        open(unit=1,file=filename,status='old',iostat=ioerr)
     else
        call tgyro_catch_error('ERROR: (TGYRO) cannot open '//filename)
     endif  ! fileex

  endif ! i_proc_global

  !-------------------------------------------------------
  ! Read/BCAST input parameters
  !
  call tgyro_readbc_int(tgyro_mode)
  call tgyro_readbc_int(tgyro_relax_iterations)
  call tgyro_readbc_real(loc_alpha_elec) 
  call tgyro_readbc_real(loc_nu_scale) 
  call tgyro_readbc_real(loc_dx)  
  call tgyro_readbc_real(loc_dx_gyro) 
  call tgyro_readbc_real(loc_dx_max) 
  call tgyro_readbc_real(loc_relax) 
  call tgyro_readbc_int(loc_lock_profile_flag) 
  call tgyro_readbc_int(loc_evolve_grad_only_flag) 
  call tgyro_readbc_int(loc_restart_flag) 
  call tgyro_readbc_int(loc_scenario) 
  call tgyro_readbc_int(loc_quasineutral_flag) 
  call tgyro_readbc_int(loc_neo_method) 
  call tgyro_readbc_int(loc_n_ion) 
  call tgyro_readbc_real(zi_vec(1)) 
  call tgyro_readbc_real(zi_vec(2)) 
  call tgyro_readbc_real(zi_vec(3)) 
  call tgyro_readbc_real(zi_vec(4)) 
  call tgyro_readbc_real(zi_vec(5)) 
  call tgyro_readbc_real(mi_vec(1)) 
  call tgyro_readbc_real(mi_vec(2)) 
  call tgyro_readbc_real(mi_vec(3)) 
  call tgyro_readbc_real(mi_vec(4)) 
  call tgyro_readbc_real(mi_vec(5)) 
  call tgyro_readbc_int(therm_flag(1)) 
  call tgyro_readbc_int(therm_flag(2)) 
  call tgyro_readbc_int(therm_flag(3)) 
  call tgyro_readbc_int(therm_flag(4)) 
  call tgyro_readbc_int(therm_flag(5)) 
  call tgyro_readbc_real(loc_betae_scale) 
  call tgyro_readbc_int(loc_chang_hinton) 
  call tgyro_readbc_real(loc_me_multiplier) 
  call tgyro_readbc_int(loc_sawtooth_model) 
  call tgyro_readbc_int(loc_bc_offset)
  call tgyro_readbc_int(tgyro_tglf_revision)
  call tgyro_readbc_int(tgyro_tglf_dump_flag)
  call tgyro_readbc_int(loc_ti_feedback_flag)
  call tgyro_readbc_int(loc_te_feedback_flag)
  call tgyro_readbc_int(loc_ne_feedback_flag)
  call tgyro_readbc_int(loc_er_feedback_flag)
  call tgyro_readbc_int(loc_zeff_flag)
  call tgyro_readbc_int(loc_pflux_method)
  call tgyro_readbc_int(loc_residual_method)
  call tgyro_readbc_int(loc_num_equil_flag)
  call tgyro_readbc_int(tgyro_neo_gv_flag)
  call tgyro_readbc_int(tglf_q_low_flag)
  call tgyro_readbc_int(tgyro_global_newton_flag)
  call tgyro_readbc_int(tgyro_backtrack_method)
  call tgyro_readbc_int(tgyro_iteration_method)
  call tgyro_readbc_real(lm_boost)
  call tgyro_readbc_real(lm_drop)
  call tgyro_readbc_int(tgyro_rotation_flag)
  call tgyro_readbc_int(tgyro_rotation_theory_method)
  call tgyro_readbc_int(tgyro_stab_nsearch)
  call tgyro_readbc_int(tgyro_stab_nky)
  call tgyro_readbc_real(tgyro_stab_kymin)
  call tgyro_readbc_real(tgyro_stab_deltaky)
  call tgyro_readbc_real(tgyro_rmin)
  call tgyro_readbc_real(tgyro_rmax)
  call tgyro_readbc_int(tgyro_global_radii)
  call tgyro_readbc_int(tgyro_expwd_flag)
  call tgyro_readbc_real(tgyro_input_den_scale)
  call tgyro_readbc_real(tgyro_input_te_scale)
  call tgyro_readbc_real(tgyro_input_ti_scale)
  call tgyro_readbc_real(tgyro_input_w0_scale)
  call tgyro_readbc_real(tgyro_input_paux_scale)
  call tgyro_readbc_int(tgyro_er_bc)
  call tgyro_readbc_int(tgyro_noturb_flag)
  call tgyro_readbc_int(tgyro_use_rho)
  call tgyro_readbc_int(tgyro_dt_method)
  call tgyro_readbc_int(tgyro_gyro_restart_flag)
  call tgyro_readbc_int(tgyro_fix_concentration_flag)
  call tgyro_readbc_int(tgyro_write_profiles_flag)
  ! ** END input read; ADD NEW PARAMETERS ABOVE HERE!!
  call tgyro_readbc_int(n_inst)
  !-------------------------------------------------------

  allocate(paths(n_inst))
  allocate(procs(n_inst))

  if (i_proc_global == 0) then

     ! Read code path/CPU count

     do i=1,n_inst

        read(1,*) ipath,procs(i)
        ind = index(ipath,' ')

        ! Append '/' to path name for use later

        paths(i) = ipath(1:ind-1)//'/'

     enddo ! i

     close(1)

     open(unit=1,file=trim(runfile),position='append')

     write(1,*)'-----------------------------------------------------------------'
     write(1,*) 'CPU Partition list'
     write(1,*)
     write(1,*) '----------   ----------------  --------'
     write(1,*) ' Instance     Directory Name     CPUs  '
     write(1,*) '----------   ----------------  --------'

     do i=1,n_inst
        write(1,'(t3,i2,t19,a,t32,i5)') i,trim(paths(i)),procs(i)
     enddo ! i

     write(1,*) '----------   ----------------  --------'
     write(1,'(t22,a,t32,i5)') 'Total:',sum(procs)
     write(1,*)

     close(1)

  endif ! i_proc_global

  !-------------------------------------------------------------------
  ! Check for gyrotest flag by reading a single integer.
  !
  ! 0: not a test
  ! 1: run in single-processor test mode
  !
  if (i_proc_global == 0) then
     open(unit=1,file='gyrotest_flag',status='old')
     read(1,*) gyrotest_flag
     close(1)
  endif
  call MPI_BCAST(gyrotest_flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (gyrotest_flag == 1) then
     ! No iterations; only one call to GYRO 
     tgyro_relax_iterations = 0
     ! One CPU per instance
     procs(:) = 1
  endif
  !------------------------------------------------------------------

end subroutine tgyro_read_input
