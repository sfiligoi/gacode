!-----------------------------------------------------------
! tgyro_iteration_driver.f90
!
! PURPOSE:
!  Main driver for local solver.
!----------------------------------------------------------

subroutine tgyro_iteration_driver

  use tgyro_globals
  use tgyro_iteration_variables

  implicit none

  n_r   = n_inst+1
  p_max = n_evolve*(n_r-1)

  flux_counter = 0

  allocate(ipiv(p_max))
  allocate(jf(p_max,p_max))
  allocate(jg(p_max,p_max))
  allocate(jfg(p_max,p_max))
  allocate(x_vec(p_max)) 
  allocate(f_vec(p_max))
  allocate(g_vec(p_max))
  allocate(x_vec0(p_max)) 
  allocate(f_vec0(p_max))
  allocate(g_vec0(p_max))
  allocate(b(p_max))

  call tgyro_allocate_globals

  !---------------------------------------
  ! Error monitoring variables
  !
  error_flag = 0
  error_msg = 'INFO: clean exit from TGYRO'
  !
  b_flag(:) = ' ' 
  gyro_exit_status(:)  = 0
  gyro_exit_message(:) = 'N/A'
  !---------------------------------------

  ! Map function from radius/field to p-index.
  p = 0
  do ip=1,n_evolve
     do i=2,n_r
        p = p+1
        pmap(i,ip) = p
     enddo
  enddo

  ! Generate ALL radial profiles.
  call tgyro_init_profiles

  !----------------------------------------------
  ! Choose flux_method based on path information.
  !
  if (lpath == "IFS/") then

     ! IFS-PPPL
     flux_method = 1

     ! Step-length for Jacobian
     dx = loc_dx/r_min

  else if (lpath == "TGLF/") then

     ! TGLF
     flux_method = 2

     ! Step-length for Jacobian
     dx = loc_dx/r_min

  else if (lpath == "QFM/") then

     ! QFM
     flux_method = 4

     ! Step-length for Jacobian
     dx = loc_dx/r_min

  else

     ! GYRO
     flux_method = 3

     ! Step-length for Jacobian
     dx = loc_dx_gyro/r_min

  endif
  !---------------------------------------------

  call tgyro_write_input

  if (tgyro_mode == 2) then
     call tgyro_stab_driver
     return
  endif

  if (loc_restart_flag == 0) then
     ! Create, but do not write to, datafiles.
     call tgyro_write_data(0)
     ! Initialize relaxation parameters to starting value.
     relax(:) = 1.0
  endif
  correct_flag = 0

  p = 0
  do i=2,n_r
     if (loc_ti_feedback_flag == 1) then
        p = p+1
        x_vec(p) = dlntidr(1,i)
     endif
     if (loc_te_feedback_flag == 1) then
        p = p+1
        x_vec(p) = dlntedr(i)
     endif
     if (loc_ne_feedback_flag == 1) then
        p = p+1
        x_vec(p) = dlnnedr(i)
     endif
     if (loc_er_feedback_flag == 1) then
        p = p+1
        x_vec(p) = w0p(i)
     endif
  enddo

  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,*) 'INFO: TGYRO starting iterations'
     close(1)
  endif


  select case (tgyro_iteration_method) 

  case (1) 

     call tgyro_iteration_standard

  case (2,3) 

     call tgyro_iteration_pppl

  case (4) 

     call tgyro_iteration_serial

  case (5) 

     call tgyro_iteration_parallel

  end select

end subroutine tgyro_iteration_driver
