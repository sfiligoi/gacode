!-----------------------------------------------------------
! tgyro_iteration_driver.f90
!
! PURPOSE:
!  Main driver for local solver.
!
!  Quantities to evolve:
!   Ti: dlntidr
!   Te: dlntedr
!   ne: dlnnedr
!   Er: f_rot = [a/c_s(0)]*r_maj(0)*gamma_p/r_maj
!----------------------------------------------------------

subroutine tgyro_iteration_driver

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables
  use tgyro_ped
  use EXPRO_interface

  implicit none

  integer :: i_ion
  integer :: i_star

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
  allocate(quant(p_max))
  allocate(b(p_max))

  call tgyro_allocate_globals

  !---------------------------------------
  ! Error monitoring variables
  !
  error_flag = 0
  error_msg = 'INFO: clean exit from TGYRO'
  !
  b_flag(:) = ' ' 
  !---------------------------------------

  ! Mapping function from radius/field to p
  if (tgyro_mode /= 2) then
     p = 0
     do i=2,n_r
        ip = 0
        if (loc_ti_feedback_flag == 1) then
           p  = p+1
           ip = ip+1
           pmap(i,ip) = p
           quant(p) = 'ti'
        endif
        if (loc_te_feedback_flag == 1) then
           p  = p+1
           ip = ip+1
           pmap(i,ip) = p
           quant(p) = 'te'
        endif
        if (loc_ne_feedback_flag == 1) then
           p  = p+1
           ip = ip+1
           pmap(i,ip) = p
           quant(p) = 'ne'
        endif
        if (loc_er_feedback_flag == 1) then
           p  = p+1
           ip = ip+1
           pmap(i,ip) = p
           quant(p) = 'er'
        endif
        if (loc_he_feedback_flag == 1) then
           p  = p+1
           ip = ip+1
           pmap(i,ip) = p
           quant(p) = 'he'
        endif
     enddo
  endif

  ! Generate ALL radial profiles.
  call tgyro_init_profiles

  !----------------------------------------------
  ! Choose flux_method based on path information.
  !
  if (tgyro_noturb_flag == 1) then

     flux_method = 5

  else if (lpath(1:3) == "FUN") then

     flux_method = 6
     dx = loc_dx/r_min

  else if (lpath(1:3) == "IFS") then

     ! IFS-PPPL
     flux_method = 1

     ! Step-length for Jacobian
     dx = loc_dx/r_min

  else if (lpath(1:4) == "TGLF") then

     ! TGLF
     flux_method = 2

     ! Step-length for Jacobian
     dx = loc_dx/r_min

  else if (lpath(1:3) == "GLF") then

     ! GLF23
     flux_method = 3

     ! Step-length for Jacobian
     dx = loc_dx/r_min

  else

     ! GYRO
     flux_method = 4

     ! Step-length for Jacobian
     dx = loc_dx_gyro/r_min

  endif
  !---------------------------------------------

  call tgyro_write_input

  ! NOTE: See gyro/src/gyro_globals.f90 for definition of transport_method

  if (tgyro_mode == 2) then
     ! Branch off to stability calculation
     transport_method = 1
     call tgyro_stab_driver
     return
  else
     ! Standard transport calculation
     if (tgyro_gyro_restart_flag == 0) then
        transport_method = 2
     else
        transport_method = 3
     endif
  endif

  if (loc_restart_flag == 0) then
     ! Create, but do not write to, datafiles.
     call tgyro_write_data(0)
     ! Initialize relaxation parameters to starting value.
     relax(:) = 1.0
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,'(t2,a)') 'INFO: (TGYRO) Survived initialization and starting iterations'
     close(1)
  endif

  correct_flag = 0

  p = 0
  do i=2,n_r
     if (loc_ti_feedback_flag == 1) then
        p = p+1
        ! Assume 1 represents the thermal ion temperature
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
        x_vec(p) = f_rot(i)
     endif
     if (loc_he_feedback_flag == 1) then
        p = p+1
        x_vec(p) = dlnnidr(i_ash,i)
     endif
  enddo

  ! Make some resets if we are in test mode
  if (gyrotest_flag == 1) then
     tgyro_iteration_method = 1
     tgyro_relax_iterations = 0
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

  case (6) 

     call tgyro_iteration_simplerelax

  end select

  !--------------------------------------------------------------------------------
  ! Rewrite input.profiles (to input.profiles.new) if flag set
  !
  if (tgyro_write_profiles_flag == 1) then
     call EXPRO_palloc(MPI_COMM_WORLD,'./',1) 
     call EXPRO_pread

     if (i_proc_global == 0) then

        ! Map data inside r < r(n_r)

        if (tgyro_ped_model > 1) then
           ! Map data past r(n_r)
           call tgyro_pedestal_map(dlnnedr(n_r),zn_top,n_top(1),nn_vec(:,2),i_star,EXPRO_ne)
           EXPRO_ne(i_star:n_exp) = EXPRO_ne(i_star:n_exp)*1e-13
           call tgyro_pedestal_map(dlntedr(n_r),zt_top,t_top(1),t_vec(:),i_star,EXPRO_te)
           EXPRO_te(i_star:n_exp) = EXPRO_te(i_star:n_exp)*1e-3
           call tgyro_pedestal_map(dlntidr(1,n_r),zt_top,t_top(1),t_vec(:),i_star,EXPRO_ti(1,:))
           EXPRO_ti(1,i_star:n_exp) = EXPRO_ti(1,i_star:n_exp)*1e-3
           if (loc_n_ion == 1) then
              EXPRO_ni(1,i_star:n_exp) = EXPRO_ne(i_star:n_exp)
           else if (loc_n_ion == 2) then
              EXPRO_ni(1,i_star:n_exp) = EXPRO_ne(i_star:n_exp)-zi_vec(2)*EXPRO_ni(2,i_star:n_exp)
              if (therm_flag(2) == 1) EXPRO_ti(2,i_star:n_exp) = EXPRO_ti(1,i_star:n_exp)
           endif
        endif

        call tgyro_expro_map(r,dlnnedr,n_r,100*EXPRO_rmin,EXPRO_ne,EXPRO_n_exp)
        call tgyro_expro_map(r,dlntedr,n_r,100*EXPRO_rmin,EXPRO_te,EXPRO_n_exp)
        do i_ion=1,loc_n_ion
           if (therm_flag(i_ion) == 1) then
              call tgyro_expro_map(r,dlnnidr(i_ion,:),n_r,100*EXPRO_rmin,EXPRO_ni(i_ion,:),EXPRO_n_exp)
              call tgyro_expro_map(r,dlntidr(i_ion,:),n_r,100*EXPRO_rmin,EXPRO_ti(i_ion,:),EXPRO_n_exp)
           endif
        enddo

        ! Write data to file
        call EXPRO_write_original(&
             1,'input.profiles',&
             2,'input.profiles.new',&
             'Profiles modified by TGYRO')
        call EXPRO_compute_derived
        call EXPRO_write_derived(1,'input.profiles.extra')

     endif

     call EXPRO_palloc(MPI_COMM_WORLD,'./',0) 

  endif
  !--------------------------------------------------------------------------------

end subroutine tgyro_iteration_driver
