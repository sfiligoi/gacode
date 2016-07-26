!-----------------------------------------------------------
! tgyro_iteration_simplerelax.f90
!
! PURPOSE:
! Trial simple relaxation approach to look for more robust 
! solution method.
!----------------------------------------------------------

subroutine tgyro_iteration_simplerelax

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables

  real :: simpledz

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,'(t2,a)') 'INFO: (TGYRO) Beginning standard iterations'
     close(1)
  endif

  !---------------------------------------------------
  ! ZEROTH ITERATION
  !
  ! One pass to get fluxes.
  !
  ! We *do* want the iteration number to increase in 
  ! the case where we restart without iterating (this
  ! can be used to compare transport models).
  !  
  if (tgyro_relax_iterations == 0 .and. loc_restart_flag == 1) i_tran=i_tran+1
  !
  call tgyro_target_vector(x_vec,g_vec)
  if (loc_restart_flag == 0 .or. tgyro_relax_iterations == 0) then
     ! Need to determine initial fluxes
     gyro_restart_method = 1
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     gyro_restart_method = 2
  else
     ! Initial fluxes already computed
     p = 0
     do i=2,n_r
        if (loc_ti_feedback_flag == 1) then
           p = p+1
           f_vec(p) = eflux_i_tot(i)
        endif
        if (loc_te_feedback_flag == 1) then
           p = p+1
           f_vec(p) = eflux_e_tot(i)
        endif
        if (loc_ne_feedback_flag == 1) then
           p = p+1
           f_vec(p) = pflux_e_tot(i)
        endif
        if (loc_er_feedback_flag == 1) then
           p = p+1
           f_vec(p) = mflux_tot(i)
        endif
        if (loc_he_feedback_flag == 1) then
           p = p+1
           f_vec(p) = pflux_he_tot(i)
        endif
     enddo
     ! GYRO restart data available
     gyro_restart_method = 2
  endif
  res0 = 0.0
  call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)

  if (loc_restart_flag == 0 .or. tgyro_relax_iterations == 0) then
     call tgyro_write_data(1)
  endif
  !----------------------------------------------------

  do i_tran_loop=1,tgyro_relax_iterations

     i_tran = i_tran+1

     ! Initialize gradients

     x_vec0 = x_vec
     g_vec0 = g_vec
     f_vec0 = f_vec
     res0   = res

     !----------------------------------------------
     ! Relax based on (flux-target)/max(target,1)
     ! dz/z = -loc_relax*(Q_tot - Qtarget)/max(Qtarget,Qtot)
     p = 0
     do i=2,n_r
        do p=1,n_evolve
           simpledz = loc_relax*(f_vec(p) - g_vec(p))/&
                max(max(abs(f_vec(p)),abs(g_vec(p))),1.0)
           if (abs(simpledz) > loc_dx_max) then
              simpledz = loc_dx_max*(simpledz/abs(simpledz))
           endif
           x_vec(p) = x_vec(p)*(1.0-simpledz)
        enddo
     enddo

     !----------------------------------------------
     ! Update profiles
     p = 0
     do i=2,n_r
        if (loc_ti_feedback_flag == 1) then
           p = p+1
           dlntidr(therm_vec(:),i) = x_vec(p)
        endif
        if (loc_te_feedback_flag == 1) then
           p = p+1
           dlntedr(i) = x_vec(p)
        endif
        if (loc_ne_feedback_flag == 1) then
           p = p+1
           dlnnedr(i) = x_vec(p)
           ! Set dlnnidr(1,i) according to quasineutrality
           call tgyro_quasigrad(ne(i),dlnnedr(i),ni(:,i),dlnnidr(:,i),zi_vec(:),loc_n_ion)
        endif
        if (loc_er_feedback_flag == 1) then
           p = p+1
           f_rot(i) = x_vec(p)
        endif
        if (loc_he_feedback_flag == 1) then
           p = p+1
           dlnnidr(i_ash,i) = x_vec(p)
        endif
     enddo

     call tgyro_profile_functions 

     !----------------------------------------------------- 
     ! Update fluxes
     call tgyro_target_vector(x_vec,g_vec)
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     !
     ! Compute initial residual
     call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
     call tgyro_write_intermediate(0,res)
     ! Output results
     call tgyro_write_data(1)

  enddo

end subroutine tgyro_iteration_simplerelax
