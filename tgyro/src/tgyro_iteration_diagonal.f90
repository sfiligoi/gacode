!-----------------------------------------------------------
! tgyro_iteration_diagonal.f90
!
! PURPOSE:
! Diagonal relaxation approach.
!----------------------------------------------------------

subroutine tgyro_iteration_diagonal

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables

  real :: dz
  real, dimension(:), allocatable :: fn0,fn
  logical :: evomain,evoaux

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
     call tgyro_flux_set(f_vec)
     ! GYRO restart data available
     gyro_restart_method = 2
  endif
  res0 = 0.0
  call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)

  if (loc_restart_flag == 0 .or. tgyro_relax_iterations == 0) then
     call tgyro_write_data(1)
  endif
  !----------------------------------------------------

  allocate(fn0(p_max))
  allocate(fn(p_max))

  evomain = (loc_ti_feedback_flag == 1) .or. (loc_te_feedback_flag == 1)
  evoaux  = (evo_e(0) == 1) .or. (loc_er_feedback_flag == 1)

  do i_tran_loop=1,tgyro_relax_iterations

     i_tran = i_tran+1

     ! Initialize gradients

     x_vec0 = x_vec
     g_vec0 = g_vec
     f_vec0 = f_vec
     res0   = res

     if (i_tran_loop == 1) then
        relax = loc_relax
     endif


     if (evomain) then

        fn0 = g_vec-f_vec
        call tgyro_flux_vector(x_vec0,f_vec,loc_dx,-1)
        fn = g_vec-f_vec

        p = 0
        do i=2,n_r
           if (loc_ti_feedback_flag == 1) then
              p = p+1
              x_vec(p) = x_vec0(p)-loc_relax*fn0(p)/(fn(p)-fn0(p))*loc_dx
           endif
           if (loc_te_feedback_flag == 1) then
              p = p+1
              x_vec(p) = x_vec0(p)-loc_relax*fn0(p)/(fn(p)-fn0(p))*loc_dx
           endif
           if (loc_er_feedback_flag == 1) then
              p = p+1
           endif
           if (evo_e(0) == 1) then 
              p = p+1
           endif
        enddo

        call tgyro_target_vector(x_vec,g_vec)
        call tgyro_flux_vector(x_vec,f_vec,0.0,0)

     endif

     if (evoaux .and. mod(i_tran_loop,2)==1) then 

        fn0 = g_vec-f_vec
        call tgyro_flux_vector(x_vec,f_vec,loc_dx,-2)
        fn = g_vec-f_vec

        p = 0
        do i=2,n_r
           if (loc_ti_feedback_flag == 1) then
              p = p+1
           endif
           if (loc_te_feedback_flag == 1) then
              p = p+1
           endif
           if (loc_er_feedback_flag == 1) then
              p = p+1
              x_vec(p) = x_vec0(p)-loc_relax*fn0(p)/(fn(p)-fn0(p))*loc_dx
           endif
           if (evo_e(0) == 1) then 
              p = p+1
              x_vec(p) = x_vec0(p)-loc_relax*fn0(p)/(fn(p)-fn0(p))*loc_dx
           endif
        enddo

        call tgyro_target_vector(x_vec,g_vec)
        call tgyro_flux_vector(x_vec,f_vec,0.0,0)

     endif

     ! Compute residual
     call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)

     ! Output results
     call tgyro_write_intermediate(0,res)
     call tgyro_write_data(1)

  enddo

  deallocate(fn0,fn)

end subroutine tgyro_iteration_diagonal
