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

  implicit none

  allocate(fn0(p_max))
  allocate(fn(p_max))

  quasifix = 0
  
  do i_tran_loop=1,tgyro_relax_iterations

     i_tran = i_tran+1

     ! Update Ti,Te

     if (loc_ti_feedback_flag == 1) then
        fn0 = g_vec-f_vec
        call tgyro_flux_vector(x_vec,f_vec,loc_dx,1)
        fn  = g_vec-f_vec

        call get_dx

        do p=1,p_max
           if (mask(p,1) == 1) x_vec(p) = x_vec(p)-b(p)
        enddo

        call tgyro_target_vector(x_vec,g_vec)
        call tgyro_flux_vector(x_vec,f_vec,0.0,0)
        
        ! Compute residual
        res0 = res
        call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
        do p=1,p_max
           if (mask(p,1) == 1 .and. res0(p) < res(p)) x_vec(p) = x_vec(p)+b(p)*0.5
        enddo

     endif

     if (loc_te_feedback_flag == 1) then
        fn0 = g_vec-f_vec
        call tgyro_flux_vector(x_vec,f_vec,loc_dx,2)
        fn  = g_vec-f_vec

        call get_dx

        do p=1,p_max
           if (mask(p,2) == 1) x_vec(p) = x_vec(p)-b(p)
        enddo

        call tgyro_target_vector(x_vec,g_vec)
        call tgyro_flux_vector(x_vec,f_vec,0.0,0)
        ! Compute residual
        res0 = res
        call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
        do p=1,p_max
           if (mask(p,2) == 1 .and. res0(p) < res(p)) x_vec(p) = x_vec(p)+b(p)*0.5
        enddo
     endif

     if (evo_e(0) == 1 .or. loc_er_feedback_flag == 1) then
        fn0 = g_vec-f_vec
        call tgyro_flux_vector(x_vec,f_vec,loc_dx,-1)
        fn  = g_vec-f_vec

        call get_dx

        do p=1,p_max
           if (loc_er_feedback_flag ==1 .and. mask(p,3) == 1) x_vec(p) = x_vec(p)-b(p)*0.5
           if (evo_e(0) == 1 .and. mask(p,4) == 1) x_vec(p) = x_vec(p)-b(p)*0.5
        enddo

        call tgyro_target_vector(x_vec,g_vec)
        call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     endif

     ! Compute residual
     res0 = res
     call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)

     ! Output results
     call tgyro_write_intermediate(0,res)
     call tgyro_write_data(1)

  enddo

  deallocate(fn0,fn)

end subroutine tgyro_iteration_diagonal

subroutine get_dx

  use tgyro_globals
  use tgyro_iteration_variables

  implicit none

  do p=1,p_max
     if (abs(fn(p)-fn0(p)) > 1e-5) then
        b(p) = loc_relax*fn0(p)/(fn(p)-fn0(p))*loc_dx
        if (abs(b(p)) > loc_dx_max) then
           b(p) = loc_dx_max*b(p)/abs(b(p))
        endif
     else
        b(p) = 0.0
     endif
  enddo

end subroutine get_dx
