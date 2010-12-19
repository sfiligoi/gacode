!----------------------------------------------------------
! gyro_alloc_orbit
!
! PURPOSE:
!  Create and destroy orbit arrays.
!
! NOTES:
!  flag=0: deallocate
!  flag=1: allocate
!
!  Need to know n_lambda,n_stack before this is called.
!
!  See: gyro_banana_operators, make_geometry_arrays.
!----------------------------------------------------------

subroutine gyro_alloc_orbit(flag)

  use gyro_globals
  use gyro_collision_private

  integer, intent(in) :: flag

  if (flag == 1 .and. allocated(theta_t)) then
     if (i_proc == 0) then
        print *,'WARNING: already allocated arrays in gyro_alloc_orbit'
     endif
     return
  endif
  if (flag == 0 .and. .not.allocated(theta_t)) then
     if (i_proc == 0) then
        print *,'WARNING: cannot deallocate arrays in gyro_alloc_orbit'
     endif
     return
  endif

  if (flag == 1) then

     ! Theta/tau grid functions

     allocate(theta_t(n_x,n_lambda,n_stack))
     allocate(tau(n_x,n_lambda,n_stack))
     allocate(m_cyc(n_class,-n_stack+1:2*n_stack,2))
     allocate(p_cyc(n_class,n_x,-n_stack+1:2*n_stack,2))
     allocate(m_phys(n_class,n_stack))
     allocate(p_phys(n_class,n_stack))
     allocate(omega(n_x,n_lambda))
     allocate(m_map(n_class,n_theta(2),2))
     allocate(theta_int(n_theta_int))
     allocate(theta_plot(n_theta_plot))
     allocate(theta_r0_plot(field_r0_grid))

     ! Geometry functions
     allocate(b0_t(n_x,n_lambda,n_stack))
     allocate(g_theta_t(n_x,n_lambda,n_stack))
     allocate(grad_r_t(n_x,n_lambda,n_stack))
     allocate(qrat_t(n_x,n_lambda,n_stack))
     allocate(cos_t(n_x,n_lambda,n_stack))
     allocate(cos_p_t(n_x,n_lambda,n_stack))   
     allocate(captheta_t(n_x,n_lambda,n_stack))   
     allocate(sin_t(n_x,n_lambda,n_stack))   
     allocate(bt_t(n_x,n_lambda,n_stack))
     allocate(bp_t(n_x,n_lambda,n_stack))
     allocate(bigr_t(n_x,n_lambda,n_stack))
     allocate(b0_plot(n_x,n_theta_plot))
     allocate(g_theta_plot(n_x,n_theta_plot))   
     allocate(usin_t(n_x,n_lambda,n_stack))   
     allocate(ucos_t(n_x,n_lambda,n_stack))   

     ! Blending arrays
     allocate(c_fluxave(n_blend,n_x))
     allocate(ff_mm(n_blend,n_blend,n_x,2))
     allocate(ff_mm_piv(n_blend,n_x,2))
     allocate(ff2_mm(2*n_blend,2*n_blend,n_x))
     allocate(ff2_mm_piv(2*n_blend,n_x))
     allocate(blend_plot(n_blend,n_theta_plot,n_x))
     allocate(blend_prime_plot(n_blend,n_theta_plot,n_x))
     allocate(blend_r0_plot(n_blend,field_r0_grid))

     if (collision_flag == 1) then
        ! Collision arrays
        allocate(nu_total(n_x,n_energy,indx_e))
        allocate(xi(n_x,n_lambda,n_stack))
        if (collision_method > 2) then
           ! for ebelli collisions
           allocate(nu_coll_d(n_x,n_energy,n_kinetic,n_kinetic))
           allocate(rs_coll_const(n_x,n_kinetic,n_kinetic))
           allocate(rs_nunu_const(n_x,n_kinetic,n_kinetic,n_kinetic))
           allocate(indx_coll(n_kinetic,n_kinetic))
        endif
     endif

     if (n_field > 1) then
        allocate(coll_vel(n_x,n_blend,n_blend))
     endif

     if(n_field == 3) then
        allocate(coll_vel_perp1(n_x,n_blend,n_blend))
        allocate(coll_vel_perp2(n_x,n_blend,n_blend))
     endif

  else 

     deallocate(theta_t)
     deallocate(tau)
     deallocate(m_cyc)
     deallocate(p_cyc)
     deallocate(m_phys)
     deallocate(p_phys)
     deallocate(omega)
     deallocate(m_map)
     deallocate(theta_int)
     deallocate(theta_plot)
     deallocate(theta_r0_plot)

     deallocate(b0_t)
     deallocate(g_theta_t)
     deallocate(grad_r_t)
     deallocate(qrat_t)
     deallocate(cos_t)
     deallocate(cos_p_t)   
     deallocate(captheta_t)   
     deallocate(sin_t)   
     deallocate(bt_t)
     deallocate(bp_t)
     deallocate(bigr_t)
     deallocate(b0_plot)
     deallocate(g_theta_plot)   
     deallocate(usin_t)
     deallocate(ucos_t)

     deallocate(c_fluxave)
     deallocate(ff_mm)
     deallocate(ff_mm_piv)
     deallocate(ff2_mm)
     deallocate(ff2_mm_piv)
     deallocate(blend_plot)
     deallocate(blend_prime_plot)
     deallocate(blend_r0_plot)

     if (allocated(nu_total)) deallocate(nu_total)
     if (allocated(xi)) deallocate(xi)
     if (allocated(nu_coll_d))     deallocate(nu_coll_d)
     if (allocated(rs_coll_const))     deallocate(rs_coll_const)
     if (allocated(rs_nunu_const))     deallocate(rs_nunu_const)
     if (allocated(indx_coll))         deallocate(indx_coll)
     if (allocated(coll_vel))       deallocate(coll_vel)
     if (allocated(coll_vel_perp1)) deallocate(coll_vel_perp1)
     if (allocated(coll_vel_perp2)) deallocate(coll_vel_perp2)

  endif

end subroutine gyro_alloc_orbit
