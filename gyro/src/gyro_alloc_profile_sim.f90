!----------------------------------------------------------
! gyro_alloc_profile_sim
!
! PURPOSE:
!  Create and destroy arrays which store profile
!  data on the simulation grid.
!
! NOTES:
!  flag=0: deallocate
!  flag=1: allocate
!
!  Need to know n_spec before this is called.
!----------------------------------------------------------

subroutine gyro_alloc_profile_sim(flag)

  use gyro_globals

  integer, intent(in) :: flag

  include 'mpif.h'

  if (flag == 1 .and. allocated(r)) then
     if (i_proc == 0) then
        print *,'WARNING: already allocated arrays in gyro_alloc_profile_sim'
     endif
     return
  endif
  if (flag == 0 .and. .not.allocated(r)) then
     if (i_proc == 0) then
        print *,'WARNING: cannot deallocate arrays in gyro_alloc_profile_sim'
     endif
     return
  endif

  if (flag == 1) then

     allocate(r(n_x))
     allocate(q(n_x))
     allocate(q_s(n_x))
     allocate(r_s(n_x))
     allocate(rmaj_s(n_x))
     allocate(shat_s(n_x))
     allocate(rhogrid_s(n_x))
     allocate(kappa_s(n_x))
     allocate(s_kappa_s(n_x))
     allocate(delta_s(n_x))
     allocate(zeta_s(n_x))
     allocate(s_delta_s(n_x))
     allocate(s_zeta_s(n_x))
     allocate(z_eff_s(n_x))
     allocate(rhosda_s(n_x))
     allocate(csda_s(n_x))
     allocate(drmaj_s(n_x))
     allocate(dzmag_s(n_x))
     allocate(beta_unit_s(n_x))
     allocate(nu_s(n_spec,n_x))
     allocate(zmag_s(n_x))
     allocate(w0_s(n_x))
     allocate(w0p_s(n_x))
     allocate(gamma_e_s(n_x))
     allocate(gamma_p_s(n_x))
     allocate(mach_s(n_x))
     allocate(dlnpdr_s(n_x))
     allocate(beta_star_s(n_x))
     allocate(omega_eb_s(n_x))

     allocate(den_s(n_spec,n_x))
     allocate(tem_s(n_spec,n_x))
     allocate(b_unit_s(n_x))

     allocate(krho_i(n_n_1,n_x))

     allocate(a_fourier_geo_s(8,0:16,n_x))

     allocate(alpha_s(n_spec,n_x))
     allocate(dlnndr_s(n_spec,n_x))
     allocate(dlntdr_s(n_spec,n_x))
     allocate(pr_s(n_spec,n_x))

     allocate(phase(n_n_1,n_x))
     allocate(angp(n_x))

     allocate(r_e(n_x))
     allocate(dr_eodr(n_x))

     allocate(n(n_n))
     allocate(n_1(n_n_1))
     allocate(z(n_spec))
     allocate(mu(n_spec)) 
     allocate(krho_collect(n_n))
     allocate(c_map(n_spec))

     r=0.
     r_e=0.
     dr_eodr=0.
     q=0.
     r_s=0.
     q_s=0.
     rmaj_s=0.
     krho_i=0.
     shat_s=0.
     ky_rho_star_max=0.
     rhogrid_s=0.
     kappa_s=0.
     s_kappa_s=0.
     delta_s=0.
     zeta_s=0.
     s_delta_s=0.
     s_zeta_s=0.
     z_eff_s=0.
     rhosda_s=0.
     csda_s=0.
     drmaj_s=0.
     dzmag_s=0.
     mach_s=0.
     beta_unit_s=0.
     nu_s=0.
     ve_hat_s=0.
     ve_hat_tor_s=0.
     flow_mom_s=0.
     pow_e_s=0.
     pow_i_s=0.
     pow_ei_s=0.
     powp_s=0.
     zmag_s=0.
     w0_s=0.
     w0p_s=0.
     grad_r0_s=0.
     dlnpdr_s=0.
     beta_star_s=0.
     nui_hh=0.
     nue_hh=0.
     den_s=0.
     tem_s=0.
     b_unit_s=0.
     a_fourier_geo_s=0.
     alpha_s=0.
     dlnndr_s=0.
     dlntdr_s=0.
     pr_s=0.
     pgamma_s=0.
     phi_gamma_s=0.
     diff_to_flow=0.
     phase=0.
     angp=0.
     n=0.
     n_1=0.
     z=0.
     mu=0.
     krho_collect=0.
     c_map=0.


     ! Energy grid
     allocate(energy_max(n_kinetic))
     allocate(energy(n_energy,n_kinetic))
     allocate(w_energy(n_energy,n_kinetic))
     energy_max=0.
     energy=0.
     w_energy=0.

     ! Fourier coefficients
     allocate(cr(n_x,n_x))
     allocate(cri(n_x,n_x))
     cr=0.
     cri=0.

     ! Nondistributed radial operators
     allocate(gyro_trace(n_gk,-m_gyro:m_gyro-i_gyro))
     allocate(w_g0(-m_gyro:m_gyro-i_gyro))
     allocate(w_gd0(-mg_dx:mg_dx-ig_dx))
     allocate(w_d0(-m_dx:m_dx-i_dx))
     allocate(w_d1(-m_dx:m_dx-i_dx))
     allocate(w_d2(-m_dx:m_dx-i_dx))
     allocate(s_d1(-m_dx:m_dx))
     allocate(i_cyc(1-n_x:2*n_x))
     allocate(i_loop(1-n_x:2*n_x))
     allocate(explicit_damp_vec(n_kinetic,n_x))
     gyro_trace=0.
     w_g0=0.
     w_gd0=0.
     w_d0=0.
     w_d1=0.
     w_d2=0.
     s_d1=0.
     i_cyc=0.
     i_loop=0.
     explicit_damp_vec=0.

     if (source_flag == 1) then
        ! Source arrays
        allocate(b_src(n_x,n_lump))
        allocate(m_src(n_lump,n_lump))
        allocate(src_piv(n_lump))

        b_src=0.
        m_src=0.
        src_piv=0.

     endif

     !SEK: There might be a better method here
     if (iohdf5out == 1) then
        allocate(nu_coarse(0:n_theta_plot,n_x))
        allocate(nu_fine(0:n_theta_plot*n_theta_mult,n_x))

        nu_coarse=0.
        nu_fine=0.

     endif

     ! Required in MPI_RECV
     allocate(recv_status(MPI_STATUS_SIZE))

  else 

     deallocate(r)
     deallocate(q)
     deallocate(q_s)
     deallocate(r_s)
     deallocate(rmaj_s)
     deallocate(shat_s)
     deallocate(rhogrid_s)
     deallocate(kappa_s)
     deallocate(s_kappa_s)
     deallocate(delta_s)
     deallocate(zeta_s)
     deallocate(s_delta_s)
     deallocate(s_zeta_s)
     deallocate(z_eff_s)
     deallocate(rhosda_s)
     deallocate(csda_s)
     deallocate(drmaj_s)
     deallocate(dzmag_s)
     deallocate(beta_unit_s)
     deallocate(nu_s)
     deallocate(zmag_s)
     deallocate(w0_s)
     deallocate(w0p_s)
     deallocate(gamma_e_s)
     deallocate(gamma_p_s)
     deallocate(mach_s)
     deallocate(dlnpdr_s)
     deallocate(beta_star_s)
     deallocate(omega_eb_s)

     deallocate(den_s)
     deallocate(tem_s)
     deallocate(b_unit_s)

     deallocate(krho_i)

     deallocate(a_fourier_geo_s)

     deallocate(alpha_s)
     deallocate(dlnndr_s)
     deallocate(dlntdr_s)
     deallocate(pr_s)

     deallocate(phase)
     deallocate(angp)

     deallocate(r_e)
     deallocate(dr_eodr)

     deallocate(n)
     deallocate(n_1)
     deallocate(z)
     deallocate(mu) 
     deallocate(krho_collect)
     deallocate(c_map)

     deallocate(energy_max)
     deallocate(energy)
     deallocate(w_energy)

     deallocate(cr)
     deallocate(cri)

     deallocate(gyro_trace)
     deallocate(w_g0)
     deallocate(w_gd0)
     deallocate(w_d0)
     deallocate(w_d1)
     deallocate(w_d2)
     deallocate(s_d1)
     deallocate(i_cyc)
     deallocate(i_loop)
     deallocate(explicit_damp_vec)

     if (allocated(b_src)) deallocate(b_src)
     if (allocated(m_src)) deallocate(m_src)
     if (allocated(src_piv)) deallocate(src_piv)

     if (iohdf5out == 1) then
        deallocate(nu_coarse)
        deallocate(nu_fine)
     endif

     deallocate(recv_status)

  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_alloc_profile_sim done]'
  endif

end subroutine gyro_alloc_profile_sim
