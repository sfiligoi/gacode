subroutine tgyro_flux_vector_dense(x_vec,f_vec)

  use tgyro_globals

  implicit none

  integer :: i
  integer :: p

  real, intent(in), dimension(p_max) :: x_vec
  real, intent(inout), dimension(p_max) :: f_vec


  p = 0
  do i=2,n_r
     if (loc_ti_feedback_flag == 1) then
        p = p+1
        dlntidr(1,i) = x_vec(p)
        if (loc_n_ion == 2) dlntidr(2,i) = dlntidr(1,i)
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
  enddo

  ! Make profile consistent with new gradients
  call tgyro_profile_functions
  call tgyro_flux
  call tgyro_comm_sync

  ! LP: make vector in erg/cm^2/sec, then normalize to pivot ion energy flux
  p = 0
  do i=2,n_r

     if (loc_ti_feedback_flag == 1) then
        p = p+1
        f_vec(p) = eflux_i_tot(i) ! * q_gb(i)
     endif

     if (loc_te_feedback_flag == 1) then
        p = p+1
        f_vec(p) = eflux_e_tot(i) ! * q_gb(i)
     endif

     if (loc_ne_feedback_flag == 1) then
        p = p+1
        f_vec(p) = pflux_e_tot(i) ! * gamma_gb(i) * k * te(i)
     endif

     if (loc_er_feedback_flag == 1) then
        p = p+1
        f_vec(p) = mflux_tot(i) ! * pi_gb(i) * c_s(i)
     endif

  enddo

end subroutine tgyro_flux_vector_dense
