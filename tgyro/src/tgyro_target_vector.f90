subroutine tgyro_target_vector(x_vec,g_vec)

  use tgyro_globals

  implicit none

  integer :: i,p

  real, dimension(p_max) :: x_vec
  real, dimension(p_max) :: g_vec

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
        call tgyro_quasigrad(ne(i),dlnnedr(i),ni(:,i),dlnnidr(:,i),zi_vec(:),loc_n_ion,dlnridr(:,i)) 
     endif
     if (loc_er_feedback_flag == 1) then
        p = p+1
        f_rot(i) = x_vec(p)
     endif
  enddo

  call tgyro_profile_functions
  call tgyro_source

  ! Build vector in erg/cm^2/s, then normalize to bc energy flux

  p = 0
  do i=2,n_r
     if (loc_ti_feedback_flag == 1) then
        p = p+1
        g_vec(p) = eflux_i_target(i)! * q_gb(i)
     endif
     if (loc_te_feedback_flag == 1) then
        p = p+1
        g_vec(p) = eflux_e_target(i)! * q_gb(i)
     endif
     if (loc_ne_feedback_flag == 1) then
        p = p+1
        g_vec(p) = pflux_e_target(i)! * gamma_gb(i) * k * te(i)
     endif
     if (loc_er_feedback_flag == 1) then
        p = p+1
        g_vec(p) = mflux_target(i) !* pi_gb(i) * c_s(i)
     endif
  enddo

end subroutine tgyro_target_vector
