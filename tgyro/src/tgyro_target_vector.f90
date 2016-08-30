subroutine tgyro_target_vector(x_vec,g_vec)

  use tgyro_globals

  implicit none

  integer :: i,p

  real, intent(in), dimension(p_max) :: x_vec
  real, intent(inout), dimension(p_max) :: g_vec

  call tgyro_profile_set(x_vec,0.0,0)
  call tgyro_profile_functions
  call tgyro_source

  ! Build vector (already converted to GB units in tgyro_source)

  p = 0
  do i=2,n_r
     if (loc_ti_feedback_flag == 1) then
        p = p+1
        g_vec(p) = eflux_i_target(i)
     endif
     if (loc_te_feedback_flag == 1) then
        p = p+1
        g_vec(p) = eflux_e_target(i)
     endif
     if (loc_er_feedback_flag == 1) then
        p = p+1
        g_vec(p) = mflux_target(i)
     endif
     if (loc_ne_feedback_flag == 1) then
        p = p+1
        g_vec(p) = pflux_e_target(i)
     endif
     if (loc_he_feedback_flag == 1) then
        p = p+1
        g_vec(p) = pflux_he_target(i)
     endif
  enddo

end subroutine tgyro_target_vector
